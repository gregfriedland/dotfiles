# Reserve K8s Pod

Reserve a compute pod on the west1 GKE cluster for the specified number of hours.

## Usage
```
/reserve-k8s-pod [hours] [--type TYPE] [--image IMAGE] [--name NAME]
```

**Parameters:**
- `hours` - Duration in hours (default: 5)
- `--type TYPE` - Pod type: `gpu` (default) or `cpu`
- `--image IMAGE` - Docker image to use (optional)
- `--name NAME` - Custom pod name (optional, auto-generated if not specified)

## Pod Types

| Type | Node | Resources | Use Case |
|------|------|-----------|----------|
| `gpu` | g2-standard-16 | 16 vCPU, 64GB RAM, 1x NVIDIA L4 | Deep learning, GPU compute |
| `rtx6000` | g4-standard-48 | 48 vCPU, 150GB RAM, 1x RTX Pro 6000 (96GB) | Blackwell GPU, large models |
| `tpu-v5e` | ct5lp-highcpu-1t | 24 vCPU, 48GB RAM, 1x TPU v5e | JAX/TPU workloads, ML training/inference |
| `cpu` | t2d-standard-60 | 60 vCPU, 240GB RAM | CPU-intensive work (Rosetta, simulations) |

**Default image:** `us-central1-docker.pkg.dev/gke-test-421317/flyte/internal-base:linux-cuda-66cad2902ac7512a4ee08efbd03a2e3944441532`

## Common Images

| Image | Description |
|-------|-------------|
| `boltz/boltz` | Boltz-2 protein structure prediction |
| `us-central1-docker.pkg.dev/gke-test-421317/flyte/internal-gnina:linux-cuda-66cad2902ac7512a4ee08efbd03a2e3944441532` | GNINA docking (has gnina pre-installed) |
| `us-central1-docker.pkg.dev/gke-test-421317/flyte/gnina-build-base:latest` | GNINA build/test environment |
| `nvcr.io/nvidia/pytorch:24.01-py3` | NVIDIA PyTorch with CUDA |
| (default) | Internal base image with CUDA |

## Examples

```bash
# Default GPU pod for 5 hours
/reserve-k8s-pod

# GPU pod for 2 hours
/reserve-k8s-pod 2 --type gpu

# Large CPU pod for 4 hours (good for Rosetta, PyRosetta)
/reserve-k8s-pod 4 --type cpu

# CPU pod with custom name
/reserve-k8s-pod 2 --type cpu --name rosetta-sampling

# GPU pod with Boltz image
/reserve-k8s-pod 4 --type gpu --image boltz/boltz

# GPU pod with GNINA build environment
/reserve-k8s-pod --image us-central1-docker.pkg.dev/gke-test-421317/flyte/gnina-build-base:latest

# RTX Pro 6000 Blackwell pod (96GB VRAM)
/reserve-k8s-pod 5 --type rtx6000

# RTX 6000 with gnina-build-base for Blackwell builds
/reserve-k8s-pod 5 --type rtx6000 --image us-central1-docker.pkg.dev/gke-test-421317/flyte/gnina-build-base:latest

# On-demand TPU v5e pod for 6 hours
/reserve-k8s-pod 6 --type tpu-v5e
```

## Instructions

When this command is invoked:

### 0. Check for Existing Pods

**ALWAYS check for existing running pods before creating a new one.** Run:

```bash
# Check for existing GPU pods
kubectl --context=west1 get pods -n development -o wide | grep -E '(gpu-|gnina-)'

# Check for existing TPU pods
kubectl --context=west1 get pods -n development -o wide | grep tpu-

# Check for existing CPU pods
kubectl --context=west1 get pods -n development -o wide | grep cpu-
```

If a suitable pod already exists (same type, Running status), **run `top` on it** to check if it's available (low CPU usage indicates it's idle). If the pod is available, use it instead of creating a new one. This avoids wasting resources by creating duplicate pods.

### 1. Generate Pod Name

Generate a unique pod name using this format:
```
{type}-{image_short}-{unique_id}
```

Where:
- `{type}` = `gpu` or `cpu`
- `{image_short}` = short name derived from image (e.g., `boltz`, `gnina`, `base`)
- `{unique_id}` = 6 character random alphanumeric string

**Image short name mapping:**
| Image | Short Name |
|-------|------------|
| `boltz/boltz` | `boltz` |
| `us-central1-docker.pkg.dev/gke-test-421317/flyte/internal-gnina:*` | `gnina` |
| `us-central1-docker.pkg.dev/gke-test-421317/flyte/gnina-build-base:*` | `gnina-build` |
| `nvcr.io/nvidia/pytorch:*` | `pytorch` |
| (default) | `base` |

Example generated names:
- `gpu-boltz-a3f8k2`
- `cpu-base-x9m4p7`

If `--name` is provided, use that name directly instead.

### 2. Create GPU Pod (--type gpu)

```bash
kubectl --context=west1 run POD_NAME --image=IMAGE_HERE --namespace=development --restart=Never --overrides='{"spec": {"nodeSelector": {"cloud.google.com/gke-accelerator": "nvidia-l4"}, "affinity": {"nodeAffinity": {"requiredDuringSchedulingIgnoredDuringExecution": {"nodeSelectorTerms": [{"matchExpressions": [{"key": "flyte.org/node-role", "operator": "In", "values": ["worker"]}, {"key": "cloud.google.com/gke-spot", "operator": "DoesNotExist"}]}]}}}, "tolerations": [{"key": "flyte.org/node-role", "operator": "Equal", "value": "worker", "effect": "NoSchedule"}, {"key": "nvidia.com/gpu", "operator": "Equal", "value": "present", "effect": "NoSchedule"}, {"key": "nvidia.com/gpu", "operator": "Exists", "effect": "NoSchedule"}], "containers": [{"name": "gpu-container", "image": "IMAGE_HERE", "command": ["sleep", "SECONDS"], "env": [{"name": "LD_LIBRARY_PATH", "value": "/usr/local/nvidia/lib64"}, {"name": "PATH", "value": "/usr/local/nvidia/bin:/usr/local/cuda/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin"}], "resources": {"requests": {"cpu": "15", "memory": "57600M", "ephemeral-storage": "30G", "nvidia.com/gpu": "1"}, "limits": {"cpu": "16", "memory": "57600M", "ephemeral-storage": "30G", "nvidia.com/gpu": "1"}}, "securityContext": {"privileged": true}}]}}' -- sleep SECONDS
```

### 2b. Create RTX Pro 6000 Blackwell Pod (--type rtx6000)

**IMPORTANT**: RTX 6000 pods use the `central1` context, not `west1`.

```bash
kubectl --context=central1 run POD_NAME --image=IMAGE_HERE --namespace=development --restart=Never --overrides='{"spec": {"nodeSelector": {"cloud.google.com/gke-accelerator": "nvidia-rtx-pro-6000"}, "affinity": {"nodeAffinity": {"requiredDuringSchedulingIgnoredDuringExecution": {"nodeSelectorTerms": [{"matchExpressions": [{"key": "flyte.org/node-role", "operator": "In", "values": ["worker"]}, {"key": "cloud.google.com/gke-spot", "operator": "DoesNotExist"}]}]}}}, "tolerations": [{"key": "flyte.org/node-role", "operator": "Equal", "value": "worker", "effect": "NoSchedule"}, {"key": "nvidia.com/gpu", "operator": "Equal", "value": "present", "effect": "NoSchedule"}, {"key": "nvidia.com/gpu", "operator": "Exists", "effect": "NoSchedule"}], "containers": [{"name": "gpu-container", "image": "IMAGE_HERE", "command": ["sleep", "SECONDS"], "env": [{"name": "LD_LIBRARY_PATH", "value": "/usr/local/nvidia/lib64:/opt/libtorch/lib:/usr/local/lib"}, {"name": "PATH", "value": "/usr/local/nvidia/bin:/usr/local/cuda/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin"}], "resources": {"requests": {"cpu": "46", "memory": "150G", "ephemeral-storage": "35G", "nvidia.com/gpu": "1"}, "limits": {"cpu": "48", "memory": "160G", "ephemeral-storage": "40G", "nvidia.com/gpu": "1"}}, "securityContext": {"privileged": true}}]}}' -- sleep SECONDS
```

**RTX 6000 Pod Notes:**
- **Cluster**: central1 (not west1)
- **GPU**: NVIDIA RTX Pro 6000 Blackwell Server Edition (96GB VRAM, sm_120)
- **Node**: g4-standard-48 (48 vCPU, 192GB RAM)
- **Ephemeral Storage**: 35-40GB
- **CUDA**: Requires CUDA 12.8+ for Blackwell (sm_120) support
- **Use case**: Large model inference, Blackwell architecture testing, gnina builds with sm_120

### 2c. Create TPU v5e Pod (--type tpu-v5e)

On-demand TPU v5e (1x1 topology, single chip). Uses the `west1` context.

```bash
kubectl --context=west1 run POD_NAME --image=IMAGE_HERE --namespace=development --restart=Never --overrides='{"spec": {"nodeSelector": {"cloud.google.com/gke-tpu-accelerator": "tpu-v5-lite-podslice", "cloud.google.com/gke-tpu-topology": "1x1"}, "affinity": {"nodeAffinity": {"requiredDuringSchedulingIgnoredDuringExecution": {"nodeSelectorTerms": [{"matchExpressions": [{"key": "flyte.org/node-role", "operator": "In", "values": ["worker"]}, {"key": "cloud.google.com/gke-spot", "operator": "DoesNotExist"}]}]}}}, "tolerations": [{"key": "flyte.org/node-role", "operator": "Equal", "value": "worker", "effect": "NoSchedule"}, {"key": "google.com/tpu", "operator": "Equal", "value": "present", "effect": "NoSchedule"}], "containers": [{"name": "tpu-container", "image": "IMAGE_HERE", "command": ["sleep", "SECONDS"], "resources": {"requests": {"cpu": "23", "memory": "43200M", "ephemeral-storage": "10G", "google.com/tpu": "1"}, "limits": {"cpu": "24", "memory": "43200M", "ephemeral-storage": "10G", "google.com/tpu": "1"}}, "securityContext": {"privileged": true}}]}}' -- sleep SECONDS
```

**TPU v5e Pod Notes:**
- **Cluster**: west1
- **TPU**: v5e (tpu-v5-lite-podslice), 1x1 topology, single chip
- **Node**: ct5lp-highcpu-1t (24 vCPU, 48GB RAM)
- **Node selector labels**:
  - `cloud.google.com/gke-tpu-accelerator: tpu-v5-lite-podslice`
  - `cloud.google.com/gke-tpu-topology: 1x1`
- **Tolerations**:
  - `flyte.org/node-role=worker:NoSchedule` (standard worker toleration)
  - `google.com/tpu=present:NoSchedule` (TPU-specific)
- **Resource request**: `google.com/tpu: 1`
- **Affinity**: On-demand only (`cloud.google.com/gke-spot: DoesNotExist`)
- **Use case**: JAX/TPU ML training and inference

### 3. Create CPU Pod (--type cpu)

For large CPU workloads (Rosetta, simulations, etc.):

**west1 cluster nodepool names:**
| Instance Type | Nodepool Name (on-demand) | Nodepool Name (spot) |
|--------------|---------------------------|----------------------|
| t2d-standard-60 | `t2dstandard6020250716193052161300000004` | `t2dstandard602025071619320337750000000a` |
| t2d-standard-48 | `t2dstandard4820250716193158564400000009` | `t2dstandard4820250716193121370400000006` |

**Taints on t2d nodepools:**
- `flyte.org/node-role=worker:NoSchedule`
- `union.ai/multithreading=disabled:NoSchedule`

```bash
kubectl --context=west1 run POD_NAME --image=IMAGE_HERE --namespace=development --restart=Never --overrides='{"spec": {"nodeSelector": {"cloud.google.com/gke-nodepool": "t2dstandard6020250716193052161300000004"}, "tolerations": [{"key": "flyte.org/node-role", "operator": "Equal", "value": "worker", "effect": "NoSchedule"}, {"key": "union.ai/multithreading", "operator": "Equal", "value": "disabled", "effect": "NoSchedule"}], "containers": [{"name": "cpu-container", "image": "IMAGE_HERE", "command": ["sleep", "SECONDS"], "env": [{"name": "PATH", "value": "/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin"}], "resources": {"requests": {"cpu": "58", "memory": "230G", "ephemeral-storage": "20G"}, "limits": {"cpu": "60", "memory": "240G", "ephemeral-storage": "30G"}}, "securityContext": {"privileged": true}}]}}' -- sleep SECONDS
```

### 4. Replace Values

- `POD_NAME` with the generated or specified pod name
- `IMAGE_HERE` with the specified image or the default image
- `SECONDS` with the number of seconds (hours * 3600):
  - 1 hour = 3600 seconds
  - 2 hours = 7200 seconds
  - 4 hours = 14400 seconds

### 5. Report Pod Info

After creating the pod, report:
1. **Pod name** - So user can reference it later
2. **Pod type** - GPU or CPU
3. **Check status command**: `kubectl --context=west1 get pod POD_NAME -n development -o wide`
4. **Exec command**: `kubectl --context=west1 exec -it POD_NAME -n development -- bash`
5. **Delete command**: `kubectl --context=west1 delete pod POD_NAME -n development`

If pending, note it may take 2-5 minutes to provision a new node.

**IMPORTANT**: Do NOT delete the pod when done. Let the sleep timer expire naturally so other workloads can reuse it. Only delete if explicitly requested by the user.

## Rsync Files to a Pod

To sync local files to a pod using rsync over kubectl:

### 1. Install rsync on the pod
```bash
kubectl --context=west1 exec POD_NAME -n development -- apt-get update -qq
kubectl --context=west1 exec POD_NAME -n development -- apt-get install -y -qq rsync
```

### 2. Create a wrapper script
Write this to `/tmp/claude/tmp/kubectl-rsync.sh`:
```bash
#!/bin/bash
# Wrapper for rsync -e to tunnel through kubectl exec
# rsync passes: $0=this_script, $1=host (ignored), rest=rsync server command
shift  # drop the "host" argument rsync passes
kubectl --context=west1 exec -i POD_NAME -n development -- "$@"
```
Replace `POD_NAME` with the actual pod name, then `chmod +x` it.

### 3. Rsync files
```bash
# Create target directory on pod
kubectl --context=west1 exec POD_NAME -n development -- mkdir -p /workspace/myproject

# Sync (the "pod:" is a dummy host label, the wrapper handles routing)
rsync -avz --blocking-io \
  --exclude='.venv' --exclude='.git' --exclude='__pycache__' \
  --exclude='.pytest_cache' --exclude='*.pyc' \
  -e /tmp/claude/tmp/kubectl-rsync.sh \
  ./ pod:/workspace/myproject/
```

## List Running Pods

```bash
# All compute pods
kubectl --context=west1 get pods -n development -l '!job-name' | grep -E '(gpu-|tpu-|cpu-)'

# GPU pods only
kubectl --context=west1 get pods -n development | grep gpu-

# TPU pods only
kubectl --context=west1 get pods -n development | grep tpu-

# CPU pods only
kubectl --context=west1 get pods -n development | grep cpu-
```

## Pod Specifications

### GPU Pod (g2-standard-16)
- **Node**: g2-standard-16
- **CPU**: 15 cores (request), 16 cores (limit)
- **Memory**: 57.6 GB
- **Ephemeral Storage**: 30 GB
- **GPU**: 1x NVIDIA L4 (24GB VRAM)
- **Cluster**: west1
- **Best for**: Deep learning, GPU-accelerated docking, structure prediction

### RTX 6000 Pod (g4-standard-48)
- **Node**: g4-standard-48
- **CPU**: 46 cores (request), 48 cores (limit)
- **Memory**: 150 GB (request), 160 GB (limit)
- **Ephemeral Storage**: 35 GB (request), 40 GB (limit)
- **GPU**: 1x NVIDIA RTX Pro 6000 Blackwell Server Edition (96GB VRAM)
- **Cluster**: central1 (NOT west1)
- **Compute Capability**: sm_120 (Blackwell)
- **CUDA**: Requires 12.8+ for Blackwell support
- **Best for**: Large model inference, Blackwell testing, gnina builds with sm_120

### TPU v5e Pod (ct5lp-highcpu-1t)
- **Node**: ct5lp-highcpu-1t
- **CPU**: 23 cores (request), 24 cores (limit)
- **Memory**: 43.2 GB
- **Ephemeral Storage**: 10 GB
- **TPU**: 1x TPU v5e chip (1x1 topology)
- **Cluster**: west1
- **Best for**: JAX-based ML training, TPU inference workloads

### CPU Pod (t2d-standard-60)
- **Node**: t2d-standard-60
- **CPU**: 58 cores (request), 60 cores (limit)
- **Memory**: 230 GB (request), 240 GB (limit)
- **Ephemeral Storage**: 50 GB
- **Best for**: Rosetta/PyRosetta, MD simulations, parallel CPU workloads

## Namespace
- **Namespace**: development
- **Cluster**: west1 (GKE)
- **Kubecontext**: `west1`

## IMPORTANT: Always Use West1 Context

**ALL kubectl commands in this skill MUST include `--context=west1`** to ensure they target the correct cluster, regardless of the user's current kubectl context.

Example:
```bash
kubectl --context=west1 get pods -n development
kubectl --context=west1 run POD_NAME ...
kubectl --context=west1 exec POD_NAME ...
```
