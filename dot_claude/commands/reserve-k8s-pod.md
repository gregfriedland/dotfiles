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
| `cpu` | t2d-standard-60 | 60 vCPU, 240GB RAM | CPU-intensive work (Rosetta, simulations) |

**Default image:** `us-central1-docker.pkg.dev/gke-test-421317/flyte/internal-base:linux-cuda-bd58cd90987a6a20913194af9a19a455551a44a5`

## Common Images

| Image | Description |
|-------|-------------|
| `boltz/boltz` | Boltz-2 protein structure prediction |
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
```

## Instructions

When this command is invoked:

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
| `us-central1-docker.pkg.dev/gke-test-421317/flyte/gnina-build-base:*` | `gnina-build` |
| `nvcr.io/nvidia/pytorch:*` | `pytorch` |
| (default) | `base` |

Example generated names:
- `gpu-boltz-a3f8k2`
- `cpu-base-x9m4p7`

If `--name` is provided, use that name directly instead.

### 2. Create GPU Pod (--type gpu)

```bash
kubectl run POD_NAME --image=IMAGE_HERE --namespace=development --restart=Never --overrides='{"spec": {"nodeSelector": {"cloud.google.com/gke-accelerator": "nvidia-l4"}, "affinity": {"nodeAffinity": {"requiredDuringSchedulingIgnoredDuringExecution": {"nodeSelectorTerms": [{"matchExpressions": [{"key": "flyte.org/node-role", "operator": "In", "values": ["worker"]}, {"key": "cloud.google.com/gke-spot", "operator": "DoesNotExist"}]}]}}}, "tolerations": [{"key": "flyte.org/node-role", "operator": "Equal", "value": "worker", "effect": "NoSchedule"}, {"key": "nvidia.com/gpu", "operator": "Equal", "value": "present", "effect": "NoSchedule"}, {"key": "nvidia.com/gpu", "operator": "Exists", "effect": "NoSchedule"}], "containers": [{"name": "gpu-container", "image": "IMAGE_HERE", "command": ["sleep", "SECONDS"], "env": [{"name": "LD_LIBRARY_PATH", "value": "/usr/local/nvidia/lib64"}, {"name": "PATH", "value": "/usr/local/nvidia/bin:/usr/local/cuda/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin"}], "resources": {"requests": {"cpu": "15", "memory": "57600M", "ephemeral-storage": "30G", "nvidia.com/gpu": "1"}, "limits": {"cpu": "16", "memory": "57600M", "ephemeral-storage": "30G", "nvidia.com/gpu": "1"}}, "securityContext": {"privileged": true}}]}}' -- sleep SECONDS
```

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
kubectl run POD_NAME --image=IMAGE_HERE --namespace=development --restart=Never --overrides='{"spec": {"nodeSelector": {"cloud.google.com/gke-nodepool": "t2dstandard6020250716193052161300000004"}, "tolerations": [{"key": "flyte.org/node-role", "operator": "Equal", "value": "worker", "effect": "NoSchedule"}, {"key": "union.ai/multithreading", "operator": "Equal", "value": "disabled", "effect": "NoSchedule"}], "containers": [{"name": "cpu-container", "image": "IMAGE_HERE", "command": ["sleep", "SECONDS"], "env": [{"name": "PATH", "value": "/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin"}], "resources": {"requests": {"cpu": "58", "memory": "230G", "ephemeral-storage": "20G"}, "limits": {"cpu": "60", "memory": "240G", "ephemeral-storage": "30G"}}, "securityContext": {"privileged": true}}]}}' -- sleep SECONDS
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
3. **Check status command**: `kubectl get pod POD_NAME -n development -o wide`
4. **Exec command**: `kubectl exec -it POD_NAME -n development -- bash`
5. **Delete command**: `kubectl delete pod POD_NAME -n development`

If pending, note it may take 2-5 minutes to provision a new node.

**IMPORTANT**: Do NOT delete the pod when done. Let the sleep timer expire naturally so other workloads can reuse it. Only delete if explicitly requested by the user.

## List Running Pods

```bash
# All compute pods
kubectl get pods -n development -l '!job-name' | grep -E '(gpu-|cpu-)'

# GPU pods only
kubectl get pods -n development | grep gpu-

# CPU pods only
kubectl get pods -n development | grep cpu-
```

## Pod Specifications

### GPU Pod (g2-standard-16)
- **Node**: g2-standard-16
- **CPU**: 15 cores (request), 16 cores (limit)
- **Memory**: 57.6 GB
- **Ephemeral Storage**: 30 GB
- **GPU**: 1x NVIDIA L4 (24GB VRAM)
- **Best for**: Deep learning, GPU-accelerated docking, structure prediction

### CPU Pod (t2d-standard-60)
- **Node**: t2d-standard-60
- **CPU**: 58 cores (request), 60 cores (limit)
- **Memory**: 230 GB (request), 240 GB (limit)
- **Ephemeral Storage**: 50 GB
- **Best for**: Rosetta/PyRosetta, MD simulations, parallel CPU workloads

## Namespace
- **Namespace**: development
- **Cluster**: west1 (GKE)
