# K8s Pod Operations

Utilities for working with Kubernetes pods, including file sync via rsync.

## Rsync to K8s Pod

Sync local files to a running K8s pod using rsync over kubectl exec.

### Usage
```
/k8s rsync <local_path> <pod_name> <remote_path> [--namespace NS] [--exclude PATTERN]
```

### How It Works

rsync doesn't natively support kubectl, so we create a wrapper script that uses `kubectl exec -i` as the remote shell transport:

```bash
# 1. Ensure rsync is installed on pod
kubectl exec "$POD_NAME" -n "$NAMESPACE" -- which rsync >/dev/null 2>&1 || {
    echo "Installing rsync on pod..."
    kubectl exec "$POD_NAME" -n "$NAMESPACE" -- apt-get update -qq
    kubectl exec "$POD_NAME" -n "$NAMESPACE" -- apt-get install -y -qq rsync
}

# 2. Create wrapper script for rsync to use kubectl exec as transport
# rsync calls: <rsh> <host> <command...>
# We ignore <host> and pass <command...> to kubectl exec
KRSYNC="/tmp/krsync-$$.sh"
cat > "$KRSYNC" << 'EOF'
#!/bin/bash
shift  # skip the dummy host argument
kubectl exec -i POD_NAME -n NAMESPACE -- "$@"
EOF
# Replace POD_NAME and NAMESPACE in the script
sed -i '' "s/POD_NAME/$POD_NAME/g; s/NAMESPACE/$NAMESPACE/g" "$KRSYNC"
chmod +x "$KRSYNC"

# 3. Make destination writable (if needed)
kubectl exec "$POD_NAME" -n "$NAMESPACE" -- chmod -R 777 "$REMOTE_PATH" 2>/dev/null || true

# 4. Run rsync with kubectl wrapper as remote shell
# "pod:" is a dummy host required by rsync syntax (ignored by wrapper)
rsync -avR --progress --checksum \
    --exclude='*.o' \
    --exclude='.git' \
    --exclude='build' \
    --exclude='__pycache__' \
    -e "$KRSYNC" \
    "$LOCAL_PATH" \
    pod:"$REMOTE_PATH"

# 5. Cleanup
rm -f "$KRSYNC"
```

### Key Flags

| Flag | Purpose |
|------|---------|
| `-a` | Archive mode (preserves permissions, timestamps, etc.) |
| `-v` | Verbose output |
| `-R` | Use relative paths (preserves directory structure) |
| `--progress` | Show transfer progress |
| `--checksum` | Compare by checksum not timestamp (more reliable over kubectl) |
| `--delete` | Remove files on remote that don't exist locally |
| `-e "$KRSYNC"` | Use our kubectl wrapper as remote shell |

### Common Excludes

```bash
--exclude='*.o'           # Object files
--exclude='*.a'           # Archives
--exclude='.git'          # Git directory
--exclude='build'         # Build directories
--exclude='build-docker'  # Docker build dirs
--exclude='__pycache__'   # Python cache
--exclude='node_modules'  # Node.js deps
--exclude='.venv'         # Python venvs
```

### Examples

```bash
# Sync current directory to pod
/k8s rsync . gpu-pod-abc123 /workspace --namespace development

# Sync specific files with exclusions
/k8s rsync ./src gpu-pod-abc123 /app/src --exclude='*.pyc' --exclude='.git'

# Sync with delete (mirror local to remote)
/k8s rsync ./code gpu-pod-abc123 /code --delete
```

### Troubleshooting

1. **"rsync: command not found"** - Install rsync on pod first:
   ```bash
   kubectl exec POD -n NS -- apt-get update && apt-get install -y rsync
   ```

2. **Permission denied** - Make destination writable:
   ```bash
   kubectl exec POD -n NS -- chmod -R 777 /destination/path
   ```

3. **Slow transfers** - Use `--checksum` instead of timestamp comparison (timestamps unreliable over kubectl)

## Other K8s Operations

### Copy single file (without rsync)
```bash
kubectl cp local/file.txt POD_NAME:/remote/path -n NAMESPACE
kubectl cp POD_NAME:/remote/file.txt local/path -n NAMESPACE
```

### Execute command on pod
```bash
kubectl exec -it POD_NAME -n NAMESPACE -- bash
kubectl exec POD_NAME -n NAMESPACE -- command arg1 arg2
```

### Get pod status
```bash
kubectl get pod POD_NAME -n NAMESPACE -o wide
kubectl describe pod POD_NAME -n NAMESPACE
```

### Stream logs
```bash
kubectl logs -f POD_NAME -n NAMESPACE
```
