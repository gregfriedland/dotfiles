# GCP Management

Common GCP management commands and utilities for the Rezo organization.

## Usage
```
/gcp [subcommand] [options]
```

## Subcommands

| Subcommand | Description |
|------------|-------------|
| `groups` | List all Google Workspace groups |
| `members GROUP` | List members of a specific group |
| `projects` | List GCP projects |
| `billing` | Show billing summary |
| `nat` | Check NAT traffic and costs |

## Instructions

### /gcp groups

List all Google Workspace groups in the organization:

```bash
gcloud identity groups search --customer=C01vjieca --labels='cloudidentity.googleapis.com/groups.discussion_forum' --view=FULL --page-size=100 --format='table(displayName,groupKey.id,description)'
```

### /gcp members GROUP

List members of a specific group. GROUP can be just the name (e.g., `agentspace_admins`) or full email.

```bash
gcloud identity groups memberships list --group-email=GROUP@rezotx.com --format='table(preferredMemberKey.id,roles.name)'
```

Example:
```bash
# List AgentSpace Admins members
gcloud identity groups memberships list --group-email=agentspace_admins@rezotx.com --format='table(preferredMemberKey.id,roles.name)'
```

### /gcp projects

List all GCP projects:
```bash
gcloud projects list --format='table(projectId,name,projectNumber)'
```

### /gcp billing

Query billing data from BigQuery (last 30 days by service):
```bash
bq query --use_legacy_sql=false --format=prettyjson '
SELECT
  service.description as service,
  ROUND(SUM(cost), 2) as total_cost,
  ROUND(SUM(usage.amount_in_pricing_units), 2) as usage_amount,
  usage.pricing_unit
FROM `gke-test-421317.all_billing_data_detailed.gcp_billing_export_resource_v1_0164FD_86FF59_10FC73`
WHERE DATE(usage_start_time) >= DATE_SUB(CURRENT_DATE(), INTERVAL 30 DAY)
GROUP BY service.description, usage.pricing_unit
ORDER BY total_cost DESC
LIMIT 20
'
```

### /gcp nat

Check Cloud NAT traffic for the CI cluster:

1. **Check NAT router status**:
```bash
gcloud compute routers describe ci-cluster-router --region=us-central1 --project=ci-tests-1234 --format='yaml(nats)'
```

2. **Check VPC flow logs** (last 10 minutes):
```bash
gcloud logging read 'resource.type="gce_subnetwork" AND logName:"vpc_flows"' --project=ci-tests-1234 --limit=100 --freshness=10m --format='csv(jsonPayload.connection.dest_ip,jsonPayload.bytes_sent)'
```

3. **Check Private Google Access DNS zones**:
```bash
gcloud dns managed-zones list --project=ci-tests-1234 --format='table(name,dnsName,visibility)'
```

## Common Tasks

### Add user to a group
```bash
gcloud identity groups memberships add --group-email=GROUP@rezotx.com --member-email=USER@rezotx.com
```

### Create a new group
```bash
gcloud identity groups create GROUP@rezotx.com --organization=rezotx.com --display-name="Display Name" --description="Description"
```

### Check IAM policy for a project
```bash
gcloud projects get-iam-policy PROJECT_ID --format='table(bindings.role,bindings.members)'
```

## Rsync to GCE Instances via IAP

To use rsync with GCE instances that require IAP tunneling:

### 1. Configure SSH for GCE

```bash
gcloud compute config-ssh --project=gke-test-421317
```

### 2. Add IAP ProxyCommand to SSH config

Edit `~/.ssh/config` and add the `ProxyCommand` and `User` lines to the instance entry:

```
Host INSTANCE_NAME.ZONE.PROJECT_ID
    HostName EXTERNAL_IP
    User greg_rezotx_com  # OS Login username
    IdentityFile ~/.ssh/google_compute_engine
    UserKnownHostsFile=~/.ssh/google_compute_known_hosts
    HostKeyAlias=compute.INSTANCE_ID
    IdentitiesOnly=yes
    CheckHostIP=no
    ProxyCommand gcloud compute start-iap-tunnel INSTANCE_NAME 22 --listen-on-stdin --zone=ZONE --project=PROJECT_ID
```

### 3. Use rsync normally

```bash
rsync -avz --progress /local/path/ INSTANCE_NAME.ZONE.PROJECT_ID:~/remote/path/
```

### Example

For instance `gnina-profile-vm` in `us-west1-a`:

```bash
# SSH config entry
Host gnina-profile-vm.us-west1-a.gke-test-421317
    User greg_rezotx_com
    ProxyCommand gcloud compute start-iap-tunnel gnina-profile-vm 22 --listen-on-stdin --zone=us-west1-a --project=gke-test-421317
    ...

# Rsync command
rsync -avz --progress --exclude='.git' /local/src/ gnina-profile-vm.us-west1-a.gke-test-421317:~/src/
```

## Organization Info

- **Customer ID**: C01vjieca
- **Organization**: rezotx.com
- **Billing Dataset**: `gke-test-421317.all_billing_data_detailed`
