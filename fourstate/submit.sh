#!/bin/bash
# Submit every unfinished run under a folder.
#   bash submit.sh RUNS [--exe vasp_ncl|vasp_std] [--job templates/job.sh] [--dry]
# Skips runs with STATUS=DONE and runs already queued (job name = folder tag).
# Set SUBMIT="sg mygroup -c" etc. if your cluster needs a wrapper around sbatch.
set -e
TOP=$(realpath "$1"); shift
EXE=vasp_ncl; JOB=$(dirname "$(realpath "$0")")/../templates/job.sh; DRY=
while [ $# -gt 0 ]; do case $1 in
  --exe) EXE=$2; shift 2;; --job) JOB=$(realpath "$2"); shift 2;; --dry) DRY=1; shift;;
  *) echo "unknown option $1"; exit 1;; esac; done
queued=$(squeue -u "$USER" -h -o '%j' 2>/dev/null || true)
n=0
for d in $(find "$TOP" -name state.json -printf '%h
' | sort); do
  [ -f "$d/STATUS" ] && grep -q DONE "$d/STATUS" && continue
  tag=$(basename "$(dirname "$(dirname "$d")")")_$(basename "$(dirname "$d")" | cut -d_ -f1)_$(basename "$d")
  echo "$queued" | grep -qx "$tag" && continue
  sed "s/^#SBATCH --job-name=.*/#SBATCH --job-name=$tag/" "$JOB" > "$d/job.sh"
  if [ -z "$DRY" ]; then
    if [ -n "$SUBMIT" ]; then (cd "$d" && $SUBMIT "sbatch --export=ALL,VASP_EXE=$EXE job.sh")
    else (cd "$d" && sbatch --export=ALL,VASP_EXE=$EXE job.sh); fi
  fi
  n=$((n+1))
done
echo "$n runs ${DRY:+(dry run, nothing submitted)}"
