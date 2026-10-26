#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'EOF'
Usage:
  neb-recal.sh JOB_DIR
      Auto-detect numeric image directories in JOB_DIR and skip endpoints.

  neb-recal.sh JOB_DIR START END
      Process the inclusive range START..END in JOB_DIR (e.g. 01 08).

For every selected image, this script:
  1. verifies that CONTCAR and POSCAR are nonempty;
  2. backs up POSCAR as POSCAR.<timestamp>.bak;
  3. copies CONTCAR to POSCAR.

Run it from the work directory containing JOB_DIR. JOB_DIR may also be an
absolute or relative path.
It does not delete VASP output files and does not modify CONTCAR.
EOF
}

die() {
    printf 'ERROR: %s\n' "$*" >&2
    exit 1
}

[[ $# -eq 1 || $# -eq 3 ]] || {
    usage >&2
    exit 2
}

job_dir=${1%/}
[[ -n "$job_dir" ]] || die 'JOB_DIR must not be empty'
[[ -d "$job_dir" ]] || die "job directory not found: $job_dir"

declare -a images=()

if [[ $# -eq 1 ]]; then
    while IFS= read -r dir; do
        images+=("$dir")
    done < <(
        for path in "$job_dir"/[0-9][0-9]; do
            [[ -d "$path" ]] && printf '%s\n' "${path##*/}"
        done | sort
    )

    ((${#images[@]} >= 3)) || \
        die 'at least three numeric image directories are required for auto-detection'

    # The first and last directories are fixed NEB endpoints.
    images=("${images[@]:1:${#images[@]}-2}")
else
    start=$2
    end=$3
    [[ $start =~ ^[0-9]+$ && $end =~ ^[0-9]+$ ]] || \
        die 'START and END must be nonnegative integers, e.g. 01 08'

    start_num=$((10#$start))
    end_num=$((10#$end))
    ((start_num <= end_num)) || die "START ($start) is greater than END ($end)"

    width=${#start}
    ((${#end} > width)) && width=${#end}

    for ((i = start_num; i <= end_num; i++)); do
        printf -v dir "%0${width}d" "$i"
        images+=("$dir")
    done
fi

printf 'Selected images:'
printf ' %s' "${images[@]}"
printf '\n'
printf 'Job directory: %s\n' "$job_dir"

# Validate the full selection before changing any file.
for dir in "${images[@]}"; do
    image_dir="$job_dir/$dir"
    [[ -d "$image_dir" ]] || die "directory not found: $image_dir"
    [[ -s "$image_dir/CONTCAR" ]] || die "missing or empty file: $image_dir/CONTCAR"
    [[ -s "$image_dir/POSCAR" ]] || die "missing or empty file: $image_dir/POSCAR"
done

timestamp=$(date +%Y%m%d_%H%M%S)

for dir in "${images[@]}"; do
    image_dir="$job_dir/$dir"
    backup="$image_dir/POSCAR.$timestamp.bak"
    cp -p -- "$image_dir/POSCAR" "$backup"
    cp -p -- "$image_dir/CONTCAR" "$image_dir/POSCAR"

    if cmp -s -- "$image_dir/CONTCAR" "$image_dir/POSCAR"; then
        printf 'Updated %-4s (backup: %s)\n' "$dir" "$backup"
    else
        die "verification failed after updating $image_dir/POSCAR"
    fi
done

printf 'NEB restart structures prepared successfully.\n'
