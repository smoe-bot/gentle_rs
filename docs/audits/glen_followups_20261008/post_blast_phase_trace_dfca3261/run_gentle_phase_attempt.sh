#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 4 ]]; then
  echo "usage: $0 CASE MODE HANDOFF MANIFEST" >&2
  exit 64
fi

case_name=$1
mode=$2
handoff=$3
manifest=$4
evidence_root=/tmp/gentle-phase-trace-dfca3261-evidence
cli=/tmp/gentle-phase-trace-dfca3261-target/debug/gentle_cli

case "$mode" in
  baseline|trace) ;;
  *) echo "unsupported mode: $mode" >&2; exit 64 ;;
esac

run_dir=$(mktemp -d "$evidence_root/$case_name-$mode.XXXXXX")
state="$run_dir/state.gentle.json"
[[ ! -e "$state" ]]
cmd=("$cli" --state "$state" primers specificity-multi-import "$handoff" "$manifest")

if [[ "$mode" == trace ]]; then
  env_args=(env GENTLE_PRIMER_IMPORT_TRACE=1)
  printf 'GENTLE_PRIMER_IMPORT_TRACE=1\n' > "$run_dir/environment.txt"
else
  env_args=(env -u GENTLE_PRIMER_IMPORT_TRACE)
  printf 'GENTLE_PRIMER_IMPORT_TRACE=<unset>\n' > "$run_dir/environment.txt"
fi

printf '%q ' /usr/bin/time -v -o "$run_dir/time.txt" timeout --signal=TERM --kill-after=5s 60s "${env_args[@]}" "${cmd[@]}" > "$run_dir/command.txt"
printf '\n' >> "$run_dir/command.txt"
printf 'absent\n' > "$run_dir/state-before.txt"

set +e
/usr/bin/time -v -o "$run_dir/time.txt" \
  timeout --signal=TERM --kill-after=5s 60s \
  "${env_args[@]}" "${cmd[@]}" \
  >"$run_dir/stdout.json" 2>"$run_dir/stderr.txt"
rc=$?
set -e
printf '%s\n' "$rc" > "$run_dir/exit-code.txt"

if [[ -e "$state" ]]; then
  stat -c 'present %s' "$state" > "$run_dir/state-after.txt"
  sha256sum "$state" >> "$run_dir/state-after.txt"
else
  printf 'absent\n' > "$run_dir/state-after.txt"
fi

pgrep -af "$state" > "$run_dir/processes-after.txt" || true
if [[ -s "$run_dir/processes-after.txt" ]]; then
  echo "run-owned process survived: $run_dir" >&2
  exit 90
fi

printf '%s\n' "$run_dir" > "$evidence_root/$case_name-$mode.latest"
printf 'run=%s rc=%s state=%s stderr_lines=%s stdout_bytes=%s\n' \
  "$run_dir" "$rc" "$(head -1 "$run_dir/state-after.txt")" \
  "$(wc -l < "$run_dir/stderr.txt")" "$(stat -c %s "$run_dir/stdout.json")"
