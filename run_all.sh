#!/usr/bin/env bash
# ==============================================================================
# run_all.sh
# Full rerun in dependency order (the project's pipeline driver; `targets` is not used).
# Interruptible and resumable: every P8 replicate and every P9
# fit is saved as its own file the moment it finishes, and a rerun skips finished work.
#
# Usage: ./run_all.sh                  all steps
#        ./run_all.sh 7                steps 7 to 10
#        ./run_all.sh 8 3 4 9 10       exactly these steps, in this order
#        ./run_all.sh status           what is saved and what is left
#        ./run_all.sh stop             ask a running job to finish its current tasks and stop
#
# Ending the day:   ./run_all.sh stop         (or start the run with CPD_STOP_AT=18:30)
#   Workers finish the tasks already in progress (a few minutes), save them and exit; nothing
#   is lost. Next day: run the same command again.
# A hard kill (Ctrl-C, pkill) is also safe but loses the tasks in progress (minutes of work).
#
# Environment: CPD_MAX_CORES (hard cap, default 24), N_WORKERS (default 16),
#              CPD_STOP_AT ("HH:MM" or "YYYY-MM-DD HH:MM:SS"), REDO_REAL=1 to rerun step 7.
#
#   1  P3 gate (Stan reference, SBC)                       ~40 min
#   2  Geweke Getting-it-Right test                        ~6 min
#   3  Long well-log run (P4 mixing and speed)             ~10 min on an idle machine
#   4  P4 gate                                             ~30 min
#   5  P5, P6, P7 gates                                    ~20 min
#   6  Full simulation study (P8)                          ~9-11 h
#   7  Real data at the converged length (+ BTC extension) skipped if output/real_data_long.rds exists
#   8  Replicated ablation and sensitivity (P9)            ~3.5 h
#   9  P8 and P9 gates
#  10  Tables and figures, P10 gate, unit tests
# ==============================================================================
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p output

p8_total=3400; p9_total=1120
count() { ls "$1"/*.rds 2>/dev/null | wc -l; }

case "${1:-}" in
  stop)
    touch STOP
    echo "STOP requested. Running jobs finish their current tasks, save them and exit (a few minutes)."
    echo "Resume later with the same command; run_all.sh clears the STOP file when it starts."
    exit 0 ;;
  status)
    echo "P8 full study        : $(count output/p8_full) / $p8_total replicate files saved"
    echo "P9 replicated        : $(count output/p9_replicated) / $p9_total fit files saved"
    for f in geweke_test_results welllog_long_run real_data_long gate_p4_results gate_p8_results gate_p9_results; do
      if [ -f "output/$f.rds" ]; then echo "  $f.rds   $(date -r output/$f.rds '+%F %H:%M')"; else echo "  $f.rds   (missing)"; fi
    done
    pgrep -f 'run_p8_full|run_p9_replicated|run_real_data|run_welllog|gate_p[0-9]|test_geweke' >/dev/null && echo "A run is active." || echo "No run is active."
    [ -f STOP ] && echo "STOP file present." || true
    exit 0 ;;
esac

# Never more than CPD_MAX_CORES of the 32 logical cores (thermal limit); BLAS/OpenMP single-threaded
# so that worker count equals core use. Measured throughput on this machine saturates at about
# 16-20 workers (24 workers are no faster), so 16 is the default.
export CPD_MAX_CORES=${CPD_MAX_CORES:-24}
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
W=${N_WORKERS:-16}
[ "$W" -gt "$CPD_MAX_CORES" ] && W=$CPD_MAX_CORES
export CPD_STOP_AT=${CPD_STOP_AT:-}

exec > >(tee -a output/run_all.log) 2>&1
log() { echo "[$(date '+%F %T')] $*"; }
[ -f STOP ] && { rm -f STOP; log "cleared an old STOP file"; }
log "run_all.sh $* | workers $W | core cap $CPD_MAX_CORES | stop at '${CPD_STOP_AT:-never}'"
if [ "$#" -le 1 ]; then ORDER=$(seq "${1:-1}" 10); else ORDER="$*"; fi

# Rscript wrapper: exit status 75 means "stopped on request" (everything saved)
rs() {
  local rc=0
  Rscript "$@" || rc=$?
  if [ "$rc" -eq 75 ]; then
    log "Stopped on request. All finished work is saved. Rerun './run_all.sh $*' to resume." >&2
    log "(resume the same command you started with)"
    exit 0
  elif [ "$rc" -ne 0 ]; then
    log "FAILED: Rscript $* (exit $rc)"; exit "$rc"
  fi
}

run_step() {
  case "$1" in
    1) log "1 P3 gate";           rs tests/gate_p3_verification.R ;;
    2) log "2 Geweke test";       GEWEKE_WORKERS=$W rs tests/test_geweke_getting_it_right.R ;;
    3) log "3 long well-log run"; rs tests/run_welllog_long.R ;;
    4) log "4 P4 gate";           rs tests/gate_p4_verification.R ;;
    5) log "5 P5-P7 gates";       rs tests/gate_p5_verification.R; rs tests/gate_p6_verification.R
                                  rs tests/gate_p7_verification.R ;;
    6) log "6 P8 full study";     rs tests/run_p8_full.R "$W" ;;
    7) if [ -f output/real_data_long.rds ] && [ -z "${REDO_REAL:-}" ]; then
         log "7 real data: output/real_data_long.rds exists, skipped (REDO_REAL=1 to rerun)"
       else
         log "7 real data"; rs tests/run_real_data_long.R; rs tests/run_real_data_btc_extend.R
       fi ;;
    8) log "8 P9 replicated";     rs tests/run_p9_replicated.R "$W" ;;
    9) log "9 P8, P9 gates";      rs tests/gate_p8_verification.R; rs tests/gate_p9_verification.R ;;
    10) log "10 outputs";         rs run_pipeline.R; rs tests/gate_p10_verification.R
        rs -e 'r <- as.data.frame(testthat::test_dir("tests/testthat", reporter = "silent", stop_on_failure = FALSE)); cat("unit tests:", sum(r$nb), "checks,", sum(r$failed), "failed,", sum(r$error), "errors\n"); if (sum(r$failed) + sum(r$error) > 0) quit(status = 1)' ;;
    *) echo "unknown step $1" >&2; exit 2 ;;
  esac
  # a STOP file or deadline also ends the sequence between steps
  if [ -f STOP ]; then log "STOP file found after step $1; resume by rerunning the same command."; exit 0; fi
}

for k in $ORDER; do run_step "$k"; done
log "done"
