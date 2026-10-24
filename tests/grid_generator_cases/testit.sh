#!/bin/bash
#
# End-to-end test of grid-generator output with an FFT-accelerated SNaC build:
# valid generated cases must run with the generated grids, and invalid grids
# must stop the solver with a clear message.
#
set -e

TESTDIR=$(pwd)
SNACDIR=$(pwd)/../..
RUNDIR=$SNACDIR/run
EXEC=snac
MPIRUN_OPTIONS='--oversubscribe'
if mpirun --version 2>&1 | grep -qi intel; then MPIRUN_OPTIONS=""; fi
MPIRUN="mpirun -n 2 $MPIRUN_OPTIONS"
CASES=$(mktemp -d)

rm -rf $RUNDIR
echo "Compiling ..."
cd $SNACDIR && make clean && make FFT_AXIS=3 -j run
python "$TESTDIR/make_cases.py" "$CASES"

# run_case <case> <expected exit: ok|error> [expected stderr message]
run_case() {
  local case=$CASES/$1
  cp "$RUNDIR/$EXEC" "$TESTDIR/dns.nml" "$case/"
  mkdir -p "$case/data"
  cd "$case"
  echo "Running SNaC on $1 ..."
  set +e
  ${MPIRUN} ./${EXEC} 1> log_file.log 2> err_log.log
  local status=$?
  set -e
  if [ "$2" = "ok" ] && [ $status -ne 0 ]; then
    cat err_log.log
    echo "FAILED: $1 should run"
    exit 1
  fi
  if [ "$2" = "error" ]; then
    if [ $status -eq 0 ]; then
      echo "FAILED: $1 should stop with an error"
      exit 1
    fi
    if ! grep -q "$3" err_log.log; then
      cat err_log.log
      echo "FAILED: $1 did not report: $3"
      exit 1
    fi
  fi
  echo "$1: as expected"
  cd "$TESTDIR"
}

run_case valid ok
python "$TESTDIR/check_grid.py" "$CASES/valid"
run_case stretched_fft error "grid spacing in the FFT direction must be uniform"
run_case legacy_geometric ok
run_case geometric_without_file error "geometric grid mappings (gt = 4 to 7) require gr > 0"

rm -rf "$CASES" $RUNDIR
echo "Passed!"
