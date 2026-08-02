#!/bin/bash
set -e

TESTDIR=$(pwd)
SNACDIR=$(pwd)/../..
RUNDIR=$SNACDIR/run
MPIRUN_OPTIONS='--oversubscribe'
if mpirun --version 2>&1 | grep -qi intel; then MPIRUN_OPTIONS=""; fi

rm -rf "$RUNDIR"
cd "$SNACDIR"
make clean
make FFT_AXIS=0 -j run
cp "$TESTDIR/dns.nml" "$RUNDIR/dns.nml"
cp "$TESTDIR/blocks.nml" "$RUNDIR/blocks.nml"
cd "$RUNDIR"
set +e
mpirun -n 8 $MPIRUN_OPTIONS ./snac > run.log 2> run.err
status=$?
set -e
if [ -s run.err ]; then
  cat run.err
  exit 1
fi
if [ $status -ne 0 ]; then
  exit $status
fi
python "$TESTDIR/test.py"
rm -rf "$RUNDIR"
