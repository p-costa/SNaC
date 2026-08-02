#!/bin/bash
set -e

TESTDIR=$(pwd)
SNACDIR=$(pwd)/../..
RUNDIR=$SNACDIR/run
EXEC=snac
MPIRUN_OPTIONS='--oversubscribe'
if mpirun --version 2>&1 | grep -qi intel; then MPIRUN_OPTIONS=""; fi

run_case() {
  DNS_FILE=$1
  CASE_NAME=$2
  BLOCKS_FILE=${3:-blocks.nml}
  NP=${4:-1}
  cp "$TESTDIR/$DNS_FILE" "$RUNDIR/dns.nml"
  cp "$TESTDIR/$BLOCKS_FILE" "$RUNDIR/blocks.nml"
  rm -rf "$RUNDIR/data"
  mkdir -p "$RUNDIR/data"
  cd "$RUNDIR"
  set +e
  mpirun -n "$NP" $MPIRUN_OPTIONS ./${EXEC} > "$CASE_NAME.log" 2> "$CASE_NAME.err"
  status=$?
  set -e
  if [ -s "$CASE_NAME.err" ]; then
    cat "$CASE_NAME.err"
    exit 1
  fi
  if [ $status -ne 0 ]; then
    exit $status
  fi
  python "$TESTDIR/test.py" "$CASE_NAME.log"
  cd "$TESTDIR"
}

rm -rf "$RUNDIR"
cd "$SNACDIR"
make clean
make FFT_AXIS=0 -j run
run_case dns_bicgstab_pfmg.nml bicgstab_pfmg
run_case dns_flexgmres_pfmg.nml flexgmres_pfmg
run_case dns_flexgmres_smg.nml flexgmres_smg
run_case dns_gmres_none.nml gmres_none

cd "$SNACDIR"
make clean
make FFT_AXIS=1 -j run
run_case dns_fft_hybrid.nml fft_hybrid

cd "$SNACDIR"
make clean
make FFT_AXIS=1 FFT_USE_SLABS=1 -j run
run_case dns_fft_hybrid.nml fft_hybrid_slab blocks_slab.nml 2

rm -rf "$RUNDIR"
