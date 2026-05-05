#!/bin/bash
set -e

TESTDIR=$(pwd)
SNACDIR=$(pwd)/../..
SRCDIR=$SNACDIR/src
RUNDIR=$SNACDIR/run
UTILSDIR=$SNACDIR/utils
EXEC=snac
MPIRUN_OPTIONS='--oversubscribe'
if mpirun --version 2>&1 | grep -qi intel; then MPIRUN_OPTIONS=""; fi
MPIRUN="mpirun -n 4 $MPIRUN_OPTIONS"

rm -rf $RUNDIR
echo "Compiling ..."
for FLAGS in ' ' '-D_FFT_X'; do
  cd $SRCDIR && make -f $TESTDIR/Makefile clean && make -f $TESTDIR/Makefile OTH="$FLAGS" -j run
  cp $TESTDIR/dns.nml $TESTDIR/blocks.nml $RUNDIR && cd $RUNDIR
  echo "Running SNaC..."
  set +e
  ${MPIRUN} ./${EXEC} 1> log_file.log 2> err_log.log
  status=$?
  set -e
  (head -n 50 log_file.log; echo ""; echo "[...output omitted...]"; echo ""; tail -n 50 log_file.log)
  if [ -s err_log.log ]; then
    echo ""
    echo "stderr from SNaC"
    cat err_log.log
  fi
  if [ $status -ne 0 ]; then
    exit $status
  fi
  cp $TESTDIR/*.* data/
  cp $UTILSDIR/read_binary_data/python/read_single_field_binary.py $RUNDIR/data/
  cd $RUNDIR/data
  echo "Running test..."
  python test.py
  rm -rf $RUNDIR
done
