#!/bin/bash
TESTDIR=$(pwd)
SNACDIR=$(pwd)/../..
SRCDIR=$SNACDIR/src
RUNDIR=$SNACDIR/run
UTILSDIR=$SNACDIR/utils
EXEC=snac
rm -rf $RUNDIR
echo "Compiling ..."
sleep 2
cp $TESTDIR/Makefile $SRCDIR && cd $SRCDIR && make clean && make -j run
for i in "_x" "_y" "_z"; do
  cp $TESTDIR/dns${i}.in $RUNDIR/dns.in && cd $RUNDIR && mkdir -p $RUNDIR/data
  echo "Running SNaC..."; sleep 2
  mpirun -n 4 ./${EXEC}
  mkdir -p $RUNDIR/data${i} && mv $RUNDIR/data/* $RUNDIR/data${i} && cd $RUNDIR/data${i}
  cp $TESTDIR/*.* $RUNDIR/data${i} && cp $UTILSDIR/read_binary_data/python/read_single_field_binary.py $RUNDIR/data${i} && cd $RUNDIR/data${i}
  echo "Running test..."; sleep 2
  pytest test.py
done
rm -rf $RUNDIR
