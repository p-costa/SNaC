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
cp $TESTDIR/dns.in $RUNDIR && cd $RUNDIR
echo "Running SNaC..."
sleep 2
mpirun -n 4 ./${EXEC}
cp $TESTDIR/*.* data/ && cp $UTILSDIR/read_binary_data/python/read_single_field_binary.py $RUNDIR/data/ && cd $RUNDIR/data/
echo "Running test..."
sleep 2
pytest test.py
rm -rf $RUNDIR
