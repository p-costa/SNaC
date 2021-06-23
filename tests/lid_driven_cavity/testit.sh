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
for FLAGS in ' ' '-D_FFT_X'; do
cp $TESTDIR/Makefile $SRCDIR && cd $SRCDIR && make clean && make OTH=$FLAGS -j run
cp $TESTDIR/dns.in $RUNDIR && cd $RUNDIR
rm -rf $RUNDIR/geo && cp -r $TESTDIR/geo $RUNDIR/geo
echo "Running SNaC..."
sleep 2
mpirun -n 4 --oversubscribe ./${EXEC}
cp $TESTDIR/*.* data/ && cp $UTILSDIR/read_binary_data/python/read_single_field_binary.py $RUNDIR/data/ && cd $RUNDIR/data/
echo "Running test..."
sleep 2
pytest test.py
rm -rf $RUNDIR
done
