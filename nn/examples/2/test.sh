#!/bin/bash

RES=256x256

if (($# > 0))
then
    RES=$1
fi

if [ ! -x ../../nnbathy ]
then
    echo "error: no executable found"
    echo 'Run "./configure" and "make" in the source directory'
    exit 1
fi

echo ""
echo -n "Linear interpolation..."
../../nnbathy -i data.txt -n $RES -P alg=l > lin.txt
echo "done"
echo -n "Natural Neighbours Sibson interpolation..."
../../nnbathy -i data.txt -n $RES > nn-inf.txt
echo "done"
echo -n "Natural Neighbours Sibson interpolation with wmin = 0..."
../../nnbathy -i data.txt -n $RES -W 0 > nn-0.txt
echo "done"
echo -n "Natural Neighbours Non-Sibsonian interpolation with wmin = 0..."
../../nnbathy -i data.txt -n $RES -W 0 -P alg=ns > nn-ns.txt
echo "done"
echo ""
echo 'To visualize, in Matlab run "viewexample"'
