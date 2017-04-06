#!/bin/bash

g++ -o example1 example1.cpp -I/home/jjzhu/src/cnpy/include -lcnpy -L/home/jjzhu/src/cnpy/lib

TDIR=/home/jjzhu/src
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$TDIR/boost_1_58_0/lib/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$TDIR/cnpy/lib/

./example1

