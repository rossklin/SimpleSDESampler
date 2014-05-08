#!/bin/bash
cd ../nlopt-2.4.1
./configure --prefix=$PWD
make
make install
cp include/nlopt.h* ../inst/include/
cp -r lib ../inst/lib
