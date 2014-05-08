#!/bin/bash
cd ..
tar xvzf nlopt-2.4.1.tar.gz
cd nlopt-2.4.1
./configure --prefix=$PWD --enable-shared
make
make install
cp include/* ../inst/include/
cp -r lib ../inst/lib
ls -R ../inst
cd ../src
