#! /bin/bash

rm -f PB_MIC
make -j 8
if [ -f PB_MIC ] ; then
    if [ -z $1 ] ; then 
        time ./PB_MIC tests/mine50.in
    else
        for i in $@ ; do
            item=tests/`basename $i`
            if [ -f $item ] ; then
                time ./PB_MIC tests/`basename $i`
                fi
            done
        fi
    fi
