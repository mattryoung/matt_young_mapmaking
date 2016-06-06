#!/bin/bash

for f in *.f90
do
    g=../make_pkvd/$f
    if [ -e $g ] ; then
	if cmp $f $g >/dev/null 2>&1
	then
	   rm $f
	fi
    fi
done
