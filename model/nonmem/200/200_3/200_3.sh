#!/bin/bash

#$ -wd /data/iu-nonmem-bayes-2023/model/nonmem/200/200_3

/opt/NONMEM/nm75/run/nmfe75 200_3.ctl  200_3.lst  -parafile=200_3.pnm -maxlim=2
