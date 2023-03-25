#!/bin/bash

#$ -wd /data/iu-nonmem-bayes-2023/model/nonmem/200/200_1

/opt/NONMEM/nm75/run/nmfe75 200_1.ctl  200_1.lst  -parafile=200_1.pnm -maxlim=2
