#!/bin/bash

#$ -wd /data/iu-nonmem-bayes-2023/model/nonmem/200/200_2

/opt/NONMEM/nm75/run/nmfe75 200_2.ctl  200_2.lst  -parafile=200_2.pnm -maxlim=2
