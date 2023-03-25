#!/bin/bash

#$ -wd /data/iu-nonmem-bayes-2023/model/nonmem/200/200_4

/opt/NONMEM/nm75/run/nmfe75 200_4.ctl  200_4.lst  -parafile=200_4.pnm -maxlim=2
