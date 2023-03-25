#!/bin/bash

#$ -wd /data/iu-nonmem-bayes-2023/model/nonmem/200

/opt/NONMEM/nm75/run/nmfe75 200.ctl  200.lst  -maxlim=2
