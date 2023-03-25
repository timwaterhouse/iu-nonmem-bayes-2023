$PROBLEM From bbr: see 200.yaml for details

$INPUT C NUM ID TIME CMT EVID MDV DV AMT RATE DOSE LDOS SEX AGE WT EGFR BILI BLQ DAY TAD ADDL II NMTIME STUDY


$DATA ../../../data/derived/study-001-002-003.csv
  IGNORE=@ 
  IGNORE=(BLQ.EQ.1)
  IGNORE=(TIME.LT.0)

$SUBROUTINE ADVAN4 TRANS4

$PK

; cap EGFR at 120
  GFRC = EGFR
  IF(EGFR.GT.120) GFRC = 120

  ; Patient vs healthy
  PAT=0
  IF(STUDY.EQ.3) PAT=1
  ; male vs female
  ISEX=0
  IF(SEX.EQ.2) ISEX=1

  WTCL = LOG(WT/70)*THETA(7)
  WTV = LOG(WT/70)*THETA(8)
  AGECL = LOG(AGE/40)*THETA(9)
  EGFRCL = LOG(GFRC/120)*THETA(10)
  PATCL = THETA(11)*PAT
  SEXCL = THETA(12)*ISEX

  MU_1 = THETA(1) + WTCL + AGECL + EGFRCL + PATCL + SEXCL
  CL = EXP(MU_1 + ETA(1))

  MU_2 = THETA(2) + WTV
  V2 = EXP(MU_2 + ETA(2))

  MU_3 = THETA(3)
  KA = EXP(MU_3 + ETA(3))

  MU_4 = THETA(4)
  D1 = EXP(MU_4 + ETA(4))

  MU_5 = THETA(5) + WTCL
  Q = EXP(MU_5 + ETA(5))

  MU_6 = THETA(6) + WTV
  V3 = EXP(MU_6 + ETA(6))

  S2 = V2

$ERROR
  IPRED = F
  Y = IPRED*(1 + EPS(1)) + EPS(2)

$THETA 
  2        ; CL
  5        ; V2
  0.5      ; KA
  0.1      ; D1
  -1       ; Q
  4        ; V3

  0.75     ; WTCL
  1.00     ; WTV
  0.01     ; AGECL
  0.01     ; EGFRCL
  0.01     ; PATCL
  0.01     ; SEXCL

$OMEGA BLOCK(2)
  0.1      ; CL
  0.05 0.1 ; V2
$OMEGA
  0.1      ; KA
  0 FIX    ; D1
  0 FIX    ; Q
  0 FIX    ; V3

$SIGMA
  0.1      ; PROP
  10       ; ADD

$PRIOR NWPRI

$THETAP 
  2   FIX     ; CL
  5   FIX     ; V2
  0.5 FIX     ; KA
  0.1 FIX     ; D1
  -1  FIX     ; Q
  4   FIX     ; V3

  0.75 FIX  ; WTCL
  1.00 FIX  ; WTV
  0.01 FIX  ; AGECL
  0.01 FIX  ; EGFRCL
  0.01 FIX  ; PATCL
  0.01 FIX  ; SEXCL

$THETAPV BLOCK(12) VALUES(4, 0) FIX

$OMEGAP BLOCK(2) FIX
  0.1      ; CL
  0.05 0.1 ; V2
$OMEGAP
  0.1 FIX  ; KA
  0 FIX    ; D1
  0 FIX    ; Q
  0 FIX    ; V3

$OMEGAPD (3 FIX) (1 FIX) (1 FIX) (1 FIX) (1 FIX)

$SIGMAP
  0.1 FIX     ; PROP
  10  FIX     ; ADD

$SIGMAPD (1 FIX) (1 FIX)

; CHAIN options:
;   CTYPE=0: initial estimates for THETA are sampled from a uniform
;     distribution between (1-IACCEPT)*THETA and (1+IACCEPT)*THETA)
;   CTYPE=2: initial estimates for THETA are from a normal distribution with
;     mean from the initial estimate in $THETA and variance from $THETAPV
;   DF=0: initial estimates for OMEGA come from Wishart distribution using
;     values in $OMEGA and degrees of freedom equal to dimensions of OMEGA
;   DFS=0: initial estimates for SIGMA come from Wishart distribution using
;     values in $SIGMA and degrees of freedom equal to dimensions of SIGMA
$EST METHOD=CHAIN FILE=200.chn NSAMPLE=4 ISAMPLE=0 SEED=1 CTYPE=0 IACCEPT=0.3 DF=10 DFS=3
;$EST METHOD=NUTS AUTO=2 SEED=1 NBURN=500 NITER=1000 PRINT=1 MSFO=./200.msf RANMETHOD=P PARAFPRINT=10000 BAYES_PHI_STORE=1

;$TABLE NUM KA D1 CL V2 Q V3 ETAS(1:LAST) EPRED IPRED NPDE EWRES NOPRINT ONEHEADER FILE=200.tab RANMETHOD=P

