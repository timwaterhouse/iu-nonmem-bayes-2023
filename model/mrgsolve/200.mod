$PLUGIN autodec nm-vars
$GLOBAL
#define CP (PLASMA/V2)
  
$PARAM  WT = 75
        AGE = 50
        SEX = 2
        EGFR = 120
        STUDY = 1

$NMEXT
run = "200_1"
project = "../nonmem/200"
root = "cppfile"

$CMT DEPOT PLASMA PERIPH

$PKMODEL ncmt=2, depot=TRUE

$PK

// cap EGFR at 120
  GFRC = EGFR;
  if(EGFR > 120) GFRC = 120;

  // Patient vs healthy
  PAT=0;
  if(STUDY == 3) PAT=1;
  // male vs female
  ISEX=0;
  if(SEX == 2) ISEX=1;

  WTCL = LOG(WT/70)*THETA(7);
  WTV = LOG(WT/70)*THETA(8);
  AGECL = LOG(AGE/40)*THETA(9);
  EGFRCL = LOG(GFRC/120)*THETA(10);
  PATCL = THETA(11)*PAT;
  SEXCL = THETA(12)*ISEX;

  MU_1 = THETA(1) + WTCL + AGECL + EGFRCL + PATCL + SEXCL;
  CL = EXP(MU_1 + ETA(1));

  MU_2 = THETA(2) + WTV;
  V2 = EXP(MU_2 + ETA(2));

  MU_3 = THETA(3);
  KA = EXP(MU_3 + ETA(3));

  MU_4 = THETA(4);
  D1 = EXP(MU_4 + ETA(4));

  MU_5 = THETA(5) + WTCL;
  Q = EXP(MU_5 + ETA(5));

  MU_6 = THETA(6) + WTV;
  V3 = EXP(MU_6 + ETA(6));

  S2 = V2;

$ERROR
  IPRED = CP;
  Y = IPRED*(1 + EPS(1)) + EPS(2);

$CAPTURE IPRED Y
