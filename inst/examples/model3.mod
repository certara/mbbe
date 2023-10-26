$PROBLEM    Simulated BE
$INPUT      ID TIME AMT DROP CMT CMTT DV WT EGFR AGE COV1 COV2 EVID
            PERIOD GROUP TRT SEQ
$DATA       data.csv IGNORE=@ REWIND
$OMEGA
 1  FIX  ;  ETA(1) CL
 1  FIX  ;  ETA(2) V2
$OMEGA  BLOCK(1) FIX
 1  ; ETA(3) BOVCL
$OMEGA  BLOCK(1) SAME
$OMEGA  BLOCK(1) SAME
$OMEGA  BLOCK(1) SAME
$SIGMA  1  FIX  ;     EPS(1)
;;;; Start subs
$SUBROUTINE ADVAN4 ;; advan4
$PK
  CWTKG	= WT/70
  CAGE	= AGE/24
  CEGFR	= EGFR/100
  CCOV1	= COV1/10
  CCOV2 = COV2/10
  IF(PERIOD.EQ.1) BOVCL = EXP(THETA(14))*ETA(3)
  IF(PERIOD.EQ.2) BOVCL = EXP(THETA(14))*ETA(4)
  IF(PERIOD.EQ.3) BOVCL = EXP(THETA(14))*ETA(5)
  IF(PERIOD.EQ.4) BOVCL = EXP(THETA(14))*ETA(6)
  ;; COVARIATE POWER MODEL WITH COVARIATE CENTERED AT 1 CAN BE +IVE OR -IVE
  TVCL	= EXP(THETA(1)) *CWTKG**THETA(15) *CAGE**THETA(16) *CEGFR**THETA(17)
  CL 	= TVCL*EXP(EXP(THETA(4))*ETA(1))*EXP(BOVCL)
  TVV	= EXP(THETA(9)) *CWTKG**THETA(18)   *CCOV2**THETA(19)
  V2 	= TVV *EXP(EXP(THETA(5))*ETA(2))
  IF(TRT.EQ.1) THEN ;; REFERENCE
     TVKA 	= EXP(THETA(6))
     F1		= 1
     ALAG1	= EXP(THETA(12))
  ELSE   ;; TEST
     TVKA 	= EXP(THETA(7))
     F1 	= EXP(THETA(8))
     ALAG1	= EXP(THETA(13))
  END IF
  S2   	= V2/1000 	; CONC IN NG/ML (MCG/L), DOSE IN MG, VOL IN L
  K	= CL/V2
  KA 	= TVKA
  K32 	= EXP(THETA(10))
  K23	= EXP(THETA(11))

$ERROR
  IPRED 	= F
  ADD  		= EXP(THETA(2)) ;; NEEDS TO BE HERE, TO EDIT IN SPPC
  PROP 		= EXP(THETA(3))
  SD 		= SQRT(PROP**2*IPRED**2 + ADD**2) ; Residual weight ADD AND P PROP IN SD AND CV UNITS, NOT VARIANCE

  Y 		= IPRED + SD*EPS(1)

  ;;;; Start EST
$ESTIMATION METHOD=0  SADDLE_RESET=1 MAX=9999 NOABORT PRINT=50
$THETA
 0.17501 		; THETA(1) LN(CL)
 -1.0088 		; THETA(2) ADD ERROR
 -2.99615 		; THETA(3) PROP ERROR
 -2.00413 		; THETA(4) LN(BSV) CL
 -2.08842 		; THETA(5) LN(BSV) V2
 -1.00151 		; THETA(6) LN(KA) reference
 -0.198877 		; THETA(7) LN(KA) test
 -0.0978107 		; THETA(8) LN(F1) TEST
 0.0795407 		; THETA(9) init for V2
 (-5,-3.01036,5) 	; THETA(10) LN(K32)
 (-5,-2.99986,5) 	; THETA(11) LN(K23)
 (-5,-0.500935,5) 	; THETA(12) LN(ALAG) REFERENCE
 (-5,-0.500276,5) 	; THETA(13) LN(ALAG) test
 (-10,-5.14509,5) 	; THETA(14) BOVCL
 0.770386 		; THETA(15) CL~WT
 0.506578 		; THETA(16) CL~AGE
 1.0409 		; THETA(17) CL~EGFR
 0.973552		; THETA(18) V~WT estimated
 0.0506798 		; THETA(19) V~COV2


$TABLE ID TIME GROUP TRT DV IPRED EVID PERIOD SEQ NOPRINT FILE=DATA.TXT NOHEADER NOAPPEND
