$PROBLEM    Simulated BE
$INPUT      ID TIME AMT DROP CMT CMTT DV WT EGFR AGE COV1 COV2 EVID
            PERIOD GROUP TRT SEQ
$DATA       data.csv IGNORE=@ REWIND
$OMEGA
 1  FIX  ;  ETA(1) CL
 1  FIX  ;  ETA(2) V2
$SIGMA  1  FIX  ;     EPS(1)
;;;; Start subs
$SUBROUTINE ADVAN4 ;; advan4
$PK
  ;OCC = PERIOD
  CWTKG	= WT/70
  CAGE	= AGE/24
  CEGFR	= EGFR/100
  CCOV1	= COV1/10
  CCOV2 = COV2/10
  ;; COVARIATE POWER MODEL WITH COVARIATE CENTERED AT 1 CAN BE +IVE OR -IVE
  TVCL	= EXP(THETA(1)) *CWTKG**THETA(14)  *CEGFR**THETA(15) *CCOV1**THETA(16) *CCOV2**THETA(17)
  CL 	= TVCL*EXP(EXP(THETA(4))*ETA(1))
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
$ESTIMATION METHOD=0   MAX=9999 NOABORT PRINT=50
$THETA
 (-0.4,0.221616,2) 		; THETA(1) LN(CL)
 -8.01311 		; THETA(2) ADD ERROR
 -2.97418 		; THETA(3) PROP ERROR
 -2.11214 		; THETA(4) LN(BSV) CL
 -2.08397 		; THETA(5) LN(BSV) V2
 (-2,-0.999677,1) 		; THETA(6) LN(KA) reference
 (-2,-1.00021,1) 		; THETA(7) LN(KA) test
 -0.00145582 	; THETA(8) LN(F1) TEST
 (-0.4,0.0908271,1) 		; THETA(9) init for V2
 (-5,3.99066,5) ; THETA(10) LN(K32)
 (-5,-4.98,5) 	; THETA(11) LN(K23)
 (-2,-0.500032,1) 	; THETA(12) LN(ALAG) REFERENCE
 (-2,-0.501285,1) 	; THETA(13) LN(ALAG) test
 0.00907517 	; THETA(14) CL~WT
 -2.09069 		; THETA(15) CL~EGFR
 0.0265038 		; THETA(16) CL~COV1
 0.167989 		; THETA(17) CL~COV2
 0.0319491 		; THETA(18) V~WT estimated
 -0.11656 		; THETA(19) V~COV2


$TABLE ID TIME GROUP TRT DV IPRED EVID PERIOD SEQ NOPRINT FILE=DATA.TXT NOHEADER NOAPPEND
