$PROBLEM    tacrolimis for ML
$INPUT      C ID TOTIME DV_ORG DOSE TRT TREATMENT=DROP TIME MDV AMT OCC
  EVID GROUP   DROP   BLQ DV SEQ
$DATA      U:\fda\mbbe\Modelaveraging\data_seq.csv IGNORE=@ REWIND 
$OMEGA  
  1  FIX  	;  ETA(1) CL
  1  FIX  	;  ETA(2) V2 
  1  FIX  	;  ETA(3) KA REFERENCE
  1  FIX  	;  ETA(4) KA TEST
  1  FIX  	;  ETA(5) ALAG1 REFERENCE
  1  FIX 	;  ETA(6) ALAG1 TEST 
  1  FIX  	;  ETA(7) D1 REFERENCE
  1  FIX 	;  ETA(8) D1 TEST 
$OMEGA  BLOCK(4)
  1  FIX ; ETA(9) BOVV
  0 1    ; ETA(10) BOVCL
  0 0 1    ; ETA(11) BOVKA
  0 0 0 1    ; ETA(12) BOVALAG1
$OMEGA   BLOCK(4) SAME ; ETAs 13,14,15,16
$OMEGA   BLOCK(4) SAME ; ETAs 17,18,19,20
$OMEGA   BLOCK(4) SAME ; ETAs 21,22,23,24
$SIGMA  
  1  FIX  ; EPS(1) PROPORTIONAL
  1  FIX  ; EPS(2) ADDITIVE
  ;;;; Start subs
$SUBROUTINE ADVAN4 ;; advan4
$PK

  IF(GROUP.EQ.1) BOVV = THETA(16)*ETA(9)
  IF(GROUP.EQ.2) BOVV = THETA(16)*ETA(13)
  IF(GROUP.EQ.3) BOVV = THETA(16)*ETA(17)
  IF(GROUP.EQ.4) BOVV = THETA(16)*ETA(21)

  IF(GROUP.EQ.1) BOVCL = THETA(17)*ETA(10)
  IF(GROUP.EQ.2) BOVCL = THETA(17)*ETA(14)
  IF(GROUP.EQ.3) BOVCL = THETA(17)*ETA(18)
  IF(GROUP.EQ.4) BOVCL = THETA(17)*ETA(22)

  BOVKA=0 
  BOVALAG1=0 

  ;; EMAIL FROM MICHAEL/SHUHUA 27FEB2023
  ;;CL  = EXP(THETA(1) + THETA(3) * ETA(1) + THETA(34) * ETA(2))
  ;;V2  = EXP(THETA(2) + THETA(34) * ETA(1) + THETA(4) * ETA(2))

  CL 	= EXP(THETA(1)+THETA(3)*ETA(1)    + BOVCL) ;; ETACLV IS CORRELATION OF ETA(V) AND ETA(CL)
  V2 	= EXP(THETA(2)+THETA(4)*ETA(2)    + BOVV) 
  IF(TRT.EQ.1) THEN ;; REFERENCE
  KA 	= EXP(THETA(5) + BOVKA  )  
  F1	= 1
  ALAG1	= EXP(THETA(12)+THETA(14)*ETA(5))
  ; NO D1
  ELSE   ;; TEST
  KA 	= EXP(THETA(6) + BOVKA )  
  F1 	= EXP(THETA(7))
  ALAG1	= EXP(THETA(13)+THETA(15)*ETA(6))
  ; NO D1
  END IF   

  S2   	= V2/1000000 	; CONC IN PG/M (NG/L), DOSE IN MG, VOL IN L 
  K=CL/V2  

  K32 	= EXP(THETA(10))
  K23	= EXP(THETA(11))

$ERROR   
  IPRED = F  
  REALOBS = DV ;;; TO PRINT OUT FOR ORG.DAT
  ADD = EXP(THETA(9))
  PROP = EXP(THETA(8))
  LLOQ = 50 
  SD = SQRT(PROP**2*IPRED**2 + ADD**2) ; Residual weight ADD AND P PROP IN SD AND CV UNITS, NOT VARIANCE
  IF (BLQ.EQ.0) THEN
  F_FLAG=0 ; ELS
  Y = IPRED + SD*EPS(1)                          ; Individual model prediction,
  ENDIF

  IF (BLQ.EQ.1) THEN
  F_FLAG=1 ; LIKELIHOOD
  Y=PHI((LLOQ-IPRED)/SD)
  ENDIF

  ;;;; Start EST

$ESTIMATION METHOD=1 LAPLAC INTER MAX=9999 PRINT=10 NOABORT

$COVARIANCE PRINT=E UNCOND PRECOND=2

$THETA  
  (-1,3.8,10) 		; THETA(1) LN(CL1)
  (-1,5.2,12) 		; THETA(2) LN(V2)
  (-4,-0.61) 		; THETA(3) LN(ETA(CL))
  (-4,-0.38)		; THETA(4) LN(ETA(V))
  (-5,-0.4,5)	 	; THETA(5) LN(KA) REFERENCE
  (-5,0.6,5)	 	; THETA(6) LN(KA) TEST 
  (-3,0.15,1)	 	; THETA(7) LN(F) TEST  
  (-6,-0.8,2)		; THETA(8) LN(PROPERROR)
  (-7,1,3)		; THETA(9) LN(ADDERROR)
  (-5,-2.53064,5) 		; THETA(10) LN(K32)
  (-5,-0.89632,5) 		; THETA(11) LN(K23)
  (-5,-2,5)		; THETA(12) LN(ALAG) REFERENCE 
  (-5,-2,5) 		; THETA(13) LN(ALAG) test
  (-5,-1,5)		; THETA(14) LN(ALAGETA) REFERENCE 
  (-5,-1,5) 		; THETA(15) LN(ALAGETA) test
  ; NO D1
  ;;no eta on ka
  (-4,-0.3,5)	; THETA(16) BOVV 
  (-4,-0.3,5)	; THETA(17) BOVCL 
  ;; NO BOVKA 
  ;; NO BOVALAG1  
  ;; NO COVARIANCE BETWEEN V AND CL  


;; Phenotype 
;; OrderedDict([('ADVAN', 2), ('KAETA', 0), ('ALAG1', 3), ('BOVALAG1', 0), ('DURATION1', 0), ('BOVV', 1), ('BOVCL', 1), ('BOVKA', 0), ('COVARVCL', 1)])
;; Genotype 
;; [1, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0, 1]
;; Num non-influential tokens = 0