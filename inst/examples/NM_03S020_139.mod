$PROBLEM    Simulated BE- Passes
$INPUT      ID TIME AMT  DROP   CMT CMTT DV WT EGFR AGE COV1 COV2 EVID PERIOD GROUP TRT SEQ

$DATA u:\fda\mbbe\inst\examples\DATA.csv IGNORE=@ REWIND
$OMEGA  
  1  FIX  	;  ETA(1) CL
  1  FIX  	;  ETA(2) V2 
  1  FIX  	;  ETA(3) KA 
  1  FIX  	;  ETA(4) ALAG1 
  1  FIX  	;  ETA(5) D1 
$OMEGA  BLOCK(2)
  1  FIX ; ETA(6) BOVV
  0 1    ; ETA(7) BOVCL 
$OMEGA   BLOCK(2) SAME ; ETAs 8,9
$OMEGA   BLOCK(2) SAME ; ETAs 10,11
$OMEGA   BLOCK(2) SAME ; ETAs 12,13
$SIGMA  
  1  FIX  ; EPS(1) PROPORTIONAL
  1  FIX  ; EPS(2) ADDITIVE
  ;;;; Start subs
$SUBROUTINE ADVAN4 ;; advan4  
$PK
  OCC = PERIOD
  CWTKG	= WT/70
  CAGE	= AGE/24
  CEGFR	= EGFR/100
  CCOV1	= COV1/10
  CCOV2 = COV2/10 
  BOVV = 0 
  BOVCL = 0   
  ;; COVARIATE POWER MODEL WITH COVARIATE CENTERED AT 1 CAN BE +IVE OR -IVE
  TVCL= EXP(THETA(1)) *CWTKG**THETA(14)  *CEGFR**THETA(15) *CCOV1**THETA(16) *CCOV2**THETA(17)
  CL 	= TVCL*EXP(EXP(THETA(4))*ETA(1))*EXP(BOVCL) *EXP(0) ;; ETACLV IS CORRELATION OF ETA(V) AND ETA(CL)
  TVV	= EXP(THETA(9)) *CWTKG**THETA(18)   *CCOV2**THETA(19)
  V2 	= TVV *EXP(EXP(THETA(5))*ETA(2))*EXP(BOVV) *EXP(0)
  IF(TRT.EQ.1) THEN ;; REFERENCE
  TVKA 	= EXP(THETA(6))  
  F1	= 1
  ALAG1	= EXP(THETA(12))
  ; NO D1
  ELSE   ;; TEST
  TVKA 	= EXP(THETA(7))     
  F1 	= EXP(THETA(8))
  ALAG1	= EXP(THETA(13))
  ; NO D1
  END IF   
  S2   	= V2/1000 	; CONC IN NG/ML (MCG/L), DOSE IN MG, VOL IN L 
  K=CL/V2  
  KA = TVKA *EXP(0)  
  K32 	= EXP(THETA(10))
  K23	= EXP(THETA(11))

$ERROR  
  IPRED 	= F  
  REALOBS 	= DV ;;; TO PRINT OUT FOR ORG.DAT
  ADD  		= EXP(THETA(2)) ;; NEEDS TO BE HERE, TO EDIT IN SPPC
  PROP 		= EXP(THETA(3))
  SD 		= SQRT(PROP**2*IPRED**2 + ADD**2) ; Residual weight ADD AND P PROP IN SD AND CV UNITS, NOT VARIANCE
  OBSWERROR	= IPRED + SD*EPS(1) ;; observation with error

  Y 	= IPRED + SD*EPS(1)    

  ;;;; Start EST 
$EST METHOD=COND INTER  SADDLE_RESET=1 MAX=9999 NOABORT PRINT = 50
$COV PRINT=E UNCOND  PRECOND=2 MATRIX=S ;; MATRIX S FOR SPEED

$THETA 	; all are log of THETA
  (0.2)	; THETA(1) LN(CL) 
  (-0.5) 	; THETA(2) ADD ERROR
  (-3) 	; THETA(3) PROP ERROR
  (-2)  	; THETA(4) LN(BSV) CL
  (-2)  	; THETA(5) LN(BSV) V2
  (-1)	; THETA(6) LN(KA) reference
  (-0.3)    ; THETA(7) LN(KA) test
  (-0.1)  ; THETA(8) LN(F1) TEST  

  (0.1) 		; THETA(9) init for V2
  (-5,-3,5) 		; THETA(10) LN(K32)
  (-5,-3,5) 		; THETA(11) LN(K23)
  (-5,-0.5,5)		; THETA(12) LN(ALAG) REFERENCE 
  (-5,-0.5,5) 		; THETA(13) LN(ALAG) test
  ; NO D1
  ;;no eta on ka
  ;; NO BOVV 
  ;; NO BOVCL  
  ;; NO COVARIANCE BETWEEN V AND CL  
  ;; covariate parameters need not be log transformed, they can be < 0

  0.75		; THETA(14) CL~WT 

  0.8		; THETA(15) CL~EGFR 
  -0.1		; THETA(16) CL~COV1 
  -0.1		; THETA(17) CL~COV2

  1		; THETA(18) V~WT estimated 


  -0.1		; THETA(19) V~COV2


$TABLE ID TIME GROUP TRT DV IPRED EVID PERIOD SEQ NOPRINT FILE=DATA.TXT NOHEADER NOAPPEND 

;; Phenotype 
;; OrderedDict([('ADVAN', 1), ('BOVV', 0), ('BOVCL', 0), ('KAETA', 0), ('ALAG1', 1), ('DURATION1', 0), ('COVARVCL', 1), ('CL~WT', 1), ('CL~AGE', 0), ('CL~EGFR', 1), ('CL~COV1', 1), ('CL~COV2', 1), ('V~WT', 1), ('V~AGE', 0), ('V~COV1', 0), ('V~COV2', 1)])
;; Genotype 
;; [1, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 0, 1]
;; Num non-influential tokens = 0