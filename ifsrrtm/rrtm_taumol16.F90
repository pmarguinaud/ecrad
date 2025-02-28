!----------------------------------------------------------------------------
SUBROUTINE RRTM_TAUMOL16 (KIDIA,KFDIA,KLEV,P_TAU,&
 & P_TAUAERL,P_FAC00,P_FAC01,P_FAC10,P_FAC11,P_FORFAC,P_FORFRAC,K_INDFOR,K_JP,K_JT,K_JT1,P_ONEMINUS,&
 & P_COLH2O,P_COLCH4,K_LAYTROP,P_SELFFAC,P_SELFFRAC,K_INDSELF,PFRAC, &
 & P_RAT_H2OCH4,P_RAT_H2OCH4_1)  

!     BAND 16:  2600-3000 cm-1 (low - H2O,CH4; high - nothing)

!     AUTHOR.
!     -------
!      JJMorcrette, ECMWF

!     MODIFICATIONS.
!     --------------
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      NEC           25-Oct-2007 Optimisations
!      JJMorcrette 20110613 flexible number of g-points
!      ABozzo 201306 updated to rrtmg v4.85
!      band 16:  2600-3250 cm-1 (low key- h2o,ch4; high key - ch4)
! ---------------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK, DR_HOOK, JPHOOK

USE PARRRTM  , ONLY : JPBAND
USE YOERRTM  , ONLY : JPGPT  ,NGS15  ,NG16
USE YOERRTWN , ONLY : NSPA,NSPB   
USE YOERRTA16, ONLY : ABSA,ABSB,FRACREFA,FRACREFB,SELFREF,FORREF 
USE YOERRTRF, ONLY : CHI_MLS

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: P_TAU(KIDIA:KFDIA,JPGPT,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_TAUAERL(KIDIA:KFDIA,KLEV,JPBAND) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_FAC00(KIDIA:KFDIA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_FAC01(KIDIA:KFDIA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_FAC10(KIDIA:KFDIA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_FAC11(KIDIA:KFDIA,KLEV) 
INTEGER(KIND=JPIM),INTENT(IN)    :: K_JP(KIDIA:KFDIA,KLEV) 
INTEGER(KIND=JPIM),INTENT(IN)    :: K_JT(KIDIA:KFDIA,KLEV) 
INTEGER(KIND=JPIM),INTENT(IN)    :: K_JT1(KIDIA:KFDIA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_ONEMINUS
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_COLH2O(KIDIA:KFDIA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_COLCH4(KIDIA:KFDIA,KLEV) 
INTEGER(KIND=JPIM),INTENT(IN)    :: K_LAYTROP(KIDIA:KFDIA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_SELFFAC(KIDIA:KFDIA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_SELFFRAC(KIDIA:KFDIA,KLEV) 
INTEGER(KIND=JPIM),INTENT(IN)    :: K_INDSELF(KIDIA:KFDIA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFRAC(KIDIA:KFDIA,JPGPT,KLEV) 

REAL(KIND=JPRB)   ,INTENT(IN)   :: P_RAT_H2OCH4(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)   :: P_RAT_H2OCH4_1(KIDIA:KFDIA,KLEV)
INTEGER(KIND=JPIM),INTENT(IN)   :: K_INDFOR(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)   :: P_FORFAC(KIDIA:KFDIA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)   :: P_FORFRAC(KIDIA:KFDIA,KLEV) 

! ---------------------------------------------------------------------------

INTEGER(KIND=JPIM) :: IG, IND0, IND1, INDS,INDF, JS,JS1,JPL,JLAY
INTEGER(KIND=JPIM) :: JLON

REAL(KIND=JPRB) ::  Z_FAC000, Z_FAC100, Z_FAC200,&
 & Z_FAC010, Z_FAC110, Z_FAC210, &
 & Z_FAC001, Z_FAC101, Z_FAC201, &
 & Z_FAC011, Z_FAC111, Z_FAC211
REAL(KIND=JPRB) :: ZP, ZP4, ZFK0, ZFK1, ZFK2

REAL(KIND=JPRB) :: ZREFRAT_PLANCK_A
REAL(KIND=JPRB) :: ZTAUFOR,ZTAUSELF,ZTAU_MAJOR,ZTAU_MAJOR1
REAL(KIND=JPRB) :: Z_FS, Z_SPECMULT, Z_SPECPARM,Z_SPECCOMB,  &
& Z_FS1, Z_SPECMULT1, Z_SPECPARM1,Z_SPECCOMB1, &
& Z_FPL, Z_SPECMULT_PLANCK, Z_SPECPARM_PLANCK,Z_SPECCOMB_PLANCK
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('RRTM_TAUMOL16',0,ZHOOK_HANDLE)

! Calculate reference ratio to be used in calculation of Planck
! fraction in lower atmosphere.

! P = 387. mb (Level 6)
      ZREFRAT_PLANCK_A = CHI_MLS(1,6)/CHI_MLS(6,6)

! Compute the optical depth by interpolating in ln(pressure), 
! temperature,and appropriate species.  Below laytrop, the water
! vapor self-continuum and foreign continuum is interpolated 
! (in temperature) separately.  

DO JLAY = 1, KLEV
  DO JLON = KIDIA, KFDIA
    IF (JLAY <= K_LAYTROP(JLON)) THEN
      Z_SPECCOMB = P_COLH2O(JLON,JLAY) + P_RAT_H2OCH4(JLON,JLAY)*P_COLCH4(JLON,JLAY)
      Z_SPECPARM = P_COLH2O(JLON,JLAY)/Z_SPECCOMB
      Z_SPECPARM = MIN(Z_SPECPARM,P_ONEMINUS)
      Z_SPECMULT = 8._JPRB*(Z_SPECPARM)
      JS = 1 + INT(Z_SPECMULT)
      Z_FS = MOD(Z_SPECMULT,1.0_JPRB)

      Z_SPECCOMB1 = P_COLH2O(JLON,JLAY) + P_RAT_H2OCH4_1(JLON,JLAY)*P_COLCH4(JLON,JLAY)
      Z_SPECPARM1 = P_COLH2O(JLON,JLAY)/Z_SPECCOMB1
      IF (Z_SPECPARM1 >= P_ONEMINUS) Z_SPECPARM1 = P_ONEMINUS
      Z_SPECMULT1 = 8._JPRB*(Z_SPECPARM1)
      JS1 = 1 + INT(Z_SPECMULT1)
      Z_FS1 = MOD(Z_SPECMULT1,1.0_JPRB)

      Z_SPECCOMB_PLANCK = P_COLH2O(JLON,JLAY)+ZREFRAT_PLANCK_A*P_COLCH4(JLON,JLAY)
      Z_SPECPARM_PLANCK = P_COLH2O(JLON,JLAY)/Z_SPECCOMB_PLANCK
      IF (Z_SPECPARM_PLANCK >= P_ONEMINUS) Z_SPECPARM_PLANCK=P_ONEMINUS
      Z_SPECMULT_PLANCK = 8._JPRB*Z_SPECPARM_PLANCK
      JPL= 1 + INT(Z_SPECMULT_PLANCK)
      Z_FPL = MOD(Z_SPECMULT_PLANCK,1.0_JPRB)

      IND0 = ((K_JP(JLON,JLAY)-1)*5+(K_JT(JLON,JLAY)-1))*NSPA(16) + JS
      IND1 = (K_JP(JLON,JLAY)*5+(K_JT1(JLON,JLAY)-1))*NSPA(16) + JS1
      INDS = K_INDSELF(JLON,JLAY)
      INDF = K_INDFOR(JLON,JLAY)

IF (Z_SPECPARM < 0.125_JPRB) THEN
            ZP = Z_FS - 1
            ZP4 = ZP**4
            ZFK0 = ZP4
            ZFK1 = 1 - ZP - 2.0_JPRB*ZP4
            ZFK2 = ZP + ZP4
            Z_FAC000 = ZFK0*P_FAC00(JLON,JLAY)
            Z_FAC100 = ZFK1*P_FAC00(JLON,JLAY)
            Z_FAC200 = ZFK2*P_FAC00(JLON,JLAY)
            Z_FAC010 = ZFK0*P_FAC10(JLON,JLAY)
            Z_FAC110 = ZFK1*P_FAC10(JLON,JLAY)
            Z_FAC210 = ZFK2*P_FAC10(JLON,JLAY)
      ELSEIF (Z_SPECPARM > 0.875_JPRB) THEN
            ZP = -Z_FS 
            ZP4 = ZP**4
            ZFK0 = ZP4
            ZFK1 = 1 - ZP - 2.0_JPRB*ZP4
            ZFK2 = ZP + ZP4
            Z_FAC000 = ZFK0*P_FAC00(JLON,JLAY)
            Z_FAC100 = ZFK1*P_FAC00(JLON,JLAY)
            Z_FAC200 = ZFK2*P_FAC00(JLON,JLAY)
            Z_FAC010 = ZFK0*P_FAC10(JLON,JLAY)
            Z_FAC110 = ZFK1*P_FAC10(JLON,JLAY)
            Z_FAC210 = ZFK2*P_FAC10(JLON,JLAY)
      ELSE
            Z_FAC000 = (1._JPRB - Z_FS) * P_FAC00(JLON,JLAY)
            Z_FAC010 = (1._JPRB - Z_FS) * P_FAC10(JLON,JLAY)
            Z_FAC100 = Z_FS * P_FAC00(JLON,JLAY)
            Z_FAC110 = Z_FS * P_FAC10(JLON,JLAY)
      ENDIF
      IF (Z_SPECPARM1 < 0.125_JPRB) THEN
            ZP = Z_FS1 - 1
            ZP4 = ZP**4
            ZFK0 = ZP4
            ZFK1 = 1 - ZP - 2.0_JPRB*ZP4
            ZFK2 = ZP + ZP4
            Z_FAC001 = ZFK0*P_FAC01(JLON,JLAY)
            Z_FAC101 = ZFK1*P_FAC01(JLON,JLAY)
            Z_FAC201 = ZFK2*P_FAC01(JLON,JLAY)
            Z_FAC011 = ZFK0*P_FAC11(JLON,JLAY)
            Z_FAC111 = ZFK1*P_FAC11(JLON,JLAY)
            Z_FAC211 = ZFK2*P_FAC11(JLON,JLAY)
      ELSEIF (Z_SPECPARM1 > 0.875_JPRB) THEN
            ZP = -Z_FS1 
            ZP4 = ZP**4
            ZFK0 = ZP4
            ZFK1 = 1 - ZP - 2.0_JPRB*ZP4
            ZFK2 = ZP + ZP4
            Z_FAC001 = ZFK0*P_FAC01(JLON,JLAY)
            Z_FAC101 = ZFK1*P_FAC01(JLON,JLAY)
            Z_FAC201 = ZFK2*P_FAC01(JLON,JLAY)
            Z_FAC011 = ZFK0*P_FAC11(JLON,JLAY)
            Z_FAC111 = ZFK1*P_FAC11(JLON,JLAY)
            Z_FAC211 = ZFK2*P_FAC11(JLON,JLAY)
      ELSE
            Z_FAC001 = (1._JPRB - Z_FS1) * P_FAC01(JLON,JLAY)
            Z_FAC011 = (1._JPRB - Z_FS1) * P_FAC11(JLON,JLAY)
            Z_FAC101 = Z_FS1 * P_FAC01(JLON,JLAY)
            Z_FAC111 = Z_FS1 * P_FAC11(JLON,JLAY)
      ENDIF


      DO IG = 1, NG16
!- - DS_990714
        ZTAUSELF = P_SELFFAC(JLON,JLAY)* (SELFREF(INDS,IG) + P_SELFFRAC(JLON,JLAY) * &
          &       (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
        ZTAUFOR = P_FORFAC(JLON,JLAY) * (FORREF(INDF,IG) + P_FORFRAC(JLON,JLAY) * &
          &       (FORREF(INDF+1,IG) - FORREF(INDF,IG))) 

     IF (Z_SPECPARM < 0.125_JPRB) THEN
               ZTAU_MAJOR = Z_SPECCOMB * &
                 &   (Z_FAC000 * ABSA(IND0,IG) + &
                 &   Z_FAC100 * ABSA(IND0+1,IG) + &
                 &   Z_FAC200 * ABSA(IND0+2,IG) + &
                 &   Z_FAC010 * ABSA(IND0+9,IG) + &
                 &   Z_FAC110 * ABSA(IND0+10,IG) + &
                 &   Z_FAC210 * ABSA(IND0+11,IG))
            ELSEIF (Z_SPECPARM > 0.875_JPRB) THEN
               ZTAU_MAJOR = Z_SPECCOMB * &
                 &   (Z_FAC200 * ABSA(IND0-1,IG) + &
                 &   Z_FAC100 * ABSA(IND0,IG) + &
                 &   Z_FAC000 * ABSA(IND0+1,IG) + &
                 &   Z_FAC210 * ABSA(IND0+8,IG) + &
                 &   Z_FAC110 * ABSA(IND0+9,IG) + &
                 &   Z_FAC010 * ABSA(IND0+10,IG))
            ELSE
               ZTAU_MAJOR = Z_SPECCOMB * &
                 &   (Z_FAC000 * ABSA(IND0,IG) + &
                 &   Z_FAC100 * ABSA(IND0+1,IG) + &
                 &   Z_FAC010 * ABSA(IND0+9,IG) + &
                 &   Z_FAC110 * ABSA(IND0+10,IG))
            ENDIF

            IF (Z_SPECPARM1 < 0.125_JPRB) THEN
               ZTAU_MAJOR1 = Z_SPECCOMB1 * &
                &    (Z_FAC001 * ABSA(IND1,IG) + &
                &    Z_FAC101 * ABSA(IND1+1,IG) + &
                &    Z_FAC201 * ABSA(IND1+2,IG) + &
                &    Z_FAC011 * ABSA(IND1+9,IG) + &
                &    Z_FAC111 * ABSA(IND1+10,IG) + &
                &    Z_FAC211 * ABSA(IND1+11,IG))
            ELSEIF (Z_SPECPARM1 > 0.875_JPRB) THEN
               ZTAU_MAJOR1 = Z_SPECCOMB1 * &
                &    (Z_FAC201 * ABSA(IND1-1,IG) + &
                &    Z_FAC101 * ABSA(IND1,IG) + &
                &    Z_FAC001 * ABSA(IND1+1,IG) + &
                &    Z_FAC211 * ABSA(IND1+8,IG) + &
                &    Z_FAC111 * ABSA(IND1+9,IG) + &
                &    Z_FAC011 * ABSA(IND1+10,IG))
            ELSE
               ZTAU_MAJOR1 = Z_SPECCOMB1 * &
                &    (Z_FAC001 * ABSA(IND1,IG) +  &
                &    Z_FAC101 * ABSA(IND1+1,IG) + &
                &    Z_FAC011 * ABSA(IND1+9,IG) + &
                &    Z_FAC111 * ABSA(IND1+10,IG))
            ENDIF


        P_TAU(JLON,NGS15+IG,JLAY) = ZTAU_MAJOR + ZTAU_MAJOR1 &
                & + ZTAUSELF + ZTAUFOR &
                & + P_TAUAERL(JLON,JLAY,16)  
        PFRAC(JLON,NGS15+IG,JLAY) = FRACREFA(IG,JPL) + Z_FPL * &
         & (FRACREFA(IG,JPL+1) - FRACREFA(IG,JPL))

        ENDDO

!-- DS_990714
    ENDIF

    IF (JLAY > K_LAYTROP(JLON)) THEN
      IND0 = ((K_JP(JLON,JLAY)-13)*5+(K_JT(JLON,JLAY)-1))*NSPB(16) + 1
      IND1 = ((K_JP(JLON,JLAY)-12)*5+(K_JT1(JLON,JLAY)-1))*NSPB(16) + 1
      DO IG = 1, NG16

        P_TAU(JLON,NGS15+IG,JLAY) =  P_COLCH4(JLON,JLAY) * &
               &  (P_FAC00(JLON,JLAY) * ABSB(IND0,IG) + &
               &  P_FAC10(JLON,JLAY) * ABSB(IND0+1,IG) + &
               &  P_FAC01(JLON,JLAY) * ABSB(IND1,IG) + &
               &  P_FAC11(JLON,JLAY) * ABSB(IND1+1,IG)) + &
               &  P_TAUAERL(JLON,JLAY,16)
        PFRAC(JLON,NGS15+IG,JLAY) = FRACREFB(IG)

      ENDDO
    ENDIF
  ENDDO
ENDDO

IF (LHOOK) CALL DR_HOOK('RRTM_TAUMOL16',1,ZHOOK_HANDLE)

END SUBROUTINE RRTM_TAUMOL16
