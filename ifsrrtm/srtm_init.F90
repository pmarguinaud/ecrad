SUBROUTINE SRTM_INIT_GPU(CDIRECTORY, NWVCONTINUUM, lacc)

!-- read in the basic coefficients to configure RRTM_SW
!- creates module YOESRTWN with BG, NSPA, NSPB, WAVENUM1, WAVENUM2,
!  DELWAVE, PREF, PREFLOG, TREF

USE PARKIND1
USE YOMHOOK

USE PARSRTM
USE YOESRTM  , ONLY : NGN
USE YOESRTWN , ONLY : NG, NGM, WT, NGC, RWGT, WTSM
!USE YOESRTWN , ONLY : NG, NGM, WT, NGC, NGN, RWGT, WTSM
!USE YOMLUN   , ONLY : NULOUT

IMPLICIT NONE

CHARACTER(LEN=*), INTENT(IN) :: CDIRECTORY

! Water vapour continuum model (0=default MT_CKD2.5, 1=CAVIAR)
INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL :: NWVCONTINUUM

! Local variables
INTEGER(KIND=JPIM) :: IGC, IGCSM, IBND, IG, IND, IPR, IPRSM
REAL(KIND=JPRB)    :: ZWTSUM

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
logical, intent (in) :: lacc

!#include "susrtmcf.intfb.h"
#include "susrtm.intfb.h"
#include "srtm_kgb16.intfb.h"
#include "srtm_kgb17.intfb.h"
#include "srtm_kgb18.intfb.h"
#include "srtm_kgb19.intfb.h"
#include "srtm_kgb20.intfb.h"
#include "srtm_kgb21.intfb.h"
#include "srtm_kgb22.intfb.h"
#include "srtm_kgb23.intfb.h"
#include "srtm_kgb24.intfb.h"
#include "srtm_kgb25.intfb.h"
#include "srtm_kgb26.intfb.h"
#include "srtm_kgb27.intfb.h"
#include "srtm_kgb28.intfb.h"
#include "srtm_kgb29.intfb.h"
!#include "susrtop.intfb.h"

#include "modify_wv_continuum.intfb.h"

#include "srtm_cmbgb16.intfb.h"
#include "srtm_cmbgb17.intfb.h"
#include "srtm_cmbgb18.intfb.h"
#include "srtm_cmbgb19.intfb.h"
#include "srtm_cmbgb20.intfb.h"
#include "srtm_cmbgb21.intfb.h"
#include "srtm_cmbgb22.intfb.h"
#include "srtm_cmbgb23.intfb.h"
#include "srtm_cmbgb24.intfb.h"
#include "srtm_cmbgb25.intfb.h"
#include "srtm_cmbgb26.intfb.h"
#include "srtm_cmbgb27.intfb.h"
#include "srtm_cmbgb28.intfb.h"
#include "srtm_cmbgb29.intfb.h"

IF (LHOOK) CALL DR_HOOK('SRTM_INIT',0,ZHOOK_HANDLE)

!CALL SUSRTMCF
CALL SUSRTM_GPU(lacc=lacc)

!-- read in the molecular absorption coefficients

CALL SRTM_KGB16(CDIRECTORY)
CALL SRTM_KGB17
CALL SRTM_KGB18
CALL SRTM_KGB19
CALL SRTM_KGB20
CALL SRTM_KGB21
CALL SRTM_KGB22
CALL SRTM_KGB23
CALL SRTM_KGB24
CALL SRTM_KGB25
CALL SRTM_KGB26
CALL SRTM_KGB27
CALL SRTM_KGB28
CALL SRTM_KGB29

IF (PRESENT(NWVCONTINUUM)) THEN
  ! Modify the shortwave water vapour continuum, if requested
  CALL MODIFY_WV_CONTINUUM(NWVCONTINUUM)
END IF

!-- read in the cloud optical properties
!- creates module YOESRTOP with EXTLIQ1, SSALIQ1, ASYLIQ1, 
!  EXTICE3, SSAICE3, ASYICE3, FDLICE3  

!-- RRTM_SW cloud optical properties are not used
!   SRTM_CLDPROP is not called
!   no need to call SUSRTOP

!CALL SUSRTOP ( -1 )


!Mike Iacono 20050804
!-- Perform g-point reduction from 16 per band (224 total points) to
!-- a band dependent number (112 total points) for all absorption
!-- coefficient input data and Planck fraction input data.
!-- Compute relative weighting for new g-point combinations.

IGCSM = 0
DO IBND = 1,JPSW
  IPRSM = 0
  IF (NGC(IBND) < JPG) THEN
    DO IGC = 1,NGC(IBND)
      IGCSM = IGCSM + 1
      ZWTSUM = 0.
      DO IPR = 1, NGN(IGCSM)
        IPRSM = IPRSM + 1
        ZWTSUM = ZWTSUM + WT(IPRSM)
      ENDDO
      WTSM(IGC) = ZWTSUM
    ENDDO

    DO IG = 1,NG(IBND+15)
      IND = (IBND-1)*JPG + IG
      RWGT(IND) = WT(IG)/WTSM(NGM(IND))
    ENDDO
  ELSE
    DO IG = 1,NG(IBND+15)
      IGCSM = IGCSM + 1
      IND = (IBND-1)*JPG + IG
      RWGT(IND) = 1.0
    ENDDO
  ENDIF
ENDDO

CALL SRTM_CMBGB16_GPU(lacc=lacc)
CALL SRTM_CMBGB17_GPU(lacc=lacc)
CALL SRTM_CMBGB18_GPU(lacc=lacc)
CALL SRTM_CMBGB19_GPU(lacc=lacc)
CALL SRTM_CMBGB20_GPU(lacc=lacc)
CALL SRTM_CMBGB21_GPU(lacc=lacc)
CALL SRTM_CMBGB22_GPU(lacc=lacc)
CALL SRTM_CMBGB23_GPU(lacc=lacc)
CALL SRTM_CMBGB24_GPU(lacc=lacc)
CALL SRTM_CMBGB25_GPU(lacc=lacc)
CALL SRTM_CMBGB26_GPU(lacc=lacc)
CALL SRTM_CMBGB27_GPU(lacc=lacc)
CALL SRTM_CMBGB28_GPU(lacc=lacc)
CALL SRTM_CMBGB29_GPU(lacc=lacc)

!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SRTM_INIT',1,ZHOOK_HANDLE)
END SUBROUTINE SRTM_INIT_GPU



SUBROUTINE SRTM_INIT_CPU(CDIRECTORY, NWVCONTINUUM)

!-- read in the basic coefficients to configure RRTM_SW
!- creates module YOESRTWN with BG, NSPA, NSPB, WAVENUM1, WAVENUM2,
!  DELWAVE, PREF, PREFLOG, TREF

USE PARKIND1
USE YOMHOOK

USE PARSRTM
USE YOESRTM  , ONLY : NGN
USE YOESRTWN , ONLY : NG, NGM, WT, NGC, RWGT, WTSM
!USE YOESRTWN , ONLY : NG, NGM, WT, NGC, NGN, RWGT, WTSM
!USE YOMLUN   , ONLY : NULOUT

IMPLICIT NONE

CHARACTER(LEN=*), INTENT(IN) :: CDIRECTORY

! Water vapour continuum model (0=default MT_CKD2.5, 1=CAVIAR)
INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL :: NWVCONTINUUM

! Local variables
INTEGER(KIND=JPIM) :: IGC, IGCSM, IBND, IG, IND, IPR, IPRSM
REAL(KIND=JPRB)    :: ZWTSUM

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!#include "susrtmcf.intfb.h"
#include "susrtm.intfb.h"
#include "srtm_kgb16.intfb.h"
#include "srtm_kgb17.intfb.h"
#include "srtm_kgb18.intfb.h"
#include "srtm_kgb19.intfb.h"
#include "srtm_kgb20.intfb.h"
#include "srtm_kgb21.intfb.h"
#include "srtm_kgb22.intfb.h"
#include "srtm_kgb23.intfb.h"
#include "srtm_kgb24.intfb.h"
#include "srtm_kgb25.intfb.h"
#include "srtm_kgb26.intfb.h"
#include "srtm_kgb27.intfb.h"
#include "srtm_kgb28.intfb.h"
#include "srtm_kgb29.intfb.h"
!#include "susrtop.intfb.h"

#include "modify_wv_continuum.intfb.h"

#include "srtm_cmbgb16.intfb.h"
#include "srtm_cmbgb17.intfb.h"
#include "srtm_cmbgb18.intfb.h"
#include "srtm_cmbgb19.intfb.h"
#include "srtm_cmbgb20.intfb.h"
#include "srtm_cmbgb21.intfb.h"
#include "srtm_cmbgb22.intfb.h"
#include "srtm_cmbgb23.intfb.h"
#include "srtm_cmbgb24.intfb.h"
#include "srtm_cmbgb25.intfb.h"
#include "srtm_cmbgb26.intfb.h"
#include "srtm_cmbgb27.intfb.h"
#include "srtm_cmbgb28.intfb.h"
#include "srtm_cmbgb29.intfb.h"

IF (LHOOK) CALL DR_HOOK('SRTM_INIT',0,ZHOOK_HANDLE)

!CALL SUSRTMCF
CALL SUSRTM_CPU

!-- read in the molecular absorption coefficients

CALL SRTM_KGB16(CDIRECTORY)
CALL SRTM_KGB17
CALL SRTM_KGB18
CALL SRTM_KGB19
CALL SRTM_KGB20
CALL SRTM_KGB21
CALL SRTM_KGB22
CALL SRTM_KGB23
CALL SRTM_KGB24
CALL SRTM_KGB25
CALL SRTM_KGB26
CALL SRTM_KGB27
CALL SRTM_KGB28
CALL SRTM_KGB29

IF (PRESENT(NWVCONTINUUM)) THEN
  ! Modify the shortwave water vapour continuum, if requested
  CALL MODIFY_WV_CONTINUUM(NWVCONTINUUM)
END IF

!-- read in the cloud optical properties
!- creates module YOESRTOP with EXTLIQ1, SSALIQ1, ASYLIQ1, 
!  EXTICE3, SSAICE3, ASYICE3, FDLICE3  

!-- RRTM_SW cloud optical properties are not used
!   SRTM_CLDPROP is not called
!   no need to call SUSRTOP

!CALL SUSRTOP ( -1 )


!Mike Iacono 20050804
!-- Perform g-point reduction from 16 per band (224 total points) to
!-- a band dependent number (112 total points) for all absorption
!-- coefficient input data and Planck fraction input data.
!-- Compute relative weighting for new g-point combinations.

IGCSM = 0
DO IBND = 1,JPSW
  IPRSM = 0
  IF (NGC(IBND) < JPG) THEN
    DO IGC = 1,NGC(IBND)
      IGCSM = IGCSM + 1
      ZWTSUM = 0.
      DO IPR = 1, NGN(IGCSM)
        IPRSM = IPRSM + 1
        ZWTSUM = ZWTSUM + WT(IPRSM)
      ENDDO
      WTSM(IGC) = ZWTSUM
    ENDDO

    DO IG = 1,NG(IBND+15)
      IND = (IBND-1)*JPG + IG
      RWGT(IND) = WT(IG)/WTSM(NGM(IND))
    ENDDO
  ELSE
    DO IG = 1,NG(IBND+15)
      IGCSM = IGCSM + 1
      IND = (IBND-1)*JPG + IG
      RWGT(IND) = 1.0
    ENDDO
  ENDIF
ENDDO

CALL SRTM_CMBGB16_CPU
CALL SRTM_CMBGB17_CPU
CALL SRTM_CMBGB18_CPU
CALL SRTM_CMBGB19_CPU
CALL SRTM_CMBGB20_CPU
CALL SRTM_CMBGB21_CPU
CALL SRTM_CMBGB22_CPU
CALL SRTM_CMBGB23_CPU
CALL SRTM_CMBGB24_CPU
CALL SRTM_CMBGB25_CPU
CALL SRTM_CMBGB26_CPU
CALL SRTM_CMBGB27_CPU
CALL SRTM_CMBGB28_CPU
CALL SRTM_CMBGB29_CPU

!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SRTM_INIT',1,ZHOOK_HANDLE)
END SUBROUTINE SRTM_INIT_CPU

