!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 microph 2006/05/18 13:07:25
!-----------------------------------------------------------------
!     ######spl
      MODULE MODI_CONDENS
!     ###################
!
INTERFACE
!
     SUBROUTINE CONDENS(HTURBDIM, PQ1, PN, PRC, PSRC)

CHARACTER*4,             INTENT(IN)  ::   HTURBDIM ! dimensionality of the
                                                   ! turbulence scheme
REAL, DIMENSION(:,:,:),  INTENT(IN)  ::   PQ1      ! Saturation
REAL, DIMENSION(:,:,:),  INTENT(OUT) ::   PN       ! Cloud fraction
REAL, DIMENSION(:,:,:),  INTENT(OUT) ::   PRC      ! Cloud water mixing ratio
                                                   ! rc/2Sigma_s
REAL, DIMENSION(:,:,:),  INTENT(OUT) ::   PSRC     ! Second-order flux 
                                                   ! s'rc'/2 Sigma_s2
!
END SUBROUTINE CONDENS
!
END INTERFACE
!
END MODULE MODI_CONDENS
!     ######spl
     SUBROUTINE CONDENS(HTURBDIM, PQ1, PN, PRC, PSRC)
!    #################################################
!
!!****    *CONDENS* - Computes the normalized cloud parameters in a subgrid 
!!                        condensation scheme
!!
!!    PURPOSE
!!    -------
!       The purpose of this routine is to compute the following normalized 
!     cloud parameters:
!                               cloud fraction,   
!                               cloud water mixing ratio,
!                               second-order moment s'rc'.
!        They are obtained by using a statistical distribution for the
!      condensation determined by its saturation  threshold and its asymmetry. 
!!**   METHOD
!!     ------
!!       The parameters are computed by empirical determination fitted
!!       from 3D numerical simulation results.  The input parameters
!!       are the degree of saturation Q1 and the asymmetry AS. In fact,
!!       The asymmetry is a prescribed function of Q1:
!!                  As = MIN(2., MAX( 0., 1-Q1 ) )  for 1D turbulence
!!                  As = 0                          for 3D turbulence
!!
!!       The outputs N, rc, src are computed according to a gamma distribution
!!       These computations are very time consumming and their results have 
!!       therefore been tabulated. The different representative curves of 
!!       N, rc and src are plotted in the scientific documentation for both 
!!       cases.  
!!
!!     EXTERNAL
!!     --------
!!       NONE
!!
!!     IMPLICIT ARGUMENTS
!!     ------------------
!!       NONE
!!   
!!     REFERENCE
!!     ---------
!!       Book 1 of the  Meso-NH documentation
!!
!!
!!     AUTHOR
!!     ------
!!       Philippe Bougeault        * Meteo-France *
!!
!!     MODIFICATIONS
!!     -------------
!!      Original       1981 ?  Fortran V
!!      Modifications: March 2, 1995 (J.M. Carriere) Fortran 90
!!                                   and doctorization
!!      Modifications: September 12, 1996 (J. Stein) tabulated version
!!
!! ----------------------------------------------------------------------
!
!*       0. DECLARATIONS
!           ------------
!
!
IMPLICIT NONE
!
!*       0.1 declarations of arguments and result
!
CHARACTER*4,             INTENT(IN)  ::   HTURBDIM ! dimensionality of the
                                                   ! turbulence scheme
REAL, DIMENSION(:,:,:),  INTENT(IN)  ::   PQ1      ! Saturation
REAL, DIMENSION(:,:,:),  INTENT(OUT) ::   PN       ! Cloud fraction
REAL, DIMENSION(:,:,:),  INTENT(OUT) ::   PRC      ! Cloud water mixing ratio
                                                   ! rc/2Sigma_s
REAL, DIMENSION(:,:,:),  INTENT(OUT) ::   PSRC     ! Second-order flux 
                                                   ! s'rc'/2 Sigma_s2
!
!*       0.2 declarations of local variables
!
INTEGER           :: JR1,JR2,JR3, IRM1,IRM2,IRM3
!
REAL, DIMENSION(-22:11) :: ZN_1D = (/                                   &
       0.           ,  0.           ,  1.7225742E-05,  0.275373E-04 ,   &
       4.5657158E-05,  0.748634E-04 ,  1.2344122E-04,  0.203788E-03 ,   &
       3.3539534E-04,  0.553310E-03 ,  9.1189146E-04,  0.150353E-02 ,   &
       2.4784803E-03,  0.408673E-02 ,  6.7381263E-03,  0.111092E-01 ,   &
       1.8315554E-02,  0.301974E-01 ,  4.9787164E-02,  0.831191E-01 ,   &
       0.1512039    ,  0.286653E+00 ,  0.5000000    ,  0.691489E+00 ,   &
       0.8413813    ,  0.933222E+00 ,  0.9772662    ,  0.993797E+00 ,   &
       0.9986521    ,  0.999768E+00 ,  0.9999684    ,  0.9999997    ,   &
       1.0000000    ,  1.000000      /)
!
REAL, DIMENSION(-22:11) :: ZRC_1D = (/                                   &
       0.           ,  0.           ,  1.1461278E-05,  0.275279E-04 ,    &
       4.3084903E-05,  0.747532E-04 ,  1.2315845E-04,  0.201069E-03 ,    &
       3.3593364E-04,  0.551618E-03 ,  9.1182487E-04,  0.150296E-02 ,    &
       2.4801120E-03,  0.408695E-02 ,  6.7372285E-03,  0.111084E-01 ,    &
       1.8315896E-02,  0.301974E-01 ,  4.9786866E-02,  0.721706E-01 ,    &
       0.1165014    ,  0.210263E+00 ,  0.3990000    ,  0.697847E+00 ,    &
       1.0833505    ,  0.152933E+01 ,  2.0084987    ,  0.250201E+01 ,    &
       3.0003829    ,  0.350006E+01 ,  4.0000072    ,  0.450000E+01 ,    &
       5.0000000    ,  5.500000     /)
!
REAL, DIMENSION(-22:11) :: ZSRC_1D =(/                                   &
       0.           ,  0.           ,  2.0094444E-04,   0.316670E-03,    &
       4.9965648E-04,  0.785956E-03 ,  1.2341294E-03,   0.193327E-02,    &
       3.0190963E-03,  0.470144E-02 ,  7.2950651E-03,   0.112759E-01,    &
       1.7350994E-02,  0.265640E-01 ,  4.0427860E-02,   0.610997E-01,    &
       9.1578111E-02,  0.135888E+00 ,  0.1991484    ,   0.230756E+00,    &
       0.2850565    ,  0.375050E+00 ,  0.5000000    ,   0.691489E+00,    &
       0.8413813    ,  0.933222E+00 ,  0.9772662    ,   0.993797E+00,    &
       0.9986521    ,  0.999768E+00 ,  0.9999684    ,   0.999997E+00,    &
       1.0000000    ,  1.000000     /)
!
REAL, DIMENSION(-22:11) :: ZN_3D = (/                                   &
       0.           ,  0.           ,  0.           ,  0.           ,   &
       0.           ,  0.           ,  0.           ,  0.           ,   &
       0.           ,  0.           ,  0.           ,  0.298023E-07 ,   &
       0.298023E-06 ,  0.339746E-05 ,  0.315905E-04 ,  0.232160E-03 ,   &
       0.134790E-02 ,  0.620306E-02 ,  0.227338E-01 ,  0.667779E-01 ,   &
       0.158619E+00 ,  0.308511E+00 ,  0.5000000    ,  0.691489E+00 ,   &
       0.8413813    ,  0.933222E+00 ,  0.9772662    ,  0.993797E+00 ,   &
       0.9986521    ,  0.999768E+00 ,  0.9999684    ,  0.9999997    ,   &
       1.0000000    ,  1.000000      /)
!
REAL, DIMENSION(-22:11) :: ZRC_3D = (/                                   &
       0.           ,  0.           ,  0.           ,  0.           ,    &
       0.           ,  0.           ,  0.           ,  0.           ,    &
       0.           ,  0.           ,  0.           ,  0.           ,    &
       0.           ,  0.648644E-06 ,  0.716466E-05 ,  0.586357E-04 ,    &
       0.382772E-03 ,  0.200666E-02 ,  0.849849E-02 ,  0.293255E-01 ,    &
       0.833505E-01 ,  0.197847E+00 ,  0.3990000    ,  0.697847E+00 ,    &
       1.0833505    ,  0.152933E+01 ,  2.0084987    ,  0.250201E+01 ,    &
       3.0003829    ,  0.350006E+01 ,  4.0000072    ,  0.450000E+01 ,    &
       5.0000000    ,  5.500000     /)
!
REAL, DIMENSION(-22:11) :: ZSRC_3D =(/                                   &
       0.           ,  0.           ,  0.           ,   0.          ,    &
       0.           ,  0.           ,  0.           ,   0.          ,    &
       0.           ,  0.           ,  0.           ,   0.298023E-07,    &
       0.298023E-06 ,  0.339746E-05 ,  0.315905E-04 ,   0.232160E-03,    &
       0.134790E-02 ,  0.620306E-02 ,  0.227338E-01 ,   0.667779E-01,    &
       0.158619E+00 ,  0.308511E+00 ,  0.5000000    ,   0.691489E+00,    &
       0.8413813    ,  0.933222E+00 ,  0.9772662    ,   0.993797E+00,    &
       0.9986521    ,  0.999768E+00 ,  0.9999684    ,   0.999997E+00,    &
       1.0000000    ,  1.000000     /)
!
REAL,     DIMENSION(SIZE(PQ1,1),SIZE(PQ1,2),SIZE(PQ1,3)) :: ZINC
INTEGER , DIMENSION(SIZE(PQ1,1),SIZE(PQ1,2),SIZE(PQ1,3)) :: INQ1
!
! ---------------------------------------------------------------------------
!
IRM1 = SIZE(PQ1(:,:,:),1)
IRM2 = SIZE(PQ1(:,:,:),2)
IRM3 = SIZE(PQ1(:,:,:),3)
!
IF (HTURBDIM=='1DIM') THEN
  DO JR3 = 1,IRM3
    DO JR2 = 1,IRM2
      DO JR1 = 1,IRM1
        !
        INQ1(JR1,JR2,JR3)= MIN( MAX(-22,FLOOR(2*PQ1(JR1,JR2,JR3)) ), 10)
        ZINC(JR1,JR2,JR3)= 2.*PQ1(JR1,JR2,JR3) - INQ1(JR1,JR2,JR3)
        !
        PN(JR1,JR2,JR3) =  MIN(  1.,                                    &
           ( 1. - ZINC(JR1,JR2,JR3) ) * ZN_1D( INQ1(JR1,JR2,JR3)  )     &
         + ZINC(JR1,JR2,JR3)          * ZN_1D( INQ1(JR1,JR2,JR3) + 1)   &
                              )
        !
        PRC(JR1,JR2,JR3) =                                                &
           ( 1. - ZINC(JR1,JR2,JR3) ) * ZRC_1D( INQ1(JR1,JR2,JR3)  )      &
         + ZINC(JR1,JR2,JR3)          * ZRC_1D( INQ1(JR1,JR2,JR3) + 1) 
        !     
        PSRC(JR1,JR2,JR3) =  MIN(  1.,                                      &
           (1. - ZINC(JR1,JR2,JR3) ) * ZSRC_1D( INQ1(JR1,JR2,JR3)  )        &
         + ZINC(JR1,JR2,JR3)         * ZSRC_1D( INQ1(JR1,JR2,JR3) + 1)      &
                                )
      END DO
    END DO
  END DO
ELSE
  DO JR3 = 1,IRM3
    DO JR2 = 1,IRM2
      DO JR1 = 1,IRM1
        !
        INQ1(JR1,JR2,JR3)= MIN( MAX(-22,FLOOR(2*PQ1(JR1,JR2,JR3)) ), 10)
        ZINC(JR1,JR2,JR3)= 2.*PQ1(JR1,JR2,JR3) - INQ1(JR1,JR2,JR3)
        !
        PN(JR1,JR2,JR3) =  MIN(  1.,                                    &
           ( 1. - ZINC(JR1,JR2,JR3) ) * ZN_3D( INQ1(JR1,JR2,JR3)  )     &
         + ZINC(JR1,JR2,JR3)          * ZN_3D( INQ1(JR1,JR2,JR3) + 1)   &
                              )
        !
        PRC(JR1,JR2,JR3) =                                                &
           ( 1. - ZINC(JR1,JR2,JR3) ) * ZRC_3D( INQ1(JR1,JR2,JR3)  )      &
         + ZINC(JR1,JR2,JR3)          * ZRC_3D( INQ1(JR1,JR2,JR3) + 1) 
        !     
        PSRC(JR1,JR2,JR3) =  MIN(  1.,                                      &
           (1. - ZINC(JR1,JR2,JR3) ) * ZSRC_3D( INQ1(JR1,JR2,JR3)  )        &
         + ZINC(JR1,JR2,JR3)         * ZSRC_3D( INQ1(JR1,JR2,JR3) + 1)      &
                                )
      END DO
    END DO
  END DO
END IF
!
!
!----------------------------------------------------------------------------
!
END SUBROUTINE CONDENS
