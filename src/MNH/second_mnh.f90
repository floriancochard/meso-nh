!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 tools 2006/05/18 13:07:25
!-----------------------------------------------------------------
SUBROUTINE SECOND_MNH(XT)

IMPLICIT NONE
REAL           :: XT
REAL           :: ZT
call cpu_time(ZT)
XT=ZT
END SUBROUTINE SECOND_MNH

