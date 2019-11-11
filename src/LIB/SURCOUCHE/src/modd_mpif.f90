!MNH_LIC Copyright 1994-2019 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!  J.Escobar 5/06/2018 : add cpp key MNH_USE_MPI_STATUSES_IGNORE for use of true MPI_STATUSES_IGNORE
!                        & bypass bug with ifort+openmpi
!  P.Wautelet 22/01/2019: replace double precision declarations by real(kind(0.0d0)) (to allow compilation by NAG compiler)
!-----------------------------------------------------------------
MODULE MODD_MPIF
#ifdef USE_MPI
  USE MPI
  IMPLICIT NONE
#else
  IMPLICIT NONE
  INCLUDE 'mpif.h'
#ifdef MNH_USE_MPI_STATUSES_IGNORE
  ! bypass ifort bug with use only MNH_STATUSES_IGNORE => MPI_STATUSES_IGNORE
  real(kind(0.0d0)) XXXXXX
  equivalence ( MPI_STATUSES_IGNORE , XXXXXX )
#endif
#endif
END MODULE MODD_MPIF
