!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
C
C----- DESCRIPTION DES "PARAMETER" DU LOGICIEL DE FICHIERS INDEXES -----
C
C     JPDBLE= PRECISION UTILISE POUR LES ENTIERS:
C     * SI JPDBLE=8 les INTEGER (KIND=JPDBLE) seront a priori en 64 BITS
C     * SI JPDBLE=4 les INTEGER (KIND=JPDBLE) seront a priori en 32 BITS
C
C     JPDBLR= PRECISION UTILISE POUR LES FLOTTANTS (REELS):
C     * SI JPDBLR=8 les REAL (KIND=JPDBLR) seront a priori en 64 BITS
C     * SI JPDBLR=4 les REAL (KIND=JPDBLR) seront a priori en 32 BITS
C
C     (les conventions peuvent dependre de la plate-forme consideree)
C
C     JP_SIMPLE_ENTIER= sous-type entier permettant de representer
C                       l'intervalle +/-(10**9 - 1)
C
      INTEGER, PARAMETER :: JPDBLE = 8, JPDBLR = 8
      INTEGER, PARAMETER :: JP_SIMPLE_ENTIER = SELECTED_INT_KIND(9)
