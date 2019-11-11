!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
C****
C            !--------------------------------------------------!
C            !        Sous-programme du logiciel LFI            !
C            ! (Logiciel de Fichiers Indexes par nom d'article) !
C            !--------------------------------------------------!
C
C       - Version originale de LFI: Octobre 1989, auteur:
C                                   Jean CLOCHARD, METEO FRANCE.
C
C       - Aout 1991: Ajout de la notion de "facteur multiplicatif"
C         (on sait traiter un fichier dont la longueur d'article
C          "physique" est multiple de la longueur elementaire JPLARD),
C         et (sur option) toute la messagerie peut etre en anglais.
C
