C/**************************************************************************
C Include file: binsolu.h  
C
CPurpose: concentration of solute (micromol solute / microgram water)
C    in a binary solution at RH = Aw = 0, 0.1, 0.2, 0.3 etc
C
CInclude dependencies: used in zsr.f
C
CRevision history: Developed by Betty Pun, December 1999, under CARB funding
C                  Reconfigured for Fortran and only 5 Type A species by
C              Rob Griffin January, 2000
C                  Reconfigured for combined philic-phobic module, April, 2001
C***********************************************************************/
C        graduation of RH scale in molalbin definition, 
C      starting with RH = 0 (first index of first dim) 
C                     ending with 1 (last index of first dim)
C
      RHGRAD = 0.1

C       binary solution molality (umol/ug water) corresponding to RH 
C       at 10% graduations --> 11 rows
C
      MOLALBIN(0,1) = 555.6
      MOLALBIN(1,1) = 0.65000
      MOLALBIN(2,1) = 0.30762
      MOLALBIN(3,1) = 0.19305
      MOLALBIN(4,1) = 0.13537
      MOLALBIN(5,1) = 0.10038
      MOLALBIN(6,1) = 0.07666
      MOLALBIN(7,1) = 0.05927
      MOLALBIN(8,1) = 0.04567
      MOLALBIN(9,1) = 0.03431
      MOLALBIN(10,1) = 0.0
      MOLALBIN(0,2) = 555.6
      MOLALBIN(1,2) = 1.98908
      MOLALBIN(2,2) = 0.95252
      MOLALBIN(3,2) = 0.60644
      MOLALBIN(4,2) = 0.43292
      MOLALBIN(5,2) = 0.32835
      MOLALBIN(6,2) = 0.25820
      MOLALBIN(7,2) = 0.20764
      MOLALBIN(8,2) = 0.16923
      MOLALBIN(9,2) = 0.13830
      MOLALBIN(10,2) = 0.0
      MOLALBIN(0,3) = 555.6
      MOLALBIN(1,3) = 63.80219
      MOLALBIN(2,3) = 31.65389
      MOLALBIN(3,3) = 20.93034
      MOLALBIN(4,3) = 15.57215
      MOLALBIN(5,3) = 12.35612
      MOLALBIN(6,3) = 10.21267
      MOLALBIN(7,3) = 8.68083
      MOLALBIN(8,3) = 7.53173
      MOLALBIN(9,3) = 6.63831
      MOLALBIN(10,3) = 0.0
      MOLALBIN(0,4) = 555.6
      MOLALBIN(1,4) = 2.81667
      MOLALBIN(2,4) = 1.37478
      MOLALBIN(3,4) = 0.89336
      MOLALBIN(4,4) = 0.65205
      MOLALBIN(5,4) = 0.50677
      MOLALBIN(6,4) = 0.40951
      MOLALBIN(7,4) = 0.33965
      MOLALBIN(8,4) = 0.28692
      MOLALBIN(9,4) = 0.24558
      MOLALBIN(10,4) = 0.0
      MOLALBIN(0,5) = 555.6
      MOLALBIN(1,5) = 1.37883
      MOLALBIN(2,5) = 0.63792
      MOLALBIN(3,5) = 0.39248
      MOLALBIN(4,5) = 0.27079
      MOLALBIN(5,5) = 0.19847
      MOLALBIN(6,5) = 0.15067
      MOLALBIN(7,5) = 0.11672
      MOLALBIN(8,5) = 0.09124
      MOLALBIN(9,5) = 0.07117
      MOLALBIN(10,5) = 0.0
      MOLALBIN(0,6) = 555.6
      MOLALBIN(1,6) = 0.49396
      MOLALBIN(2,6) = 0.22574
      MOLALBIN(3,6) = 0.13569
      MOLALBIN(4,6) = 0.09018
      MOLALBIN(5,6) = 0.06232
      MOLALBIN(6,6) = 0.04319
      MOLALBIN(7,6) = 0.0289
      MOLALBIN(8,6) = 0.0174
      MOLALBIN(9,6) = 0.00755
      MOLALBIN(10,6)= 0.0
      MOLALBIN(0,7) = 555.6
      MOLALBIN(1,7) = 0.45742
      MOLALBIN(2,7) = 0.21651
      MOLALBIN(3,7) = 0.13549
      MOLALBIN(4,7) = 0.09439
      MOLALBIN(5,7) = 0.06918
      MOLALBIN(6,7) = 0.05184
      MOLALBIN(7,7) = 0.03891
      MOLALBIN(8,7) = 0.02848
      MOLALBIN(9,7) = 0.01921
      MOLALBIN(10,7)= 0.0
      MOLALBIN(0,8) = 555.6
      MOLALBIN(1,8) = 0.86272
      MOLALBIN(2,8) = 0.40057
      MOLALBIN(3,8) = 0.24605
      MOLALBIN(4,8) = 0.16855
      MOLALBIN(5,8) = 0.1216
      MOLALBIN(6,8) = 0.08988
      MOLALBIN(7,8) = 0.06671
      MOLALBIN(8,8) = 0.0486
      MOLALBIN(9,8) = 0.03329
      MOLALBIN(10,8)= 0.0
      MOLALBIN(0,9) = 555.6
      MOLALBIN(1,9) = 0.55832
      MOLALBIN(2,9) = 0.263
      MOLALBIN(3,9) = 0.16412
      MOLALBIN(4,9) = 0.11418
      MOLALBIN(5,9) = 0.08375
      MOLALBIN(6,9) = 0.06308
      MOLALBIN(7,9) = 0.0478
      MOLALBIN(8,9) = 0.03574
      MOLALBIN(9,9) = 0.02543
      MOLALBIN(10,9)= 0.0
      MOLALBIN(0,10) = 555.6
      MOLALBIN(1,10) = 0.92947
      MOLALBIN(2,10) = 0.42461
      MOLALBIN(3,10) = 0.25726
      MOLALBIN(4,10) = 0.17392
      MOLALBIN(5,10) = 0.12435
      MOLALBIN(6,10) = 0.09149
      MOLALBIN(7,10) = 0.06809
      MOLALBIN(8,10) = 0.05035
      MOLALBIN(9,10) = 0.03601
      MOLALBIN(10,10)= 0.0
