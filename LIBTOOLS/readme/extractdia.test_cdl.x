#! /bin/sh
FILE=${1:-Bret45.99082200dg.Z}
#DIRLFI=${2:-.}
export DIRLFI
#
ARCH=LXNAGf95
B=32
#
rm ${FILE}*zc*
/mesonh/MAKE/tools/diachro/${ARCH}_${B}/extractdia << EOF
$FILE
ZCDL 
5
1,10,1,10
1,1,1,1,1,1
3
1500 3000 5000
LALO
LAT
ALT
LON
END
EOF
