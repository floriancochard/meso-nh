#! /bin/sh
FILE=${1:-Bret45.99082200dg.Z}
DIRLFI=${2:-.}
export DIRLFI
#
ARCH=LXNAGf95
B=32
#
rm ${FILE}LLHV
/mesonh/MAKE/tools/diachro/${ARCH}_${B}/extractdia << EOF
$FILE
LLHV 
0
1,10,1,10,2,5
1,1,1,1,1,1
FF
THM
DD
END
EOF
