#! /bin/sh
FILE=${1:-16J36.1.00A12.001dg.Z}
DIRLFI=${2:-DATA}
export DIRLFI
#
ARCH=LXNAGf95
B=32
#
rm $(basename $FILE .Z)2.lfi
/mesonh/MAKE/tools/diachro/${ARCH}_${B}/extractdia << EOF
$FILE
DIAC 
1
30,50,20,40,0,0
1,1,1,1
FF
THM
DD
ALT
END
EOF
