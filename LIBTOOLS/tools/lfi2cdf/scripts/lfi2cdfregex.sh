#!/bin/sh
#
#
usage(){
    cat >&2 <<EOF
Usage : 

  ${0##*/} '~/pattern/'  infile.lfi : select articles that match regex 'pattern'.
  ${0##*/} '!~/pattern/' infile.lfi : select articles that doesn't match regex 'pattern'.
    
Example :
  - Select all COVER articles :
      ${0##*/} '~/^COVER/' infile.lfi 

EOF
    exit 1
}

[ -z "$2" ] && usage  

REGEXP=$1
INFILE=$2


VARLIST=$(lfi2cdf -l $INFILE | awk -F\" '$2 && gsub("[[:space:]]+","",$2)+1 && $2 '$REGEXP' {printf("%s,",$2)}')
[ -n "$VARLIST" ] && VARLIST="-v$VARLIST" 
CMD="lfi2cdf $VARLIST $INFILE"
echo $CMD
#$CMD


