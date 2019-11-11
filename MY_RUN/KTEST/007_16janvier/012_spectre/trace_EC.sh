EXP=${1:-16JAN}
ECH=${2:-300}
var=${3:-EC}
levd=${4:-10}
levf=${5:-10}
res=${6:-9}




fileracU="spectra_16JAN.1.12B18.001_U"
fileracV="spectra_16JAN.1.12B18.001_V"
fileracW="spectra_16JAN.1.12B18.001_W"
fileracTheta="spectra_16JAN.1.12B18.001_Theta"
fileracRv="spectra_16JAN.1.12B18.001_Rv"
fileracLSU="spectra_16JAN.1.12B18.001_LSU"
fileracLSV="spectra_16JAN.1.12B18.001_LSV"
fileracLSW="spectra_16JAN.1.12B18.001_LSW"
fileracLSTheta="spectra_16JAN.1.12B18.001_LSTheta"
fileracLSRv="spectra_16JAN.1.12B18.001_LSRv"


tronc=`wc -l $fileracU | awk '{print $1}'`

echo "***************************"
echo "On traite les fichiers  :"
ls $fileracU
ls $fileracV
echo " nombres d'onde :"
echo $tronc
echo "entre les niveaux : " $levd " et " $levf
echo "parametre : " $var
echo "***************************"

sleep 2

nivd=`expr $levd + 2`
nivf=`expr $levf + 2`


#echo $nivd
#echo $nivf

TRUNC=`expr $tronc - 1`

#echo $TRUNC
# on moyenne sur les niveua xverticaux choisis
case "$var" in 
 U) 
 tail -$TRUNC $fileracU  | awk -v d=$nivd -v f=$nivf '{s=0; for (i=d; i<=f; i++) s=s+$i; print $1, $2, $3, s/(f - d + 1), s}' > fileU ;;
 V) 
 tail -$TRUNC $fileracV  | awk -v d=$nivd -v f=$nivf '{s=0; for (i=d; i<=f; i++) s=s+$i; print $1, $2, $3, s/(f - d + 1), s}' > fileV ;;
 W) 
 tail -$TRUNC $fileracW  | awk -v d=$nivd -v f=$nivf '{s=0; for (i=d; i<=f; i++) s=s+$i; print $1, $2, $3, s/(f - d + 1), s}' > fileW ;;
 EC )
 tail -$TRUNC $fileracU  | awk -v d=$nivd -v f=$nivf '{s=0; for (i=d; i<=f; i++) s=s+$i; print $1, $2, $3, s/(f - d + 1), s}' > fileU 
 tail -$TRUNC $fileracV  | awk -v d=$nivd -v f=$nivf '{s=0; for (i=d; i<=f; i++) s=s+$i; print $1, $2, $3, s/(f - d + 1), s}' > fileV 
 paste fileU fileV > interm
 awk '{print $1, $2, $3,($4+$9)/2}' interm > fileEC
 rm fileU fileV interm ;;
 Theta)
 tail -$TRUNC $fileracTheta  | awk -v d=$nivd -v f=$nivf '{s=0; for (i=d; i<=f; i++) s=s+$i; print $1, $2, $3, s/(f - d + 1), s}' > fileTheta ;; 
 Rv)
 tail -$TRUNC $fileracRv  | awk -v d=$nivd -v f=$nivf '{s=0; for (i=d; i<=f; i++) s=s+$i; print $1, $2, $3, s/(f - d + 1), s}' > fileRv ;; 
 LSU)
 tail -$TRUNC $fileracLSU  | awk -v d=$nivd -v f=$nivf '{s=0; for (i=d; i<=f; i++) s=s+$i; print $1, $2, $3, s/(f - d + 1), s}' > fileLSU ;; 
 LSV)
 tail -$TRUNC $fileracLSV  | awk -v d=$nivd -v f=$nivf '{s=0; for (i=d; i<=f; i++) s=s+$i; print $1, $2, $3, s/(f - d + 1), s}' > fileLSV ;;
 LSEC )
 tail -$TRUNC $fileracLSU  | awk -v d=$nivd -v f=$nivf '{s=0; for (i=d; i<=f; i++) s=s+$i; print $1, $2, $3, s/(f - d + 1), s}' > fileLSU 
 tail -$TRUNC $fileracLSV  | awk -v d=$nivd -v f=$nivf '{s=0; for (i=d; i<=f; i++) s=s+$i; print $1, $2, $3, s/(f - d + 1), s}' > fileLSV 
 paste fileLSU fileLSV > interm
 awk '{print $1, $2, $3,($4+$9)/2}' interm > fileLSEC
 rm fileLSU fileLSV interm ;;
 LSW)
 tail -$TRUNC $fileracLSW  | awk -v d=$nivd -v f=$nivf '{s=0; for (i=d; i<=f; i++) s=s+$i; print $1, $2, $3, s/(f - d + 1), s}' > fileLSW ;; 
 LSTheta)
 tail -$TRUNC $fileracLSTheta  | awk -v d=$nivd -v f=$nivf '{s=0; for (i=d; i<=f; i++) s=s+$i; print $1, $2, $3, s/(f - d + 1), s}' > fileLSTheta ;; 
 LSrv)
 tail -$TRUNC $fileracLSRv  | awk -v d=$nivd -v f=$nivf '{s=0; for (i=d; i<=f; i++) s=s+$i; print $1, $2, $3, s/(f - d + 1), s}' > fileLSRv ;;
esac

pos1=`echo " " | awk -v a=5 -v b=$tronc -v c=$res '{print 2*c*b/a}'`
pos2=`echo " " | awk -v a=10 -v b=$tronc -v c=$res '{print 2*c*b/a}'`
pos3=`echo " " | awk -v a=15 -v b=$tronc -v c=$res '{print 2*c*b/a}'`
pos4=`echo " " | awk -v a=20 -v b=$tronc -v c=$res '{print 2*c*b/a}'`
pos5=`echo " " | awk -v a=50 -v b=$tronc -v c=$res '{print 2*c*b/a}'`
pos6=`echo " " | awk -v a=100 -v b=$tronc -v c=$res '{print 2*c*b/a}'`
pos7=`echo " " | awk -v a=200 -v b=$tronc -v c=$res '{print 2*c*b/a}'`
pos8=`echo " " | awk -v a=500 -v b=$tronc -v c=$res '{print 2*c*b/a}'`
pos9=`echo " " | awk -v a=1000 -v b=$tronc -v c=$res '{print 2*c*b/a}'`
pos10=`echo " " | awk -v a=25 -v b=$tronc -v c=$res '{print 2*c*b/a}'`
pos11=`echo " " | awk -v a=30 -v b=$tronc -v c=$res '{print 2*c*b/a}'`
pos12=`echo " " | awk -v a=40 -v b=$tronc -v c=$res '{print 2*c*b/a}'`
pos13=`echo " " | awk -v a=60 -v b=$tronc -v c=$res '{print 2*c*b/a}'`
pos14=`echo " " | awk -v a=70 -v b=$tronc -v c=$res '{print 2*c*b/a}'`
pos15=`echo " " | awk -v a=80 -v b=$tronc -v c=$res '{print 2*c*b/a}'`
pos16=`echo " " | awk -v a=90 -v b=$tronc -v c=$res '{print 2*c*b/a}'`
pos17=`echo " " | awk -v a=300 -v b=$tronc -v c=$res '{print 2*c*b/a}'`
pos18=`echo " " | awk -v a=400 -v b=$tronc -v c=$res '{print 2*c*b/a}'`
pos19=`echo " " | awk -v a=600 -v b=$tronc -v c=$res '{print 2*c*b/a}'`
pos20=`echo " " | awk -v a=700 -v b=$tronc -v c=$res '{print 2*c*b/a}'`
pos21=`echo " " | awk -v a=800 -v b=$tronc -v c=$res '{print 2*c*b/a}'`
pos22=`echo " " | awk -v a=900 -v b=$tronc -v c=$res '{print 2*c*b/a}'`
pos23=`echo " " | awk -v a=6 -v b=$tronc -v c=$res '{print 2*c*b/a}'`
pos24=`echo " " | awk -v a=7 -v b=$tronc -v c=$res '{print 2*c*b/a}'`
pos25=`echo " " | awk -v a=8 -v b=$tronc -v c=$res '{print 2*c*b/a}'`
pos26=`echo " " | awk -v a=9 -v b=$tronc -v c=$res '{print 2*c*b/a}'`
pos27=`echo " " | awk -v a=4 -v b=$tronc -v c=$res '{print 2*c*b/a}'`
pos28=`echo " " | awk -v a=3 -v b=$tronc -v c=$res '{print 2*c*b/a}'`
pos29=`echo " " | awk -v a=2 -v b=$tronc -v c=$res '{print 2*c*b/a}'`
pos30=`echo " " | awk -v a=1 -v b=$tronc -v c=$res '{print 2*c*b/a}'`
pos31=`echo " " | awk -v a=0.5 -v b=$tronc -v c=$res '{print 2*c*b/a}'`
pos32=`echo " " | awk -v a=0.2 -v b=$tronc -v c=$res '{print 2*c*b/a}'`
pos33=`echo " " | awk -v a=0.6 -v b=$tronc -v c=$res '{print 2*c*b/a}'`
pos34=`echo " " | awk -v a=0.7 -v b=$tronc -v c=$res '{print 2*c*b/a}'`
pos35=`echo " " | awk -v a=0.8 -v b=$tronc -v c=$res '{print 2*c*b/a}'`
pos36=`echo " " | awk -v a=0.9 -v b=$tronc -v c=$res '{print 2*c*b/a}'`






echo $pos1
echo $pos2
echo $pos3
echo $pos4
echo $pos5
echo $pos6
echo $pos7
echo $pos8
echo $pos9




gnuplot<<EOF
set encoding iso_8859_1
set terminal postscript eps enhanced color  dashed dashlength 3.0  "Times-Roman" 18
set output "spectre.eps"
set logscale y
set logscale x
set logscale x2
set grid
set xtics nomirror
set xlabel "Nombre d'onde"
set ylabel "Energie cinetique (m^2.s^-^2)"
set x2label "Longueur d'onde (km)"

set xrange[1:2*$tronc]
set yrange[1e-6:100]


set x2tics ("5" $pos1, "10" $pos2,"15" $pos3,"20" $pos4,"50"  $pos5,"100" $pos6,"200" $pos7,"500" $pos8,"1000" $pos9,"25" $pos10,"30" $pos11,"40" $pos12,"" $pos13,"" $pos14,"" $pos15,"" $pos16,"" $pos17,"" $pos18,"" $pos19,"" $pos20,"" $pos21,"" $pos22,"" $pos23,"" $pos24,"" $pos25,"" $pos26,"" $pos27,"" $pos28,"" $pos29,"1" $pos30,"0.5" $pos31,"" $pos32,"" $pos33,"" $pos34,"" $pos35,"" $pos36)


set style line 1 lt 3 lc rgb "dark-gray" lw 1
set style line 2 lt 5 lc rgb "dark-gray" lw 1
set style line 3 lt 1 lc rgb "red" lw 5





plot "file${var}" using 1:4  ls 3 axis x1y1 title "${EXP} : $var" with lines smooth bezier, 10*x**(-5./3.) ls 1 title "-5/3" , 50*x**(-3) ls 2 title "-3"


EOF

mv spectre.eps spectre_${EXP}_${var}_${ECH}_n${levd}_${levf}.eps
rm -f file${var}
