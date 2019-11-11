#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

int main(int argc, char** argv) {
// PRG=ascii2llv ; gcc -lm -Wall -o $PRG ${PRG}.c && chmod u+x $PRG
// convertit un fichier ascii r2 en fichier llv
// arg 1: fichier asc ; 2 fichier llv
  int i,j;
  int lat,lon;
  float **val;
  char s1[20],s2[20],s3[20];
  FILE *fin,*fout;
  
  if((fin=fopen(argv[1],"r"))==NULL) printf("Failed to open %s!!!\n",argv[1]);
  if((fout=fopen(argv[2],"w"))==NULL) printf("Failed to open %s!!!\n",argv[2]);
  
// lecture entête
  fscanf(fin,"%s %s %s %s %s %s %s %s\n",s1,s1,s1,s1,s1,s1,s2,s3);

  lat=atoi(s2);
  lon=atoi(s3);
  printf("%d %d\n",lat,lon);
  
// allocation valeurs
  val=(float **)malloc(lon*sizeof(float*));
  
// allocation et lecture valeurs
  for(i=0;i<lon;i++) {
    val[i]=(float *)malloc(lat*sizeof(float));
    for(j=0;j<lat;j++) {
      fscanf(fin,"%s",s1);
      val[i][j]=atof(s1);      
    }
  }
  
// lecture/écriture latlon
  for(i=0;i<lon;i++) {
    for(j=0;j<lat;j++) {
      fscanf(fin,"%s %s",s1,s2);
      fprintf(fout,"%s %s %f\n",s1,s2,val[i][j]);
    }
    free(val[i]);    
  }
  free(val);
  
  fclose(fin);
  fclose(fout);

  return 1;
} /* main */
