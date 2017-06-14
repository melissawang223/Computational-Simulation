#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define nbin 50
#define nkvec 10 
#define PI 3.14159265359

int main(void)
{
  int nprocs, a, b, i, j, k, nstep, nframe, is, ii, il; 
  int *nary, *tary, nglob, currCount;
  double *rary, *vary;
  double rho, fac, gl[3], RC, DRC, rr[3], qq[3], dist, gr, cn, dk[3]; 
  double rrnrm, qqnrm, theta, cosine, sine; 
  double histgram[nbin]={0.0};
  double sq[nkvec][nkvec][nkvec]={0.0};

  FILE *fp, *fout; 
  char fname[20];

  printf("Input: nprocs, first, interval, last\n"); 
  scanf("%d %d %d %d", &nprocs, &is, &ii, &il);
  printf("nprocs = %d, first = %d, interval = %d, last = %d\n",nprocs, is,ii,il); 


  for(nframe=0, nstep=is; nstep<=il; nstep+=ii, nframe++) {
     sprintf(fname, "data/pmd%06d",nstep);

     fp = fopen(fname,"ro");
     if(fp) {
       fscanf(fp,"%d %d %lf %lf %lf\n",&nglob,&currCount,&gl[0],&gl[1],&gl[2]);
     } else {
       fclose(fp); 
       printf("ERROR -----------------------------------------------\n");
       printf("ERROR %s doesn't exsit. Please check 'data' dir again\n"); 
       printf("ERROR -----------------------------------------------\n");
       return 0;
     }
 
     // initialization, compute cutoff length, allocate memory
     if(nstep==is) {

       RC  = gl[0]; 
       RC  = gl[1] > RC ? gl[1] : RC; 
       RC  = gl[2] > RC ? gl[2] : RC; 
       RC  = 0.5*RC; 
       DRC = RC/nbin; 
       rho = nglob/(gl[0]*gl[1]*gl[2]);

       printf("-----------------------------------------------------\n");
       printf("nglob = %d currCount = %d \n",nglob,currCount);
       printf("gl[0:2] = (%lf, %lf, %lf)\n",gl[0],gl[1],gl[2]);
       printf("RC = %f, DRC = %f \n", RC, DRC);
       printf("nbin = %d, rho = %f\n", nbin, rho);
       printf("-----------------------------------------------------\n");

       nary=malloc(sizeof(long int)*nprocs);
       tary=malloc(sizeof(int)*nglob);
       rary=malloc(sizeof(double)*3*nglob);
       vary=malloc(sizeof(double)*3*nglob);

     // get reciplocal lattice vector 
       dk[0]=2.0*PI/gl[0]/nkvec;
       dk[1]=2.0*PI/gl[1]/nkvec;
       dk[2]=2.0*PI/gl[2]/nkvec;
       printf("dk[0:2] %lf %lf %lf\n", dk[0],dk[1],dk[2]);
     }

     //printf("--- %s ---\n",fname);

     printf("reading data from %s \n", fname);

     for(a=0; a<nprocs; a++) fscanf(fp,"%d",&nary[a]);
     for(a=0; a<nglob; a++) {
        fscanf(fp,"%d %lf %lf %lf %lf %lf %lf\n",
        &tary[a], &rary[3*a],&rary[3*a+1],&rary[3*a+2],&vary[3*a],&vary[3*a+1],&vary[3*a+2]);
     }

     // 1. compute g(r)
     for(i=0; i<nglob; i++){
        for(j=0; j<nglob; j++) {
           if(i!=j){
           for(dist=0.0, a=0; a<3; a++) {
              // compute interatomic distance
              rr[a]=rary[3*j+a]-rary[3*i+a]; 

              // apply periodic boundary condition
              if(rr[a]>=0.5*gl[a]) rr[a]-=gl[a];
              if(rr[a]<-0.5*gl[a]) rr[a]+=gl[a];

              dist+=rr[a]*rr[a];
            }
            // get interatomic distance for i-j pair
            dist=sqrt(dist);
            // if dist is less than the cutoff, increment histgram array
            if(dist<RC) {
               b = (int) (dist/DRC); 
               histgram[b]+=1.0;
            }
          }}
     }

     // 2. compute s(q)
     for(i=-nkvec; i<nkvec; i++) {
     for(j=-nkvec; j<nkvec; j++) {
     for(k=-nkvec; k<nkvec; k++) {
        qq[0]=dk[0]*i*5; 
        qq[1]=dk[1]*j*5; 
        qq[2]=dk[2]*k*5; 
        cosine=0.0; sine=0.0;
        for(a=0; a<nglob; a++) {
           rr[0]=rary[3*a];
           rr[1]=rary[3*a+1]; 
           rr[2]=rary[3*a+2]; 
           theta=rr[0]*qq[0]+rr[1]*qq[1]+rr[2]*qq[2];
           cosine+=cos(theta); 
           sine+=sin(theta);
/*
           sq[i][j][k]=cosine*cosine+sine*sine; 
           printf("i,j,k=(%d %d %d), sq %lf cosine %lf sine %lf rr[0:2]=%lf %lf %lf, qq[0:2]=%lf %lf %lf\n", 
           i,j,k, sq[i][j][k], cosine, sine, rr[0],rr[1],rr[2],qq[0],qq[1],qq[2]);
*/
//           printf("theta %lf cosine %lf sine %lf \n", theta, cos(theta), sin(theta));
        }
        //sq[i][j][k]=cosine*cosine+sine*sine; 
        printf("%d %d %d %lf\n", i,j,k,cosine*cosine+sine*sine);
     }}}
 
     fclose(fp);
  }

  // take average
  for(a=1; a<nbin; a++)
     histgram[a]=histgram[a]/(nframe*nglob); 

  fout = fopen("sq.d","wo");
  for(i=0; i<nkvec; i++) 
  for(j=0; j<nkvec; j++) 
  for(k=0; k<nkvec; k++) {
     sq[i][j][k]=sq[i][j][k]/(nframe*nglob);
     fprintf(fout,"%d %d %d %lf\n", i,j,k,sq[i][j][k]);
  }
  fclose(fout);

  // finalize, print gr, deallocate arrays 
  fout = fopen("gr.d","wo");
  for(cn=0.0,a=1; a<nbin; a++){
     dist=a*DRC; 
     fac=4.0*PI*dist*dist*DRC*rho;
     gr=histgram[a]/fac;
     cn+=histgram[a];
     fprintf(fout,"r= %f gr= %f cn= %f\n", dist, gr, cn); 
  }
  fclose(fout);

  free(nary); 
  free(tary); 
  free(rary); 
  free(vary); 
  exit(0); 
}
