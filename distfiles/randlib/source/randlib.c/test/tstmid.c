#include "randlib.h"
#include <stdio.h>
#include <math.h>
#include <string.h>
#define min(a,b) ((a) <= (b) ? (a) : (b))

extern void pr_ans( float av, float avtr, float var, float vartr,
                    float xmin, float xmax);
extern void statistics(float *x,long n,float *av,
                 float *var,float *xmin,float *xmax);
extern void trstatistics(char pr_type[4],float *parin,float *av,float *var);
extern void pr_ln_ar(long *array, long larray);
extern void pr_fl_ar(float *array, long larray);

void statistics(float *x,long n,float *av,float *var,
                 float *xmin,float *xmax)
/*
**********************************************************************
     SUBROUTINE STATISTICS( X, N, AV, VAR)
               compute STATistics
                              Function
     Computes AVerage and VARiance of array X(N).
**********************************************************************
*/
{
static long i;
static float sum;

    *xmin = *xmax = *x;
    sum = 0.0;
    for(i=0; i<n; i++) {
        sum += *(x+i);
        if(*(x+i) < *xmin) *xmin = *(x+i);
        if(*(x+i) > *xmax) *xmax = *(x+i);
    }
    *av = sum/(float)n;
    sum = 0.0;
    for(i=0; i<n; i++) sum += pow(*(x+i)-*av,2.0);
    *var = sum/(float)(n-1);
}
void main(int argc,char* argv[])
/*
     Interactive test for PHRTSD
*/
{
#define mxwh 15L
#define MXNCAT 100L
static long K1 = 1;
static long i,is1,is2,itmp,iwhich,mxint,nperm,nrep,ntot,perm[500],ntry,ncat;
static float av,avtr,var,vartr,xmin,xmax,array[1000],param[3],pevt,psum,rtry,
             prob[MXNCAT];
static char pr_type[4],phrase[100];
static long iarray[1000];

    puts(" Tests most generators of specific distributions.");
    puts(" Generates 1000 deviates: reports mean and variance.");
    puts(" Also reports theoretical mean and variance.");
    puts(" If theoretical mean or var doesn't exist prints -1.");
    puts(" For permutations, generates one permutation of 1..n");
    puts("     and prints it.");
    puts(" For uniform integers asks for upper bound, number of");
    puts("     replicates per integer in 1..upper bound.");
    puts("     Prints table of num times each integer generated.");
    puts(" For multinomial asks for number of events to be");
    puts("     classified, number of categories in which they");
    puts("     are to be classified, and the probabilities that");
    puts("     an event will be classified in the categories,");
    puts("     for all but the last category.  Prints table of");
    puts("     number of events by category, true probability");
    puts("     associated with each category, and observed");
    puts("     proportion of events in each category.");
/*
     Menu for choosing tests
*/
S10:
    puts(" Enter number corresponding to choice:");
    puts("      (0) Exit this program");
    puts("      (1) Generate Chi-Square deviates");
    puts("      (2) Generate noncentral Chi-Square deviates");
    puts("      (3) Generate F deviates");
    puts("      (4) Generate noncentral F  deviates");
    puts("      (5) Generate random permutation");
    puts("      (6) Generate uniform integers");
    puts("      (7) Generate uniform reals");
    puts("      (8) Generate beta deviates");
    puts("      (9) Generate binomial outcomes");
    puts("     (10) Generate Poisson outcomes");
    puts("     (11) Generate Exponential deviates");
    puts("     (12) Generate gamma deviates");
    puts("     (13) Generate multinomial outcomes");
    puts("     (14) Generate normal deviates");
    puts("     (15) Generate negative binomial outcomes");
    scanf("%ld",&iwhich);
    if(!(iwhich < 0 || iwhich > mxwh)) goto S20;
    printf(" Choices are 1..%12ld - try again.\n",mxwh);
    goto S10;
S20:
    if(iwhich == 0) {
        puts(" Normal termination rn tests");
        exit(0);
      }
    puts(" Enter phrase to initialize rn generator");
    scanf("%s",phrase);
    phrtsd(phrase,&is1,&is2);
    setall(is1,is2);
    if(1 != iwhich) goto S40;
/*
     Chi-square deviates
*/
    strcpy(pr_type,"chis");
    puts(" Enter (real) df for the chi-square generation");
    scanf("%f",param);
    for(i=0; i<1000; i++) *(array+i) = genchi(*param);
    statistics(array,1000L,&av,&var,&xmin,&xmax);
    trstatistics(pr_type,param,&avtr,&vartr);
    pr_ans( av, avtr, var, vartr, xmin, xmax );
    goto S320;
S40:
    if(2 != iwhich) goto S60;
/*
     Noncentral Chi-square deviates
*/
    strcpy( pr_type, "ncch");
    puts(" Enter (real) df");
    puts("       (real) noncentrality parameter");
    scanf("%f%f",param,(param+1));
    for(i=0; i<1000; i++) *(array+i) = gennch(*param,*(param+1));
    statistics(array,1000L,&av,&var,&xmin,&xmax);
    trstatistics(pr_type,param,&avtr,&vartr);
    pr_ans( av, avtr, var, vartr, xmin, xmax );
    goto S320;
S60:
    if(3 != iwhich) goto S80;
/*
     F deviates
*/
    strcpy( pr_type, "f   ");
    puts(" Enter (real) df of the numerator");
    puts("       (real) df of the denominator");
    scanf("%f%f",param,(param+1));
    for(i=0; i<1000; i++) *(array+i) = genf(*param,*(param+1));
    statistics(array,1000L,&av,&var,&xmin,&xmax);
    trstatistics(pr_type,param,&avtr,&vartr);
    pr_ans( av, avtr, var, vartr, xmin, xmax );
    goto S320;
S80:
    if(4 != iwhich) goto S100;
/*
     Noncentral F deviates
*/
    strcpy( pr_type, "ncf ");
    puts(" Enter (real) df of the numerator");
    puts("       (real) df of the denominator");
    puts("       (real) noncentrality parameter");
    scanf("%f%f%f",param,(param+1),(param+2));
    for(i=0; i<1000; i++) *(array+i) = gennf(*param,*(param+1),*(param+2));
    statistics(array,1000L,&av,&var,&xmin,&xmax);
    trstatistics(pr_type,param,&avtr,&vartr);
    pr_ans( av, avtr, var, vartr, xmin, xmax );
    goto S320;
S100:
    if(5 != iwhich) goto S140;
/*
     Random permutation
*/
S110:
    puts(" Enter size of permutation - range 1 to 500");
    scanf("%ld",&nperm);
    if(!(nperm < 1 || nperm > 500)) goto S120;
    printf(" Permutation size must be between 1 and 500 - try again!\n");
    goto S110;
S120:
    printf("       Random Permutation Generated - Size%12ld\n",nperm);
    for(i=1; i<=nperm; i++) *(perm+i-1) = i;
    genprm(perm,nperm);
    puts(" Perm Generated");
    pr_ln_ar(perm, nperm);
    goto S320;
S140:
    if(6 != iwhich) goto S170;
/*
     Uniform integer
*/
    puts(" Enter maximum uniform integer - upper limit 1000");
    scanf("%ld",&mxint);
    if (mxint < 1 || mxint > 1000){
        puts(" Maximum integer must be between 1 and 1000");
        goto S140;
      }
    puts(" Enter number of replications per integer");
    scanf("%ld",&nrep);
    for(i=1; i<=1000; i++) *(iarray+i-1) = 0;
    ntot = mxint*nrep;
    for(i=1; i<=ntot; i++) {
        itmp = ignuin(K1,mxint);
        *(iarray+itmp-1) += 1;
    }
    puts("         Counts of Integers Generated");
    pr_ln_ar(iarray,mxint);
    goto S320;
S170:
    if(7 != iwhich) goto S190;
/*
     Uniform real
*/
    strcpy( pr_type, "unif");
    puts(" Enter Low then High bound for uniforms");
    scanf("%f%f",param,(param+1));
    for(i=1; i<=1000; i++) *(array+i-1) = genunf(*param,*(param+1));
    statistics(array,1000L,&av,&var,&xmin,&xmax);
    trstatistics(pr_type,param,&avtr,&vartr);
    pr_ans( av, avtr, var, vartr, xmin, xmax );
    goto S320;
S190:
    if(8 != iwhich) goto S210;
/*
     Beta deviate
*/
    strcpy( pr_type, "beta");
    puts(" Enter A, B for Beta deviate");
    scanf("%f%f",param,(param+1));
    for(i=1; i<=1000; i++) *(array+i-1) = genbet(*param,*(param+1));
    statistics(array,1000L,&av,&var,&xmin,&xmax);
    trstatistics(pr_type,param,&avtr,&vartr);
    pr_ans( av, avtr, var, vartr, xmin, xmax );
    goto S320;
S210:
    if(9 != iwhich) goto S240;
/*
     Binomial outcomes
*/
    strcpy( pr_type, "bin ");
    printf(" Enter number of trials, Prob event for binomial outcomes\n");
    scanf("%ld%f",&ntry,&pevt);
    for(i=1; i<=1000; i++) *(iarray+i-1) = ignbin(ntry,pevt);
    for(i=1; i<=1000; i++) *(array+i-1) = *(iarray+i-1);
    statistics(array,1000L,&av,&var,&xmin,&xmax);
    *param = ntry;
    *(param+1) = pevt;
    trstatistics(pr_type,param,&avtr,&vartr);
    pr_ans( av, avtr, var, vartr, xmin, xmax );
    goto S320;
S240:
    if(10 != iwhich) goto S250;
/*
     Poisson outcomes
*/
    strcpy( pr_type, "pois");
    puts(" Enter mean for Poisson generation");
    scanf("%f",param);
    for(i=1; i<=1000; i++) *(iarray+i-1) = ignpoi(*param);
    for(i=1; i<=1000; i++) *(array+i-1) = *(iarray+i-1);
    statistics(array,1000L,&av,&var,&xmin,&xmax);
    trstatistics(pr_type,param,&avtr,&vartr);
    pr_ans( av, avtr, var, vartr, xmin, xmax );
    goto S320;
S250:
    if(11 != iwhich) goto S270;
/*
     Exponential outcomes
*/
    strcpy( pr_type, "exp ");
    puts(" Enter mean for Exponential generation");
    scanf("%f",param);
    for(i=1; i<=1000; i++) *(array+i-1) = genexp(*param);
    statistics(array,1000L,&av,&var,&xmin,&xmax);
    trstatistics(pr_type,param,&avtr,&vartr);
    pr_ans( av, avtr, var, vartr, xmin, xmax );
    goto S320;

S270:
    if(12 != iwhich) goto S280;
/*
     Gamma deviates
*/
    strcpy( pr_type, "gamm");
    puts(" Enter (real) A, (real) R for Gamma deviate");
    scanf("%f%f",param,(param+1));
    for(i=1; i<=1000; i++) *(array+i-1) = gengam(*param, *(param+1));
    statistics(array,1000L,&av,&var,&xmin,&xmax);
    trstatistics(pr_type,param,&avtr,&vartr);
    pr_ans( av, avtr, var, vartr, xmin, xmax );
    goto S320;

S280:
    if(13 != iwhich) goto S290;
/*
     Multinomial outcomes
*/
    puts(" Enter (int) number of observations: ");
    scanf("%ld",&ntry);
    printf(" Enter (int) num. of categories: <= %ld\n",MXNCAT);
    scanf("%ld",&ncat);
    while(ncat > MXNCAT) {
      printf(" number of categories must be <= %ld\n",MXNCAT);
      puts(" Try again ... ");
      printf(" Enter (int) num. of categories: <= %ld\n",MXNCAT);
      scanf("%ld",&ncat);
    }
    printf(" Enter (real) prob. vector of length %ld\n",ncat-1L);
    for(i=1L;i<ncat;i++) scanf("%f",prob + i - 1);
    genmul(ntry,prob,ncat,iarray);
    ntot = 0L;
    if (ntry > 0L) {
      rtry = (float) ntry;
      for(i=0L;i<ncat;i++) {
	ntot += *(iarray + i);
	*(array + i) = *(iarray + i) / rtry;
      }
    } else {
      for(i=0L;i<ncat;i++) {
	ntot += *(iarray + i);
	*(array + i) = 0.0;
      }
    }
    psum = 0.0;
    for(i=1L;i<ncat;i++) psum += *(prob + i - 1);
    *(prob + ncat - 1) = 1.0 - psum;
    printf(" Total number of observations: %ld\n",ntot);
    puts(" Total observations by category: ");
    pr_ln_ar(iarray,ncat);
    puts(" True probabilities by category: ");
    pr_fl_ar(prob,ncat);
    puts(" Observed proportions by category: ");
    pr_fl_ar(array,ncat);
    goto S320;

S290:
    if(14 != iwhich) goto S300;
/*
     Normal deviates
*/
    strcpy( pr_type, "norm");
    puts(" Enter (real) AV, (real) SD for Normal");
    scanf("%f%f",param,(param+1));
    for(i=1; i<=1000; i++) *(array+i-1) = gennor(*param, *(param+1));
    statistics(array,1000L,&av,&var,&xmin,&xmax);
    trstatistics(pr_type,param,&avtr,&vartr);
    pr_ans( av, avtr, var, vartr, xmin, xmax );
    goto S320;

S300:
    if(15 != iwhich) goto S310;
/*
     Negative Binomial outcomes
*/
    strcpy( pr_type, "nbin");
    puts(" Enter required (int) Number of events, then");
    puts(" (real) Prob event for negative binomial outcomes");
    scanf("%ld%f",&ntry,&pevt);
    for(i=1; i<=1000; i++) *(iarray+i-1) = ignnbn(ntry,pevt);
    for(i=1; i<=1000; i++) *(array+i-1) = *(iarray+i-1);
    statistics(array,1000L,&av,&var,&xmin,&xmax);
    *param = ntry;
    *(param+1) = pevt;
    trstatistics(pr_type,param,&avtr,&vartr);
    pr_ans( av, avtr, var, vartr, xmin, xmax );
    goto S320;

S310:
S320:
    goto S10;
#undef mxwh
#undef MXNCAT
}

void trstatistics(char pr_type[4],float *parin,float *av,float *var)
/*
**********************************************************************
     SUBROUTINE TRSTATISTICS( TYPE, PARIN, AV, VAR )
               TRue STATistics
     Returns mean and variance for a number of statistical distribution
     as a function of their parameters.
                              Arguments
     TYPE --> Character string indicating type of distribution
             'chis' chisquare
             'ncch' noncentral chisquare
             'f   ' F (variance ratio)
             'ncf ' noncentral f
             'unif' uniform
             'beta' beta distribution
             'bin ' binomial
             'pois' poisson
             'exp ' exponential
             'gamm' gamma
             'norm' normal
             'nbin' negative binomial
                         CHARACTER*(4) TYPE
     PARIN --> Array containing parameters of distribution
              chisquare
               PARIN(1) is df
              noncentral chisquare
               PARIN(1) is df
               PARIN(2) is noncentrality parameter
              F (variance ratio)
               PARIN(1) is df numerator
               PARIN(2) is df denominator
              noncentral F
               PARIN(1) is df numerator
               PARIN(2) is df denominator
               PARIN(3) is noncentrality parameter
              uniform
               PARIN(1) is LOW bound
               PARIN(2) is HIGH bound
              beta
               PARIN(1) is A
               PARIN(2) is B
              binomial
               PARIN(1) is Number of trials
               PARIN(2) is Prob Event at Each Trial
              poisson
               PARIN(1) is Mean
              exponential
               PARIN(1) is Mean
              gamma
               PARIN(1) is A
               PARIN(2) is R
              normal
               PARIN(1) is Mean
               PARIN(2) is Standard Deviation
              negative binomial
               PARIN(1) is required Number of events
               PARIN(2) is Probability of event
                         REAL PARIN(*)

     AV <-- Mean of specified distribution with specified parameters
                         REAL AV
     VAR <-- Variance of specified distribution with specified paramete
                         REAL VAR
                              Note
     AV and Var will be returned -1 if mean or variance is infinite
*/
{
static float a,b,range;

    if(strcmp("chis",pr_type) !=0 ) goto S10;
    *av = *parin;
    *var = 2.0**parin;
    goto S210;
S10:
    if(strcmp("ncch",pr_type) != 0) goto S20;
    a = *parin+*(parin+1);
    b = *(parin+1)/a;
    *av = a;
    *var = 2.0*a*(1.0+b);
    goto S210;
S20:
    if(strcmp("f   ",pr_type) != 0) goto S70;
    if(!(*(parin+1) <= 2.0001)) goto S30;
    *av = -1.0;
    goto S40;
S30:
    *av = *(parin+1)/(*(parin+1)-2.0);
S40:
    if(!(*(parin+1) <= 4.0001)) goto S50;
    *var = -1.0;
    goto S60;
S50:
    *var = 2.0*pow(*(parin+1),2.0)*(*parin+*(parin+1)-2.0)/(*parin*pow(*(parin
      +1)-2.0,2.0)*(*(parin+1)-4.0));
S60:
    goto S210;
S70:
    if(strcmp("ncf ",pr_type) != 0) goto S120;
    if(!(*(parin+1) <= 2.0001)) goto S80;
    *av = -1.0;
    goto S90;
S80:
    *av = *(parin+1)*(*parin+*(parin+2))/((*(parin+1)-2.0)**parin);
S90:
    if(!(*(parin+1) <= 4.0001)) goto S100;
    *var = -1.0;
    goto S110;
S100:
    a = pow(*parin+*(parin+2),2.0)+(*parin+2.0**(parin+2))*(*(parin+1)-2.0);
    b = pow(*(parin+1)-2.0,2.0)*(*(parin+1)-4.0);
    *var = 2.0*pow(*(parin+1)/ *parin,2.0)*(a/b);
S110:
    goto S210;
S120:
    if(strcmp("unif",pr_type) != 0) goto S130;
    range = *(parin+1)-*parin;
    *av = *parin+range/2.0;
    *var = pow(range,2.0)/12.0;
    goto S210;
S130:
    if(strcmp("beta",pr_type) != 0) goto S140;
    *av = *parin/(*parin+*(parin+1));
    *var = *av**(parin+1)/((*parin+*(parin+1))*(*parin+*(parin+1)+1.0));
    goto S210;
S140:
    if(strcmp("bin ",pr_type) != 0) goto S150;
    *av = *parin**(parin+1);
    *var = *av*(1.0-*(parin+1));
    goto S210;
S150:
    if(strcmp("pois",pr_type) != 0) goto S160;
    *av = *var = *parin;
    goto S210;
S160:
    if(strcmp("exp ",pr_type) != 0) goto S170;
    *av = *parin;
    *var = *av * *av;
    goto S210;
S170:
    if(strcmp("gamm",pr_type) != 0) goto S180;
    *av = *(parin + 1) / *parin;
    *var = *av / *parin;
    goto S210;
S180:
    if(strcmp("norm",pr_type) != 0) goto S190;
    *av = *parin;
    *var = *(parin + 1) * *(parin + 1);
    goto S210;
S190:
    if(strcmp("nbin",pr_type) != 0) goto S200;
    *av = *parin * (1.0 - *(parin + 1)) / *(parin + 1);
    *var = *av / *(parin + 1);
    goto S210;
S200:
    puts("Unimplemented pr_type in TRSTATISTICS");
    exit(1);
S210:
    return;
}

void pr_ans( float av, float avtr, float var, float vartr,
             float xmin, float xmax) 
{

    printf(" Mean Generated          %15.7g     True          %15.7g\n",
           (double)av, (double)avtr);
    printf(" Variance Gener          %15.7g     True          %15.7g\n",
           (double)var, (double)vartr);
    printf(" Minimum                 %15.7g     Maximum       %15.7g\n",
           (double)xmin, (double)xmax);
}

void pr_ln_ar(long *array, long larray)
/*
  Print array of longs; length of array is larray.  Print 10 to a
  line using 7 characters for each field.
*/
{
    int lo, hi, i, up;
    up = larray-1;
    for (lo=0; lo <=up; lo +=10) {
        hi = min(lo+9,up);
    	for (i=lo; i<=hi; i++)
            printf(" %7ld",*(array+i));
    puts("");    
    }
}    

void pr_fl_ar(float *array, long larray)
/*
  Print array of floats; length of array is larray.  Print 8 to a
  line using 9 characters for each field.
*/
{
    int lo, hi, i, up;
    up = larray-1;
    for (lo=0; lo <=up; lo +=8) {
        hi = min(lo+7,up);
    	for (i=lo; i<=hi; i++)
            printf(" %.7f",*(array+i));
    puts("");    
    }
}    
