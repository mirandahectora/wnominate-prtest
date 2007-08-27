/*
 * nomprob.c:  This function takes a matrix of yea positions,
 *             a matrix of nay positions, and matrix of ideal points
 *             and returns a matrix indicting the probability that
 *             each voter votes yea on each vote.
 * inputs:
 *         yea (double):   m x d matrix of yea positions
 *         no   (double):   m x d matrix of no positions
 *         ideal (double):  n x d matrix of ideal points
 *         beta (double):   "beta" parameter in NOMINATE
 *         weight (double): "Dimension" weight in NOMINATE
 *         votes (int):     number of rollcall votes (m)
 *         members (int):   number of legislators (n)
 *         dims (int):      number of dimensions (d)
 *
 *  returns:
 *         yeaProbs (double): n x m matrix of yea probabilities
 *
 */

#include "R.h"
#include "Rmath.h"

void nomprob(double *yea, double *no, double *ideal, double *beta, double *w,
             int *votes, int *members, int *dims, double *yeaProb, int *normal)
{
   int i, j, k;
   double distYea, distNo;
   const int n = *members;
   const int m = *votes;
   const int d = *dims;
   const double b = *beta;
   double (*cdf)(double,double,double,int,int);
   cdf=&pnorm; 


   if(*normal==1)
    cdf=&pnorm;
   else 
    cdf=plogis;

   for (i=0;i<m;i++) {
      for (j=0;j<n;j++) {
         distYea = 0.0;
         distNo   = 0.0;
         for (k=0;k<d;k++) {
            distYea -= w[k]*(ideal[j*d+k]-yea[i*d+k])*(ideal[j*d+k]-yea[i*d+k]);
            distNo   -= w[k]*(ideal[j*d+k]-no[i*d+k])*(ideal[j*d+k]-no[i*d+k]);
         }
         yeaProb[i*n + j] = (*cdf)( b*(exp(distYea)-exp(distNo)) , 0.0, 1.0, 1, 0 );
      }
   }
}
