#include "mex.h"
#include "matrix.h"
/*
 * aligns s1 to s2 via simple-gap NW global alignment with ends-free
 * gapping. returns the substitutions from s1 to s2 in this alignment via a
 * vector of positions along s1 and nts
 *
 * motivation:
 * nwalign, which accompanies MATLAB bioinfo package 3.7 has fast
 * alignments, but all the type-checking and conversion to and from nt to
 * int are slow, adding 2-3 orders of magnitude of time when performing
 * short alignments.
 *
 * in:
 * s1,s2: uint8, sequence vectors (a=1 c=2 g=3 t=4)
 * m: double, score matrix
 * g: double, gap penalty (open and extension)
 * h: double, homopolymer gap penalty
 * out:
 * pos: uint16, positions of substitutions on alignment of s1 to s2
 * nt: uint8, identities of substitutions (a=0 c=1 g=2 t=3)
 */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if(nrhs != 5)
        mexErrMsgTxt("This function requires 5 inputs: 2 uint8 sequence vectors, 1 double score matrix, 1 double gap penalty, 1 double homopolymer gap penalty");
    else if(nlhs > 2)
        mxErrMsgTxt("Too many output arguments.");
    else if(!mxIsUint8(prhs[0]) || !mxIsUint8(prhs[1]))
        mxErrMsgTxt("Input sequences must be of type uint8");
    else if(!mxIsDouble(prhs[2]))
        mxErrMsgTxt("Input score matrix must be of type double");
    else if(!mxIsDouble(prhs[3]) || !mxIsDouble(prhs[4]))
        mxErrMsgTxt("Input gap penality scores must be of type double");
    else if(mxGetM(prhs[0]) != 1 || mxGetM(prhs[1]) != 1)
        mxErrMsgTxt("Input sequences must be row vectors");
    else if(mxGetM(prhs[2]) != 4 || mxGetM(prhs[2]) != 4)
        mxErrMsgTxt("Score matrix must be 4x4");
    
    mwSize i,j,s,l1,l2,nsubs=0;
    l1 = mxGetN(prhs[0]);
    l2 = mxGetN(prhs[1]);
    unsigned char *s1, *s2, al[l1];
    
    s1 = mxGetData(prhs[0]);
    s2 = mxGetData(prhs[1]);
    
    int P[l1+1][l2+1],homo[l1];
    /*find locations where s1 has homopolymer and put 1s in "homo" vector*/
    for (s=0, i=0; i<l1; i++) {
        if (i == l1 - 1 || s1[i] != s1[i+1]) {
            if (i - s > 1)
                for (j=s; j<=i; j++)
                    homo[j] = 1;
            else
                for (j=s; j<=i; j++)
                    homo[j] = 0;
            s = i+1;
        }
    }
    
    double D[l1+1][l2+1],m[4][4],diag,left,up,
            G=mxGetScalar(prhs[3]),GH=mxGetScalar(prhs[4]);
    double *mp = mxGetPr(prhs[2]);
    for (i=0; i<4; i++) {
        for (j=0; j<4; j++) {
            m[j][i] = *(mp++); /*matrix is stored column-wise*/
        }
    }
    
    plhs[0] = mxCreateNumericMatrix(1,l1,mxINT16_CLASS,mxREAL);
    plhs[1] = mxCreateNumericMatrix(1,l1,mxUINT8_CLASS,mxREAL);
    
    short *pos = mxGetData(plhs[0]);
    unsigned char *nt = mxGetData(plhs[1]);
    /*mwIndex msub[2]; /*vector for holding subscripts when indexing m*/
    /*Fill out DP matrix, D, and pointback matrix, P.
     *P has 1 for diag, 2 for left, and 3 for up
     */
    for (i=0; i<=l1; i++) {
        D[i][0] = 0;
        P[i][0] = 3;
    }
    for (j=0; j<=l2; j++) {
        D[0][j] = 0;
        P[0][j] = 2;
    }
    for (i = 1; i <= l1; i++) {
        for (j = 1; j <= l2; j++) {
            if (j == l2)
                up = D[i-1][j];
            else if (homo[i-1])
                up = D[i-1][j] + GH;
            else
                up = D[i-1][j] + G;
            
            if (i == l1)
                left = D[i][j-1];
            else if (homo[i-1] && s2[j-1] == s1[i-1])
                left = D[i][j-1] + GH;
            else
                left = D[i][j-1] + G;
            
            diag = D[i-1][j-1] + m[s1[i-1] - 1][s2[j-1] - 1];
            
            /* note: the following comparisons match the MATLAB nwalign
             * conventions about what to do in cases of ties. the
             * preference ordering of ties is diag > up > left whereas in
             * older MATLAB releases it was diag > left > up.
             */
            if (up >= diag && up >= left) {
                D[i][j] = up;
                P[i][j] = 3;
                
            } else if (left >= diag) {
                D[i][j] = left;
                P[i][j] = 2;
                
            } else {
                D[i][j] = diag;
                P[i][j] = 1;
            }
        }
    }

    /*traceback over DP matrix and fill out al, */
    i = l1;
    j = l2;
    while ( i > 0 ) {
        switch ( P[i][j] ) {
            case 1:
                al[--i] = s2[--j];
                break;
            case 2:
                j--;
                break;
            case 3:
                al[--i] = 5;
                break;
        }
    }
    
    /*traverse the alignment and look for substitutions*/
    for (i = 0; i < l1; i++) {
        if ( ( s1[i] != al[i] ) && ( s1[i] != 5 ) && ( al[i] != 5 ) )
        {
            pos[nsubs] = i + 1;
            nt[nsubs] = al[i];
            nsubs++;
        }
    }
    /*resize the outputs to be 1 by nsubs*/
    mwSize dims[2];
    dims[0] = 1;
    dims[1] = nsubs;
    mxSetDimensions(plhs[0],dims,2);
    mxSetDimensions(plhs[1],dims,2);
    /*reallocate the memory to contain only nsubs entries*/
    mxSetData(plhs[0],mxRealloc(mxGetData(plhs[0]),nsubs*mxGetElementSize(plhs[0])));
    mxSetData(plhs[1],mxRealloc(mxGetData(plhs[1]),nsubs*mxGetElementSize(plhs[1])));
    return;
}