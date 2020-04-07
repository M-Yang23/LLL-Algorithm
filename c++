#include <stdio.h>
#include <math.h>
#include <stdlib.h>

/*
*    Description    : this function is used to give the round interger of the input double variable
*    Input        : a double variable
*    Output        : a round interger of the given variable
*/
int roundinteger(double var)
{
    return (var > 0) ? (int)(var + 0.5) : (int)(var - 0.5);
}

/*
*    Description    : this function is used to perform inner product of two given vector or array
*    Input        : pointer of vector 'b1' and 'b2' with their length 'm'
*    Output        : pointer of the output result
*/

double innerproduct(double *b1, double *b2, int m)
{
    double sum = 0;

    while (m--)
    {
        sum += b1[m] * b2[m];
    }
    return sum;
}

/*
*    Description    : this function performs Gram-Schmidt on input Matrix
*    Input        : pointer of the input Matrix(in array format) 'bin',pointer of the output Matrix 'bout',
*                        the projection between inner-vectors,and dimension 'N*M'('N' for cloumn,'M' for row)
*    Output        : the GSO(Gram-Schmidt Othogonal) Matrix pointer
*/

void GramSchmidt(double *bin, double *bout, double *u, double *B, int M, int N)
{
    int i, j, k;
    double uTemp, projection;

    /* copy bin to bout */
    for (j = 0; j < M*N; j++)
    {
        bout[j] = bin[j];
    }

    /* Euclid Square of the first column */
    B[0] = innerproduct(bout, bout, M);

    for (i = 0; i < N; i++)
    {
        for (k = 0; k < i; k++)
        {
            projection = innerproduct(bout + k*M, bin + i*M, M);
            uTemp = projection / B[k];
            u[(i - 1)*(N - 1) + k] = uTemp;

            for (j = 0; j < M; j++)
            {
                bout[i*M + j] -= bout[k*M + j] * uTemp;
            }
        }
        /* calculate the next B[i]*/
        B[i] = innerproduct(bout + i*M, bout + i*M, M);
    }
}

/*
*    Description    : this function is used the do the reduction procedure
*    Input        : the column vector 'b',
*    Output        : update the ...
*/
void reduction(double *b, double *u, int k, int l, int M, int N)
{
    int j, i, r;
    double uTemp, uAbs;

    uTemp = u[(k - 1)*(N - 1) + l];
    uAbs = fabs(*u);
    //uAbs = (uTemp >= 0) ? uTemp : (-uTemp);
    //uAbs = (uTemp > 0) ? uTemp : (-uTemp);

    if (uAbs > 0.5)
    {
        r = roundinteger(uTemp);

        /* update u(k,k-1) <= u(k,k-1) - r */
        u[(k - 1)*(N - 1) + l] -= r;

        /* update b(k) <= b(k) -r*b(k-1) */
        for (i = 0; i < M; i++)
        {
            b[k*M + i] -= r*b[l*M + i];
        }

        /* for u(k,j) with j<k-1, update u(k,j) <= u(k,j) - r*u(k-1,j) */
        for (j = 0; j <= l - 1; j++)
        {
            u[(k - 1)*(N - 1) + j] -= r*u[(l - 1)*(N - 1) + j];
        }
    }
}


/*
*    Description    : this function is used the do the check procedure
*    Input        : 'B' the input column Euclid Squre, 'b' the input Matrix element
'u' the projection , 'k' the current subscript  and dimension 'M*N'
*    Output        : update the ...
*/
void check(double *B, double *b, double *u, int k, int l, int M, int N)
{
    int i, j;
    double uTemp, Btemp;
    double tmp;

    uTemp = u[(k - 1)*(N - 1) + k - 1];        /* this is u(k,k-1) */

    /* c(k-1) <= b(k)
    * c(k)   <= b(k-1)
    * c(i)   <= b(i) for i != k,k-1
    */
    Btemp = B[k] + uTemp * uTemp * B[k - 1];    // Btemp = c*(k-1)

    /* update u(k,k-1) <= u(k,k-1)B[k-1]/C[k-1] £¬ this is v(k,k-1) */
    u[(k - 1)*(N - 1) + k - 1] = uTemp*B[k - 1] / Btemp;

    /* update B[k] <= C[k] */
    B[k] = B[k - 1] * B[k] / Btemp;
    /* update B[k-1] <= C[k-1] */
    B[k - 1] = Btemp;

    /* exchange b[k] <=> b[k] */
    for (j = 0; j < M; j++)
    {
        tmp = b[k*M + j];
        b[k*M + j] = b[(k - 1)*M + j];
        b[(k - 1)*M + j] = tmp;
    }

    /* for j<k-1 ,i.e  j = [1 to k-2]
    * u(k-1,j) <= u(k,j)
    * u(k,j)   <= u(k-1,j)
    */
    for (j = 0; j < k - 1; j++)
    {
        tmp = u[(k - 2)*(N - 1) + j];                            // u(k-1,j)
        u[(k - 2)*(N - 1) + j] = u[(k - 1)*(N - 1) + j];
        u[(k - 1)*(N - 1) + j] = tmp;                            // u(k,j)
    }

    /* for j>k
    *
    * u(i,k)    <= u(i,k-1) - u(i,k)*u(k,k-1)
    * u(i,k-1) <= u(k,k-1) -uTemp*u(i,k)
    */

    for (i = k + 1; i < N; i++)
    {
        tmp = u[(i - 1)*(N - 1) + k];                        // u(i,k)
        /* v(i,k) <= u(i,k-1) - u(i,k)*u(k,k-1) */
        u[(i - 1)*(N - 1) + k] = u[(i - 1)*(N - 1) + (k - 1)] - u[(i - 1)*(N - 1) + k] * uTemp;
        /* v(i,k-1)  <= u(i,k-1)*v() */                        //¸üÐÂv(i,k-1)£¬¿´ÎÄÏ×Page521Ò³
        u[(i - 1)*(N - 1) + k - 1] = u[(k - 1)*(N - 1) + (k - 1)] * u[(i - 1)*(N - 1) + (k - 1)] + tmp*(1 - uTemp*u[(k - 1)*(N - 1) + (k - 1)]);
    }

}


/*
*    Description    : LLL (Lattice Reduction Alogorithm
*    Input        : the input Matrix in array 'bin',the dimension of the Matrix 'N*M',the output 'HLR'
*/
double RLLL(double *bin, int M, int N, double *HLR)
{
    int i, k, l;

    double *u;
    double *B;
    double *H;
    
    u = (double *)calloc(N*M, sizeof(double));
    B = (double *)calloc(N, sizeof(double));
    H = (double *)calloc(N*M, sizeof(double));

    for (i = 0; i < M*N; i++)
    {
        H[i] = bin[i];
    }


    GramSchmidt(H, HLR, u, B, M, N);

    k = 1;

    while (1)
    {
        l = k - 1;
        reduction(H, u, k, l, M, N);

        /* iteration procedure */
        if (B[k]<(0.75 - u[(k - 1)*(N - 1) + k - 1] * u[(k - 1)*(N - 1) + k - 1])*B[k - 1])
        {
            check(B, H, u, k, l, M, N);
            if (k>1)
            {
                k--;
            }
        }
        else
        {
            for (l = k - 2; l >= 0; l--)
            {
                reduction(H, u, k, l, M, N);
            }

            if (k == N - 1)    break;
            k++;
        }
    }


    for (i = 0; i < M*N; i++)
    {
        HLR[i] = H[i];
    }
    
    return *HLR;
    
/*    free(u);
    free(H);
    free(B);
*/
}

int main()
{
    int M = 30;
    int N = 30;                        

    int i, j;                    

    double Arr[M*N];               

   

                                                

    for (i = 0; i < M; ++i)              

    {                                    

        for (j = 0; j < N; ++j)          

        {                                

            Arr[i*j] = rand()/10.00;   

        }                                

    }         

    double hlr[M*N];

    RLLL(Arr, N, M, hlr);

}

