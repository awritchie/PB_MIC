#include "cpp_blas.hpp"

/* blas single precision matrix-vector multiplication prototype */
void cblas_sgemv(char trans, int m, int n, float alpha, float *a, int lda,
           float *x, int incx, float beta, float *y, int incy)
{
    return sgemv_(&trans,&m,&n,&alpha,a,&lda,x,&incx,&beta,y,&incy);
}

/* blas double precision dot product */
double cblas_ddot(int n, double *x, int incx, double *y, int incy)
{
    return ddot_(&n,x,&incx,y,&incy);
}

/* blas single precision dot product */
float cblas_sdot(int n, float *x, int incx, float *y, int incy)
{
    return sdot_(&n,x,&incx,y,&incy);
}

/* blas single prevision matrix matrix multiplication */
void cblas_sgemm(char transa, char transb, int m, int n, int k, float alpha, float *a, int lda, float *b, int ldb, float beta, float *c, int ldc)
{
    return sgemm_(&transa, &transb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
}

