#include "electrostatics.hpp"

    /* Compute parameters: Copied from APBS vpbe.c
     *
     * kappa^2 = (8 pi N_A e_c^2) I_s / (1000 eps_w k_B T)
     * kappa   = 0.325567 * I_s^{1/2}   angstroms^{-1}
     * deblen  = 1 / kappa
     *         = 3.071564378 * I_s^{1/2}   angstroms
     * \bar{kappa}^2 = eps_w * kappa^2
     * zmagic  = (4 * pi * e_c^2) / (k_B T)   (we scale the diagonal later)
     *         = 7046.528838
     */

float getKappa(t_options options)
{
    return sqrt( options.I * 1.0e-16 * ( (8.0 * M_PI * N_A * e_c * e_c)/(1000.0 * options.sdie * k_B * options.temp) ) );
}

float getPhi(float &r, float *q)
{
    return q[0] / r;
}

float getPhi_Kappa(float &r, float *q, float &kappa, float &ionR, float &atomR, float &sdie)
{
    return exp(kappa*(ionR+atomR))/(1+kappa*(ionR+atomR)) * exp(-1*kappa*r)/sdie * q[0] / r;
}

float calc_r( int size, float *A, float *B )
{
#ifdef __ICC
    float C[size];
    C[:] = A[0:size] - B[0:size];
    return pow( cblas_sdot(3, &C[0], 1, &C[0], 1), 0.5);
#else
    // cblas_sdot doesn't work on my personal machine, so I'm using this crude hack...
    double C[size];
    for (int i=0;i<size;i++)
    {
        C[i] = A[i] - B[i];
    }
    return (float)pow( cblas_ddot(3, &C[0], 1, &C[0], 1), 0.5);
#endif
}


/* Kq = Rv (1.22)
 * S_ij = 1/(s_i - s_j)
 * S_ii = C * (4 pi/a_i)**.5 * f_i^shape, C ~ 1.06, a_i = area of point i, f_i = often omitted shape factor
 *
 * D_ij = -dot(n_j, s_j - s_i) * (s_j - s_i)**(-3), n_j = outward pointing surface normal vector
 * D_ii = -(2 pi + sum_(k!=i)[ D_ik * a_k] / a_i), a_i = area of point i
 *
 * A_ii = a_i
 * A_ij = 0 (?)
 *
 * Method       |   Matrix K            |   Matrix R
 * ----------------------------------------------------------
 * SS(V)PE/     |                       |
 * IEF-PCM      | [S]-(fe/2pi)[D][A][S] | -fe(1-1/(2pi)[D][A]
 *              |                       |
 * C-PCM/       |                       |
 * GCOSMO       |      [S]              | -(pdie-1)/pdie*[1]
 *              |                       |
 * DESMO        |      [S]              | -[1]+1/pdie[M]
 *
 * fe = (pdie-1)(pdie+1)
 * Mii = PhiK(s_i)/Phi0(si)
 * Mij = 0
 */
void desmo(t_options &options, std::vector<t_grid> &surface, float *K, float *R)
{
    if (options.print_method_name == 1)
    {
        std::cout << "Starting DESMO...\n";
    }
    
    const std::size_t n = surface.size();
    const std::size_t N = pow(n,2);
    
    const float C = 1.06;
    
    int NTHREADS = options.nthreads;
    #pragma omp parallel for num_threads(NTHREADS)
    for (std::size_t i=0; i<n; i++)
    {
        std::size_t ij = 0;
        float sij = 0;
        // Bottom half
        for (std::size_t j=0; j<i; j++)
        {
            ij = j + n*i;
            sij = calc_r(3,&surface[i].xyz[0], &surface[j].xyz[0]);
            K[ij] = 1. / sij;
            R[ij] = 0;
        }
        
        // Top half
        for (std::size_t j=i+1; j<n; j++)
        {
            ij = j + n*i;
            sij = calc_r(3,&surface[i].xyz[0], &surface[j].xyz[0]);
            K[ij] = 1. / sij;
            R[ij] = 0;
        }
        
        // Diagonal
        ij = i + n*i;
        K[ij] = C * pow(4 * M_PI / surface[i].a, 0.5);
        R[ij] = surface[i].PhiK / surface[i].Phi0 / options.sdie - 1.;
    }
    
    return;
}

void gcosmo(t_options &options, std::vector<t_grid> &surface, float *K, float *R)
{
    if (options.print_method_name == 1)
    {
        std::cout << "Starting GCOSMO...\n";
    }
    
    const std::size_t n = surface.size();
    const std::size_t N = pow(n,2);
    
    const float C = 1.06;
    
    int NTHREADS = options.nthreads;
    #pragma omp parallel for num_threads(NTHREADS)
    for (std::size_t i=0; i<n; i++)
    {
        std::size_t ij = 0;
        float sij = 0;
        // Bottom half
        for (std::size_t j=0; j<i; j++)
        {
            ij = j + n*i;
            sij = calc_r(3,&surface[i].xyz[0], &surface[j].xyz[0]);
            K[ij] = 1. / sij;
            R[ij] = 0;
        }
        
        // Top half
        for (std::size_t j=i+1; j<n; j++)
        {
            ij = j + n*i;
            sij = calc_r(3,&surface[i].xyz[0], &surface[j].xyz[0]);
            K[ij] = 1. / sij;
            R[ij] = 0;
        }
        
        // Diagonal
        ij = i + n*i;
        K[ij] = C * pow(4 * M_PI / surface[i].a, 0.5);
        R[ij] = ( 1 - options.sdie ) / options.sdie;
    }
    
    return;
}

void ief(t_options &options, std::vector<t_grid> &surface, float *K, float *R)
{
    if (options.print_method_name == 1)
    {
        std::cout << "Starting IEF-PCM...\n";
    }
    
    const std::size_t n = surface.size();
    const std::size_t N = pow(n,2);
    
    const float fe = (options.sdie -1)/(options.sdie + 1);
    const float C = 1.06;
    
    float *D = new float [N];
    float *A = new float [N];
    float *Dij = new float[n];
    
    
    int NTHREADS = options.nthreads;
    #pragma omp parallel for num_threads(NTHREADS)
    for (std::size_t i=0; i<n; i++)
    {
        std::size_t ij = 0;
        float sij = 0;
        Dij[i] = 0.0f;
        // Bottom half
        for (std::size_t j=0; j<i; j++)
        {
            ij = j + n*i;
            
#ifdef __ICC
            float rij[3];
            rij[:] = surface[j].xyz[:] - surface[i].xyz[:];
            D[ij] = -1 * cblas_sdot(3, &surface[j].N[0], 1, &rij[0], 1) * pow(sij,-3);
#else
            double rij[3],sN[3];
            for (int k=0; k<3; k++)
            {
                rij[k] = surface[j].xyz[k] - surface[i].xyz[k];
                sN[k] = surface[j].N[k];
            }
            D[ij] = -1 * (float)cblas_ddot(3, &sN[0], 1, &rij[0], 1) * pow(sij,-3);
#endif
            Dij[i] += D[ij] * surface[j].a;
            
            A[ij] = 0.0f;
            
            sij = calc_r(3,&surface[i].xyz[0], &surface[j].xyz[0]);
            K[ij] = 1. / sij;
            
            R[ij] = 0.0f;
        }
        
        // Top half
        for (std::size_t j=i+1; j<n; j++)
        {
            ij = j + n*i;
            
#ifdef __ICC
            float rij[3];
            rij[:] = surface[j].xyz[:] - surface[i].xyz[:];
            D[ij] = -1 * cblas_sdot(3, &surface[j].N[0], 1, &rij[0], 1) * pow(sij,-3);
#else
            double rij[3],sN[3];
            for (int k=0; k<3; k++)
            {
                rij[k] = surface[j].xyz[k] - surface[i].xyz[k];
                sN[k] = surface[j].N[k];
            }
            D[ij] = -1 * (float)cblas_ddot(3, &sN[0], 1, &rij[0], 1) * pow(sij,-3);
#endif
            Dij[i] += D[ij] * surface[j].a;
            
            A[ij] = 0.0f;
            
            sij = calc_r(3,&surface[i].xyz[0], &surface[j].xyz[0]);
            K[ij] = 1. / sij;
            
            R[ij] = 0.0f;
        }
        
        // Diagonal
        ij = i + n*i;
        D[ij] = -1. * (2. * M_PI + Dij[i]) / surface[i].a;
        A[ij] = surface[i].a;
        K[ij] = C * pow(4. * M_PI / surface[i].a, 0.5);
        R[ij] = 1.;
    }
    delete[] Dij;
    /*
     * Method       |   Matrix K            |   Matrix R
     * ----------------------------------------------------------
     * SS(V)PE/     |                       |
     * IEF-PCM      | [S]-(fe/2pi)[D][A][S] | -fe(1-1/(2pi)[D][A]
     *              |                       |
     */
    if (options.print_method_name == 1)
    {
        std::cout << "Beginning matrix-matrix multiplication...\n";
    }
    float *DA_neg_fe_ovr_2pi = new float [N];
    omp_set_num_threads(NTHREADS);
#ifdef __ICC
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n, fe / (-2. * M_PI), &D[0], n, &A[0], n, 0, &DA_neg_fe_ovr_2pi[0], n);
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n, fe / (+2. * M_PI), &D[0], n, &A[0], n, -1*fe, &R[0], n);
    delete[] D;
    delete[] A;
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1, &DA_neg_fe_ovr_2pi[0], n, &K[0], n, 1, &K[0], n);
#else
    cblas_sgemm('N', 'N', n, n, n, fe / (-2. * M_PI), &D[0], n, &A[0], n, 0, &DA_neg_fe_ovr_2pi[0], n);
    cblas_sgemm('N', 'N', n, n, n, fe / (+2. * M_PI), &D[0], n, &A[0], n, -1*fe, &R[0], n);
    delete[] D;
    delete[] A;
    cblas_sgemm('N', 'N', n, n, n, 1, &DA_neg_fe_ovr_2pi[0], n, &K[0], n, 1, &K[0], n);
#endif
    delete[] DA_neg_fe_ovr_2pi;
    
    return;
}

