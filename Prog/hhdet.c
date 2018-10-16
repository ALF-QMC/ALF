#include <plasma.h>
#include <plasma_descriptor.h>
#include <plasma_types.h>

#define T(m, n) (plasma_complex64_t*)plasma_tile_addr(*T, m, n)

/** Calculate the determinant of a householder object
 */
void hhdet(plasma_desc_t* T, plasma_complex64_t* phase, int* nvar, int* N_size)
{
    /* PLASMA matrices have a tile based layout that consists of gmt * gnt tiles.
     * We ignore a possible submatrix structure that plasma supports.
     * On top of that the QR algorithms that is applied to each tile
     * has inner blocking that has to be considered.
     */
    int upper = (T->gnt < T->gmt) ? T->gnt : T->gmt; //min(T->gnt, T->gmt)
    int ntausleft = *N_size;
    for (int i = 0; i < upper; ++i)//loop over the diagonal tiles
    {
        for (int b = 0; b < (T->nb / T->mb); ++b )//loop over the blocks used by the QR decomposition
        {//Obviously we also assume that a tile has dimensions that are multiples of the inner block size...
            int ul = (T->mb < ntausleft)? T->mb : ntausleft; // min(T->mb, ntausleft)
            if (ul > 0)
            {
                for (int k = 0; k < ul; ++k)
                {
                    plasma_complex64_t tau = ( T(i, i) )[b*T->mb*T->mb + k * T->mb + k];
                    if (creal(tau) != 0.0 || cimag(tau) != 0.0)
                    {
                    /*here we calculate the determinant of a single householder reflector:
                     * det(1 - tau * v v* ) = 1 - tau * v^* v
                     * In lapack the scalar tau and the vector v are scaled such that |tau|^2 |v|^2 = 2 Re(tau)
                     * The complete determinant det(Q) is the product of all reflectors. See http://www.netlib.org/lapack/lug/node128.html
                     * In the case of the triangular factors T the scalars tau can be recovered from their main diagonal
                     */
                        if(*nvar == 1) tau = conj(tau);
                        double x = cabs(tau);
                        plasma_complex64_t z = 1.0 - 2.0*(tau/x) * (creal(tau)/x);
                        (*phase) *= z/ cabs(z) ;
                    }
                }
                ntausleft -= T->mb;//ntausleft now can become negative
            }
        }
    }
}

#include <stdio.h>
/**
 * 
 * @param maxj
 * @param maxi
 * @param nb corresponds to the leading dimension of the tile is either the tile block width for full tiles or the remainder for border tiles
 * */
void tile_kernel(plasma_complex64_t* tpup, plasma_complex64_t* rhs, plasma_complex64_t* dl, plasma_complex64_t* dr, plasma_complex64_t* dup, int maxj, int maxi, const int nb)
{
    for(int ji = 0; ji < maxj; ++ji)
            {
                for(int ii = 0; ii < maxi; ++ii)
                {
                    if (creal(dl[ji]) <= 1.0)
                    {
                        plasma_complex64_t dlj = dl[ji];
                        if (creal( dr[ii] ) <= 1.0)
                            tpup[ji*nb + ii] = rhs[ji*nb + ii] + dr[ii] * dl[ji] * tpup[ji*nb + ii];
                        else
                            tpup[ji*nb + ii] = dup[ii] * rhs[ji*nb + ii] + dlj * tpup[ji*nb + ii];
                    }
                    else
                    {
                        plasma_complex64_t dlj = 1.0/dl[ji];
                        if (creal( dr[ii] ) <= 1.0)
                            tpup[ji*nb + ii] = dlj * rhs[ji*nb + ii] + dup[ii] * tpup[ji*nb + ii];
                        else
                            tpup[ji*nb + ii] = rhs[ji*nb + ii]/dr[ii]/dl[ji] + tpup[ji*nb + ii];
                    }
                }
            }
}

#define TPUP(m, n) (plasma_complex64_t*)plasma_tile_addr_general(*tpup, m, n)
#define RHS(m, n) (plasma_complex64_t*)plasma_tile_addr_general(*rhs, m, n)

/* Apply left and right scales.
 * tpup are assumed to be identical square matrices
 */
void applylrscales(plasma_complex64_t* dl, plasma_complex64_t* dr, plasma_complex64_t* dup, plasma_desc_t* tpup, plasma_desc_t* rhs, const int nsize)
{
    /* Original Fortran Code as reference
     *  DO J = 1,N_size
        DO I = 1,N_size
          If( dble(udvl%D(J))<=1.d0) then
            DLJ=udvl%D(J)
              If( dble(udvr%D(I))<=1.d0 ) then
                TPUP(I,J) = RHS(I,J)+udvr%D(I)*udvl%D(J)*TPUP(I,J)
              else
                TPUP(I,J) = DUP(I)*RHS(I,J) + DLJ*TPUP(I,J)
              endif
          else
            DLJ=1.d0/udvl%D(J)
              If( dble(udvr%D(I))<=1.d0 ) then
                TPUP(I,J) = DLJ*RHS(I,J)+DUP(I)*TPUP(I,J)
              else
                TPUP(I,J) = RHS(I,J)/udvr%D(I)/udvl%D(J)+TPUP(I,J)
              endif
          endif
        ENDDO
        ENDDO
     */

const int nb = tpup->nb;
int fb = nsize/nb;
plasma_complex64_t* tpupptr = tpup->matrix;
plasma_complex64_t* rhsptr = rhs->matrix;
const int nt = nsize/nb;
const int ld = nsize - fb*nb;
printf("ndim = %i, block size %i\n", nsize, nb);
// #pragma omp task depend(in:dup[0:nsize]) depend(inout:tpupptr[0:nsize*nsize]) depend(inout:rhsptr[0:nsize*nsize])
{
    for(int jt = 0; jt < nt; ++jt) // jt:= j-tile, it := i-tile, ji: = j-inner, ii := i-inner
    {
        for(int it = 0; it < nt; ++it)
        {
            plasma_complex64_t* trhs = (RHS(jt, it));
            plasma_complex64_t* ttpup = (TPUP(jt, it));
            printf("IJ jt = %i, it = %i\n", jt, it);
            printf("IJ maxjdim = %i\n", jt*nb + nb);
            printf("IJ maxidim = %i\n", it*nb + nb);
            tile_kernel(ttpup, trhs, dl + jt*nb, dr + it*nb, dup + it*nb, nb, nb, nb);
        }
        //remainder loop along the i-axis
        if(ld > 0)
        {
            plasma_complex64_t* trhs = (RHS(jt, fb));
            plasma_complex64_t* ttpup = (TPUP(jt, fb));
            printf("iJ jt = %i, it = %i\n", jt, fb);
            printf("iJ maxjdim = %i\n", jt*nb + nb);
            printf("iJ maxidim = %i\n", fb*nb + ld);
            tile_kernel(ttpup, trhs, dl + jt*nb, dr + fb*nb, dup + fb*nb, nb, ld, ld);
        }
    }
    //remainder loop along the j axis
    if(ld > 0)
    {
        for(int it = 0; it < nt; ++it)
        {
            plasma_complex64_t* trhs = (RHS(fb, it));
            plasma_complex64_t* ttpup = (TPUP(fb, it));
            printf("Ji jt = %i, it = %i\n", fb, it);
            printf("Ji maxjdim = %i\n", fb*nb + ld);
            printf("Ji maxidim = %i\n", it*nb + nb);
            tile_kernel(ttpup, trhs, dl + fb*nb, dr + it*nb, dup + it*nb, nb, ld, ld);
        }
        //remainder loop along the i and j axis
        plasma_complex64_t* trhs = (RHS(fb, fb));
        plasma_complex64_t* ttpup = (TPUP(fb, fb));
                printf("ji jt = %i, it = %i\n", fb, fb);
            printf("ji maxjdim = %i\n", fb*nb + ld);
            printf("ji maxidim = %i\n", fb*nb + ld);
        tile_kernel(ttpup, trhs, dl + fb*nb, dr + fb*nb, dup + fb*nb, ld, ld, ld);
    }
}
}

// for(int jt = 0; jt < nsize/nb; ++jt )
// {
//     for(int ji = 0; ji < nb; ++ji)
//     {
//         if(creal(dl[jt*nb + ji]) <= 1.0 )
//         {
//             for(int it = 0; it < fb; ++it)
//             {
//                 plasma_complex64_t* trhs = (RHS(jt, it));
//                 plasma_complex64_t* ttpup = (TPUP(jt, it));
// // #pragma omp task depend(in:trhs[ji*nb:ji*nb+nb]) \
// //                  depend(in:dup[it*nb:(it*nb+nb)]) \
// //                  depend(inout:ttpup[ji*nb:ji*nb+nb])
//                  {
//             plasma_complex64_t dlj = dl[jt*nb + ji]; // move it here for easier dependencies
//                 for(int ii = 0; ii < nb; ++ii)
//                 {
//                     if ( creal( dr[it*nb + ii] ) <= 1.0)
//                         (TPUP(jt, it))[ji*nb + ii] = (RHS(jt, it))[ji*nb + ii] + dr[it*nb + ii] * dl[jt*nb + ji] * (TPUP(jt, it))[ji*nb + ii];
//                     else
//                         (TPUP(jt, it))[ji*nb + ii] = dup[it*nb + ii] * (RHS(jt, it))[ji*nb + ii] + dlj * (TPUP(jt, it))[ji*nb + ii];
//                         
//                 }
//                  }
//             }
//             //REMAINDER
//                 plasma_complex64_t* trhs = (RHS(jt, fb));
//                 plasma_complex64_t* ttpup = (TPUP(jt, fb));
//                 int ld = nsize - fb*nb;
// // #pragma omp task depend(in:trhs[ji*ld:ji*ld + ld]) \
// //                  depend(in:dup[fb*nb:(fb*nb+ld)]) \
// //                  depend(inout:ttpup[ji*ld:ji*ld + ld])
//             {
//                 plasma_complex64_t dlj = dl[jt*nb + ji]; // move it here for easier dependencies
//                 for(int ii = 0; ii < ld; ++ii)
//                 {
//                     if ( creal( dr[fb*nb + ii] ) <= 1.0)
//                         (TPUP(jt, fb))[ji*ld + ii] = (RHS(jt, fb))[ji*ld + ii] + dr[fb*nb + ii] * dl[jt*nb + ji] * (TPUP(jt, fb))[ji*ld + ii];
//                     else
//                         (TPUP(jt, fb))[ji*ld + ii] = dup[fb*nb + ii] * (RHS(jt, fb))[ji*ld + ii] + dlj * (TPUP(jt, fb))[ji*ld + ii];
//                         
//                 }
//             }
//         }
//         else
//         {
//             for(int it = 0; it < fb; ++it)
//             {
//                 plasma_complex64_t* trhs = (RHS(jt, it));
//                 plasma_complex64_t* ttpup = (TPUP(jt, it));
// #pragma omp task depend(in:trhs[ji*nb:ji*nb+nb]) \
//                  depend(in:dup[it*nb:(it*nb+nb)]) \
//                  depend(inout:ttpup[ji*nb:ji*nb+nb])
//                  {
//             plasma_complex64_t dlj = 1.0/dl[jt*nb + ji];
//                 for(int ii = 0; ii < tpup->nb; ++ii)
//                 {
//                     if ( creal( dr[it*nb + ii] ) <= 1.0)
//                         (TPUP(jt, it))[ji*nb + ii] = dlj * (RHS(jt, it))[ji*nb + ii] + dup[it*nb + ii] * (TPUP(jt, it))[ji*nb + ii];
//                     else
//                         (TPUP(jt, it))[ji*nb + ii] = (RHS(jt, it))[ji*nb + ii]/dr[it*nb + ii]/dl[jt*nb + ji] + (TPUP(jt, it))[ji*nb + ii];
//                 }
//                  }
//             }
//             //REMAINDER
//                 plasma_complex64_t* trhs = (RHS(jt, fb));
//                 plasma_complex64_t* ttpup = (TPUP(jt, fb));
//                 int ld = nsize - fb*nb;
// // #pragma omp task depend(in:trhs[ji*ld:ji*ld + ld]) \
// //                  depend(in:dup[fb*nb:(fb*nb+ld)]) \
// //                  depend(inout:ttpup[ji*ld:ji*ld + ld])
//             {
//             plasma_complex64_t dlj = 1.0/dl[jt*nb + ji];
//                 int ld = nsize - fb*nb;
//                 for(int ii = 0; ii < ld; ++ii)
//                 {
//                     if ( creal( dr[fb*nb + ii] ) <= 1.0)
//                         (TPUP(jt, fb))[ji*ld + ii] = dlj * (RHS(jt, fb))[ji*ld + ii] + dup[fb*nb + ii] * (TPUP(jt, fb))[ji*ld + ii];
//                     else
//                         (TPUP(jt, fb))[ji*ld + ii] = (RHS(jt, fb))[ji*ld + ii]/dr[fb*nb + ii]/dl[jt*nb + ji] + (TPUP(jt, fb))[ji*ld + ii];
//                 }
//             }
//         }
//     }
// }
// //Loop remainder
// {
//     for(int ji = 0; ji < nsize - fb*nb; ++ji)
//     {
//         if(creal(dl[fb*nb + ji]) <= 1.0 )
//         {
//             for(int it = 0; it < fb; ++it)
//             {
//                 plasma_complex64_t* trhs = (RHS(fb, it));
//                 plasma_complex64_t* ttpup = (TPUP(fb, it));
// // #pragma omp task depend(in:trhs[ji*nb:ji*nb+nb]) \
// //                  depend(in:dup[it*nb:(it*nb+nb)]) \
// //                  depend(inout:ttpup[ji*nb:ji*nb+nb])
//                  {
//             plasma_complex64_t dlj = dl[fb*nb + ji];
//                 for(int ii = 0; ii < nb; ++ii)
//                 {
//                     if ( creal( dr[it*nb + ii] ) <= 1.0)
//                         (TPUP(fb, it))[ji*nb + ii] = (RHS(fb, it))[ji*nb + ii] + dr[it*nb + ii] * dl[fb*nb + ji] * (TPUP(fb, it))[ji*nb + ii];
//                     else
//                         (TPUP(fb, it))[ji*nb + ii] = dup[it*nb + ii] * (RHS(fb, it))[ji*nb + ii] + dlj * (TPUP(fb, it))[ji*nb + ii];
//                         
//                 }
//                  }
//             }
//             //REMAINDER
//              plasma_complex64_t* trhs = (RHS(fb, fb));
//                 plasma_complex64_t* ttpup = (TPUP(fb, fb));
//                 int ld = nsize - fb*nb;
// // #pragma omp task depend(in:trhs[ji*ld:ji*ld + ld]) \
// //                  depend(in:dup[fb*nb:(fb*nb+ld)]) \
// //                  depend(inout:ttpup[ji*ld:ji*ld + ld])
//             {
//             plasma_complex64_t dlj = dl[fb*nb + ji];
//                 for(int ii = 0; ii < nsize - fb*nb; ++ii)
//                 {
//                     if ( creal( dr[fb*nb + ii] ) <= 1.0)
//                         (TPUP(fb, fb))[ji*ld + ii] = (RHS(fb, fb))[ji*ld + ii] + dr[fb*nb + ii] * dl[fb*nb + ji] * (TPUP(fb, fb))[ji*ld + ii];
//                     else
//                         (TPUP(fb, fb))[ji*ld + ii] = dup[fb*nb + ii] * (RHS(fb, fb))[ji*ld + ii] + dlj * (TPUP(fb, fb))[ji*ld + ii];
//                         
//                 }
//             }
//         }
//         else
//         {
//             for(int it = 0; it < fb; ++it)
//             {
//                 plasma_complex64_t* trhs = (RHS(fb, it));
//                 plasma_complex64_t* ttpup = (TPUP(fb, it));
// // #pragma omp task depend(in:trhs[ji*nb:ji*nb+nb]) \
// //                  depend(in:dup[it*nb:(it*nb+nb)]) \
// //                  depend(inout:ttpup[ji*nb:ji*nb+nb])
//                  {
//             plasma_complex64_t dlj = 1.0/dl[fb*nb + ji];
//                 for(int ii = 0; ii < tpup->nb; ++ii)
//                 {
//                     if ( creal( dr[it*nb + ii] ) <= 1.0)
//                         (TPUP(fb, it))[ji*nb + ii] = dlj * (RHS(fb, it))[ji*nb + ii] + dup[it*nb + ii] * (TPUP(fb, it))[ji*nb + ii];
//                     else
//                         (TPUP(fb, it))[ji*nb + ii] = (RHS(fb, it))[ji*nb + ii]/dr[it*nb + ii]/dl[fb*nb + ji] + (TPUP(fb, it))[ji*nb + ii];
//                 }
//                  }
//             }
//             //REMAINDER
//                 plasma_complex64_t* trhs = (RHS(fb, fb));
//                 plasma_complex64_t* ttpup = (TPUP(fb, fb));
//                 int ld = nsize - fb*nb;
// // #pragma omp task depend(in:trhs[ji*ld:ji*ld + ld]) \
// //                  depend(in:dup[fb*nb:(fb*nb+ld)]) \
// //                  depend(inout:ttpup[ji*ld:ji*ld + ld])
//             {
//                 plasma_complex64_t dlj = 1.0/dl[fb*nb + ji];
//                 for(int ii = 0; ii < ld; ++ii)
//                 {
//                     if ( creal( dr[fb*nb + ii] ) <= 1.0)
//                         (TPUP(fb, fb))[ji*ld + ii] = dlj * (RHS(fb, fb))[ji*ld + ii] + dup[fb*nb + ii] * (TPUP(fb, fb))[ji*ld + ii];
//                     else
//                         (TPUP(fb, fb))[ji*ld + ii] = (RHS(fb, fb))[ji*ld + ii]/dr[fb*nb + ii]/dl[fb*nb + ji] + (TPUP(fb, fb))[ji*ld + ii];
//                 }
//             }
//         }
//     }
// }
