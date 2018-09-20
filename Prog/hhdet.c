#include <plasma.h>
#include <plasma_descriptor.h>
#include <plasma_types.h>

#include <stdio.h>

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
