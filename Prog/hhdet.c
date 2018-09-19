#include <plasma.h>
#include <plasma_descriptor.h>
#include <plasma_types.h>

#include <stdio.h>

#define T(m, n) (plasma_complex64_t*)plasma_tile_addr(*T, m, n)

/** Calculate the determinant of a householder object
 */
void hhdet(plasma_desc_t* T, plasma_complex64_t* phase, int nvar)
{
    /* PLASMA matrices have a tile based layout that consists of gmt * gnt tiles.
     * We ignore a possible submatrix structure that plasma supports.
     * On top of that the QR algorithms that is applied to each tile
     * has inner blocking that has to be considered.
     */
    for (int i = 0; i < T->gmt; ++i)//loop over the tile rows
    {
        for(int j = 0; j < T->gnt; ++j)//loop over the tile columns
        {
            printf("Tile coordinates%i %i\n", i, j);
/*            
            for (int x = 0; x < T->nb; ++x)
            {
                for(int y =  0; y < T->mb; ++y)
                {
            plasma_complex64_t tau = ( T(i,j) )[x * T->mb + y];
            printf("%f ", creal(tau));
                }
                printf("\n");
            }*/
            for (int b = 0; b < (T->nb / T->mb); ++b )//loop over the blocks used by the QR decomposition
            {//Obviously we also assume that a tile has dimensions that are multiples of the inner block size...
                for (int k = 0; k < T->mb; ++k)
                {
                    plasma_complex64_t tau = ( T(i,j) )[b*T->mb*T->mb + k * T->mb + k];
                    printf("(C) %f%+fi\n", creal(tau  ), cimag(tau) );
                    if (creal(tau) != 0.0 || cimag(tau) != 0.0)
                    {
                        if(nvar == 1) tau = conj(tau);
                        double x = cabs(tau);
                        plasma_complex64_t z = 1.0 - 2.0*(tau/x) * (creal(tau)/x);
                        (*phase) *= z/ cabs(z) ;
                    }
//                 plasma_complex64_t temp = ( T(i,j) )[k * T->mb + k];
//                 printf("(C) %f%+fi\n", crealf(temp  ), cimagf(temp) );
                }
            }
        }
    }
}
