#include <plasma.h>
#include <plasma_descriptor.h>
#include <plasma_types.h>

#include <stdio.h>

#define T(m, n) (plasma_complex64_t*)plasma_tile_addr(*T, m, n)

/** Calculate the determinant of a householder object
 */
void hhdet(plasma_desc_t* T, plasma_complex64_t* phase, int nvar)
{
    /*We have a tile based layout that consists of gmt * gnt tiles.
     * We ignore a possible submatrix structure that plasma supports.
     */
    for (int i = 0; i < T->gmt; ++i)
    {
        for(int j = 0; j < T->gnt; ++j)
        {
            printf("%i %i", i, j);
            for (int k = 0; k < T->mb; ++k)
            {
                plasma_complex64_t tau = ( T(i,j) )[k * T->mb + k];
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
