//
// Decimator
//

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "decim_internal.h"
#include "firdes.h"

// Print
void decim_print(decim _d)
{
    printf("decimator object:\n");
    printf("\tD\t: %u\n", _d->D);
    printf("\tfc\t: %5.2f\n", _d->fc);
    printf("\tb\t: %5.2f\n", _d->b);
    printf("\tt\t: %5.2f\n", _d->t);
    printf("\tslsl\t: %5.2f\n", _d->slsl);
    printf("\th\t: %u taps\n", _d->h_len);
}

// Debug print
void decim_debug_print(decim _d)
{
    decim_print(_d);
    unsigned int i;
    for (i=0; i<_d->h_len; i++)
        printf("  h(%u) = %E;\n", i+1, _d->h[i]);
}

// Execute
void xdecim_execute(decim _d, float * _x, unsigned int _x_len, float * _y, unsigned int _y_len)
{
    unsigned int i, n=0;
    for (i=0; i<(_x_len-_d->h_len); i+=2) {
        decim_dotprod(_d->h, _d->h_len, &_x[i], &_y[n]);
        n++;
    }
}

