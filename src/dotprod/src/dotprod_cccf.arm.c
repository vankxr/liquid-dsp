/*
 * Copyright (c) 2007 - 2021 Joseph Gaeddert
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

// 
// Floating-point dot product (ARM DSP)
//

#if __ARM_FEATURE_DSP == 1
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "liquid.internal.h"

#include "arm_math.h"

// basic dot product (ordinal calculation)
int dotprod_cccf_run(float complex * _h,
                     float complex * _x,
                     unsigned int    _n,
                     float complex * _y)
{
    float complex r = 0;
    unsigned int i;
    for (i=0; i<_n; i++)
        r += _h[i] * _x[i];
    *_y = r;
    return LIQUID_OK;
}

// basic dot product (ordinal calculation) with loop unrolled
int dotprod_cccf_run4(float complex * _h,
                      float complex * _x,
                      unsigned int    _n,
                      float complex * _y)
{
    float complex r = 0;

    // t = 4*(floor(_n/4))
    unsigned int t=(_n>>2)<<2; 

    // compute dotprod in groups of 4
    unsigned int i;
    for (i=0; i<t; i+=4) {
        r += _h[i]   * _x[i];
        r += _h[i+1] * _x[i+1];
        r += _h[i+2] * _x[i+2];
        r += _h[i+3] * _x[i+3];
    }

    // clean up remaining
    for ( ; i<_n; i++)
        r += _h[i] * _x[i];

    *_y = r;
    return LIQUID_OK;
}


//
// structured ARM DSP dot product
//

struct dotprod_cccf_s {
    float complex * h;             // coefficients array
    unsigned int n;     // length
};

dotprod_cccf dotprod_cccf_create(float complex * _h,
                                 unsigned int    _n)
{
    dotprod_cccf q = (dotprod_cccf)malloc(sizeof(struct dotprod_cccf_s));
    q->n = _n;

    // allocate memory for coefficients
    q->h = (float complex*) malloc((q->n)*sizeof(float complex));

    // move coefficients
    memmove(q->h, _h, (q->n)*sizeof(float complex));

    // return object
    return q;
}

dotprod_cccf dotprod_cccf_create_rev(float complex * _h,
                                     unsigned int    _n)
{
    dotprod_cccf q = (dotprod_cccf)malloc(sizeof(struct dotprod_cccf_s));
    q->n = _n;

    // allocate memory for coefficients
    q->h = (float complex*) malloc((q->n)*sizeof(float complex));

    // copy coefficients in time-reversed order
    unsigned int i;
    for (i=0; i<_n; i++)
        q->h[i] = _h[_n-i-1];

    // return object
    return q;
}

// re-create the structured dotprod object
dotprod_cccf dotprod_cccf_recreate(dotprod_cccf    _q,
                                   float complex * _h,
                                   unsigned int    _n)
{
    // check to see if length has changed
    if (_q->n != _n) {
        // set new length
        _q->n = _n;

        // re-allocate memory
        _q->h = (float complex*) realloc(_q->h, (_q->n)*sizeof(float complex));
    }

    // move new coefficients
    memmove(_q->h, _h, (_q->n)*sizeof(float complex));

    // return re-structured object
    return _q;
}

// re-create the structured dotprod object, coefficients reversed
dotprod_cccf dotprod_cccf_recreate_rev(dotprod_cccf    _q,
                                       float complex * _h,
                                       unsigned int    _n)
{
    // check to see if length has changed
    if (_q->n != _n) {
        // set new length
        _q->n = _n;

        // re-allocate memory
        _q->h = (float complex*) realloc(_q->h, (_q->n)*sizeof(float complex));
    }

    // copy coefficients in time-reversed order
    unsigned int i;
    for (i=0; i<_n; i++)
        _q->h[i] = _h[_n-i-1];

    // return re-structured object
    return _q;
}

int dotprod_cccf_destroy(dotprod_cccf _q)
{
    free(_q->h);    // free coefficients memory
    free(_q);       // free main object memory
    return LIQUID_OK;
}

int dotprod_cccf_print(dotprod_cccf _q)
{
    printf("dotprod [arm-dsp, %u coefficients]:\n", _q->n);
    unsigned int i;
    for (i=0; i<_q->n; i++) {
        printf("  %4u: %12.8f + j*%12.8f\n", i, crealf(_q->h[i]), cimagf(_q->h[i]));
    }
    return LIQUID_OK;
}

// execute structured dot product
//  _q      :   dotprod object
//  _x      :   input array
//  _y      :   output sample
int dotprod_cccf_execute(dotprod_cccf    _q,
                         float complex * _x,
                         float complex * _y)
{
    float re, im;
    arm_cmplx_dot_prod_f32((float *)_q->h, (float *)_x, _q->n, &re, &im);
    *_y = re + im * _Complex_I;
    return LIQUID_OK;
}

#else
#include "dotprod_cccf.c"
#endif
