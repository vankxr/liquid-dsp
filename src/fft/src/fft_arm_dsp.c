/*
 * Copyright (c) 2007 - 2020 Joseph Gaeddert
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
// fft_arm.c : definitions for ARM DSP library FFTs
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#if __ARM_FEATURE_DSP == 1
#include "arm_math.h"
#include "arm_const_structs.h"
#endif
#include "liquid.internal.h"

// create FFT plan for ARM DSP library FFT
//  _nfft   :   FFT size
//  _x      :   input array [size: _nfft x 1]
//  _y      :   output array [size: _nfft x 1]
//  _dir    :   fft direction: {LIQUID_FFT_FORWARD, LIQUID_FFT_BACKWARD}
//  _method :   fft method
FFT(plan) FFT(_create_plan_arm_dsp)(unsigned int _nfft,
                                TC *         _x,
                                TC *         _y,
                                int          _dir,
                                int          _flags)
{
#if __ARM_FEATURE_DSP == 1
    // allocate plan and initialize all internal arrays to NULL
    FFT(plan) q = (FFT(plan)) malloc(sizeof(struct FFT(plan_s)));

    q->nfft      = _nfft;
    q->x         = _x;
    q->y         = _y;
    q->flags     = _flags;
    q->type      = (_dir == LIQUID_FFT_FORWARD) ? LIQUID_FFT_FORWARD : LIQUID_FFT_BACKWARD;
    q->direction = (_dir == LIQUID_FFT_FORWARD) ? LIQUID_FFT_FORWARD : LIQUID_FFT_BACKWARD;
    q->method    = LIQUID_FFT_METHOD_ARM_DSP;

    q->execute   = FFT(_execute_arm_dsp);

    if      (q->nfft == 16)   q->data.arm.instance = &arm_cfft_sR_f32_len16;
    else if (q->nfft == 32)   q->data.arm.instance = &arm_cfft_sR_f32_len32;
    else if (q->nfft == 64)   q->data.arm.instance = &arm_cfft_sR_f32_len64;
    else if (q->nfft == 128)  q->data.arm.instance = &arm_cfft_sR_f32_len128;
    else if (q->nfft == 256)  q->data.arm.instance = &arm_cfft_sR_f32_len256;
    else if (q->nfft == 512)  q->data.arm.instance = &arm_cfft_sR_f32_len512;
    else if (q->nfft == 1024) q->data.arm.instance = &arm_cfft_sR_f32_len1024;
    else if (q->nfft == 2048) q->data.arm.instance = &arm_cfft_sR_f32_len2048;
    else if (q->nfft == 4096) q->data.arm.instance = &arm_cfft_sR_f32_len4096;

    return q;
#else
    return LIQUID_EUMODE;
#endif
}

// destroy FFT plan
int FFT(_destroy_plan_arm_dsp)(FFT(plan) _q)
{
#if __ARM_FEATURE_DSP == 1
    // free main object memory
    free(_q);
    return LIQUID_OK;
#else
    return LIQUID_EUMODE;
#endif
}

// execute FFT
int FFT(_execute_arm_dsp)(FFT(plan) _q)
{
#if __ARM_FEATURE_DSP == 1
    memcpy(_q->y, _q->x, _q->nfft * sizeof(TC));

    // Executes in place
    arm_cfft_f32(_q->data.arm.instance, (float *)_q->y, _q->direction == LIQUID_FFT_FORWARD ? 0 : 1, 1);

    return LIQUID_OK;
#else
    return LIQUID_EUMODE;
#endif
}

