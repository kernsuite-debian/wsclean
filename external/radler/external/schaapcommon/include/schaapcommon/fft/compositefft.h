// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef SCHAAPCOMMON_FFT_COMPOSITEFFT_H_
#define SCHAAPCOMMON_FFT_COMPOSITEFFT_H_

#include <fftw3.h>
#include <aocommon/staticfor.h>

namespace schaapcommon {
namespace fft {
void FftR2CComposite(fftwf_plan plan_r2c, fftwf_plan plan_c2c,
                     size_t image_height, size_t image_width, const float* in,
                     fftwf_complex* out, aocommon::StaticFor<size_t>& loop);

void FftC2RComposite(fftwf_plan plan_c2c, fftwf_plan plan_c2r,
                     size_t image_height, size_t image_width,
                     const fftwf_complex* in, float* out,
                     aocommon::StaticFor<size_t>& loop);
}  // namespace fft
}  // namespace schaapcommon

#endif