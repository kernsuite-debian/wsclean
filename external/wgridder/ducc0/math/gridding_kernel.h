/*
 *  This file is part of the MR utility library.
 *
 *  This code is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This code is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this code; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

/* Copyright (C) 2020-2021 Max-Planck-Society
   Author: Martin Reinecke */

#ifndef DUCC0_GRIDDING_KERNEL_H
#define DUCC0_GRIDDING_KERNEL_H

#include <algorithm>
#include <array>
#include <cstddef>
#include <functional>
#include <vector>
#include <memory>
#include <cmath>
#include <type_traits>
#include "ducc0/infra/useful_macros.h"
#include "ducc0/infra/error_handling.h"
#include "ducc0/infra/threading.h"
#include "ducc0/math/gl_integrator.h"
#include "ducc0/math/constants.h"

namespace ducc0 {

namespace detail_gridding_kernel {

using namespace std;

vector<double> getCoeffs(size_t W, size_t D, const function<double(double)> &func)
  {
  vector<double> coeff(W*(D+1));
  vector<double> chebroot(D+1);
  for (size_t i=0; i<=D; ++i)
    chebroot[i] = cos((2*i+1.)*pi/(2*D+2));
  vector<double> y(D+1), lcf(D+1), C((D+1)*(D+1)), lcf2(D+1);
  for (size_t i=0; i<W; ++i)
    {
    double l = -1+2.*i/double(W);
    double r = -1+2.*(i+1)/double(W);
    // function values at Chebyshev nodes
    double avg = 0;
    for (size_t j=0; j<=D; ++j)
      {
      y[j] = func(chebroot[j]*(r-l)*0.5 + (r+l)*0.5);
      avg += y[j];
      }
    avg/=(D+1);
    for (size_t j=0; j<=D; ++j)
      y[j] -= avg;
    // Chebyshev coefficients
    for (size_t j=0; j<=D; ++j)
      {
      lcf[j] = 0;
      for (size_t k=0; k<=D; ++k)
        lcf[j] += 2./(D+1)*y[k]*cos(j*(2*k+1)*pi/(2*D+2));
      }
    lcf[0] *= 0.5;
    // Polynomial coefficients
    fill(C.begin(), C.end(), 0.);
    C[0] = 1.;
    C[1*(D+1) + 1] = 1.;
    for (size_t j=2; j<=D; ++j)
      {
      C[j*(D+1) + 0] = -C[(j-2)*(D+1) + 0];
      for (size_t k=1; k<=j; ++k)
        C[j*(D+1) + k] = 2*C[(j-1)*(D+1) + k-1] - C[(j-2)*(D+1) + k];
      }
    for (size_t j=0; j<=D; ++j) lcf2[j] = 0;
    for (size_t j=0; j<=D; ++j)
      for (size_t k=0; k<=D; ++k)
        lcf2[k] += C[j*(D+1) + k]*lcf[j];
    lcf2[0] += avg;
    for (size_t j=0; j<=D; ++j)
      coeff[j*W + i] = lcf2[D-j];
    }
  return coeff;
  }


/*! A GriddingKernel is considered to be a symmetric real-valued function
    defined on the interval [-1; 1].
    This range is subdivided into W equal-sized parts. */
class GriddingKernel
  {
  public:
    virtual ~GriddingKernel() {}

    virtual size_t support() const = 0;

    /* Computes the correction function at a given coordinate.
       Useful coordinates lie in the range [0; 0.5]. */
    virtual double corfunc(double x) const = 0;

    /* Computes the correction function values at a coordinates
       [0, dx, 2*dx, ..., (n-1)*dx]  */
    virtual vector<double> corfunc(size_t n, double dx, int nthreads=1) const = 0;
  };

class KernelCorrection
  {
  protected:
    vector<double> x, wgtpsi;
    size_t supp;

  public:
    /* Compute correction factors for gridding kernel
       This implementation follows eqs. (3.8) to (3.10) of Barnett et al. 2018 */
    double corfunc(double v) const
      {
      double tmp=0;
      for (size_t i=0; i<x.size(); ++i)
        tmp += wgtpsi[i]*cos(pi*supp*v*x[i]);
      return 1./tmp;
      }
    /* Compute correction factors for gridding kernel
       This implementation follows eqs. (3.8) to (3.10) of Barnett et al. 2018 */
    vector<double> corfunc(size_t n, double dx, int nthreads=1) const
      {
      vector<double> res(n);
      execStatic(n, nthreads, 0, [&](auto &sched)
        {
        while (auto rng=sched.getNext()) for(auto i=rng.lo; i<rng.hi; ++i)
          res[i] = corfunc(i*dx);
        });
      return res;
      }
  };

class GLFullCorrection: public KernelCorrection
  {
  public:
    GLFullCorrection(size_t W, const function<double(double)> &func)
      {
      supp = W;
      size_t p = size_t(1.5*W)+2;
      GL_Integrator integ(2*p);
      x = integ.coordsSymmetric();
      wgtpsi = integ.weightsSymmetric();
      for (size_t i=0; i<x.size(); ++i)
        wgtpsi[i] *= func(x[i])*supp*0.5;
      }
  };

class HornerKernel: public GriddingKernel
  {
  private:
    size_t W, D;
    vector<double> coeff;
    KernelCorrection corr;

  public:
    HornerKernel(size_t W_, size_t D_, const function<double(double)> &func,
      const KernelCorrection &corr_)
      : W(W_), D(D_),
        coeff(getCoeffs(W_, D_, func)),
        corr(corr_)
      {}

    virtual size_t support() const { return W; }

    virtual double corfunc(double x) const { return corr.corfunc(x); }

    /* Computes the correction function values at a coordinates
       [0, dx, 2*dx, ..., (n-1)*dx]  */
    virtual vector<double> corfunc(size_t n, double dx, int nthreads=1) const
      { return corr.corfunc(n, dx, nthreads); }

    const vector<double> &Coeff() const { return coeff; }
    size_t degree() const { return D; }
  };

template<size_t W, typename Tsimd> class TemplateKernel
  {
  private:
    static constexpr auto D=W+3;
    using T = typename Tsimd::value_type;
    static constexpr auto vlen = Tsimd::size();
    static constexpr auto nvec = (W+vlen-1)/vlen;

    std::array<Tsimd,(D+1)*nvec> coeff;
    const T *scoeff;
    static constexpr auto sstride = nvec*vlen;
 
    void transferCoeffs(const vector<double> &input)
      {
      for (size_t j=0; j<=D; ++j)
        {
        for (size_t i=0; i<W; ++i)
          coeff[j*nvec + i/vlen][i%vlen] = T(input[j*W+i]);
        for (size_t i=W; i<vlen*nvec; ++i)
          coeff[j*nvec + i/vlen][i%vlen] = T(0);
        }
      }

  public:
    TemplateKernel(const HornerKernel &krn)
      : scoeff(reinterpret_cast<T *>(&coeff[0]))
      {
      MR_assert(W==krn.support(), "support mismatch");
      MR_assert(D==krn.degree(), "degree mismatch");
      transferCoeffs(krn.Coeff());
      }

    constexpr size_t support() const { return W; }

    [[gnu::always_inline]] void eval2s(T x, T y, T z, size_t nth, Tsimd * DUCC0_RESTRICT res) const
      {
      z = (z-nth)*2+(W-1);
      if constexpr (nvec==1)
        {
        auto tvalx = coeff[0];
        auto tvaly = coeff[0];
        auto tvalz = coeff[0];
        for (size_t j=1; j<=D; ++j)
          {
          tvalx = tvalx*x + coeff[j];
          tvaly = tvaly*y + coeff[j];
          tvalz = tvalz*z + coeff[j];
          }
        res[0] = tvalx*T(tvalz[nth]);
        res[1] = tvaly;
        }
      else
        {
        auto ptrz = scoeff+nth;
        auto tvalz = *ptrz;
        for (size_t j=1; j<=D; ++j)
          tvalz = tvalz*z + ptrz[j*sstride];
        for (size_t i=0; i<nvec; ++i)
          {
          auto tvalx = coeff[i];
          auto tvaly = coeff[i];
          for (size_t j=1; j<=D; ++j)
            {
            tvalx = tvalx*x + coeff[j*nvec+i];
            tvaly = tvaly*y + coeff[j*nvec+i];
            }
          res[i] = tvalx*tvalz;
          res[i+nvec] = tvaly;
          }
        }
      }
    [[gnu::always_inline]] void eval2(T x, T y, Tsimd * DUCC0_RESTRICT res) const
      {
      if constexpr (nvec==1)
        {
        auto tvalx = coeff[0];
        auto tvaly = coeff[0];
        for (size_t j=1; j<=D; ++j)
          {
          tvalx = tvalx*x + coeff[j];
          tvaly = tvaly*y + coeff[j];
          }
        res[0] = tvalx;
        res[1] = tvaly;
        }
      else
        {
        for (size_t i=0; i<nvec; ++i)
          {
          auto tvalx = coeff[i];
          auto tvaly = coeff[i];
          for (size_t j=1; j<=D; ++j)
            {
            tvalx = tvalx*x + coeff[j*nvec+i];
            tvaly = tvaly*y + coeff[j*nvec+i];
            }
          res[i] = tvalx;
          res[i+nvec] = tvaly;
          }
        }
      }
    [[gnu::always_inline]] void eval3(T x, T y, T z, Tsimd * DUCC0_RESTRICT res) const
      {
      if constexpr (nvec==1)
        {
        auto tvalx = coeff[0];
        auto tvaly = coeff[0];
        auto tvalz = coeff[0];
        for (size_t j=1; j<=D; ++j)
          {
          tvalx = tvalx*x + coeff[j];
          tvaly = tvaly*y + coeff[j];
          tvalz = tvalz*z + coeff[j];
          }
        res[0] = tvalx;
        res[1] = tvaly;
        res[2] = tvalz;
        }
      else
        {
        for (size_t i=0; i<nvec; ++i)
          {
          auto tvalx = coeff[i];
          auto tvaly = coeff[i];
          auto tvalz = coeff[i];
          for (size_t j=1; j<=D; ++j)
            {
            tvalx = tvalx*x + coeff[j*nvec+i];
            tvaly = tvaly*y + coeff[j*nvec+i];
            tvalz = tvalz*z + coeff[j*nvec+i];
            }
          res[i] = tvalx;
          res[i+nvec] = tvaly;
          res[i+2*nvec] = tvalz;
          }
        }
      }
  };

struct KernelParams
  {
  size_t W;
  double ofactor, epsilon, beta, e0, corr_range;
  };

const vector<KernelParams> KernelDB {
{ 4, 1.15,   0.025654879, 1.3873426689, 0.5436851297, 10.3501},
{ 4, 1.20,   0.013809249, 1.3008419165, 0.5902137484, 8.02317},
{ 4, 1.25,  0.0085840685, 1.3274088935, 0.5953499486, 6.25302},
{ 4, 1.30,  0.0057322498, 1.3617063353, 0.5965631622, 5.0918},
{ 4, 1.35,  0.0042494419, 1.3845499880, 0.5990241291, 4.32675},
{ 4, 1.40,  0.0033459552, 1.4405325088, 0.5924776015, 3.71334},
{ 4, 1.45,  0.0028187359, 1.4635220066, 0.5929442711, 3.30688},
{ 4, 1.50,  0.0023843943, 1.5539689162, 0.5772217314, 2.9249},
{ 4, 1.55,  0.0020343796, 1.5991008653, 0.5721765215, 2.66733},
{ 4, 1.60,  0.0017143851, 1.6581546365, 0.5644747137, 2.45114},
{ 4, 1.65,  0.0014730848, 1.7135331415, 0.5572788589, 2.27875},
{ 4, 1.70,  0.0012554492, 1.7464330378, 0.5548742415, 2.14472},
{ 4, 1.75,  0.0010610904, 1.7887326906, 0.5509877716, 2.02741},
{ 4, 1.80, 0.00090885567, 1.8122309426, 0.5502273972, 1.93387},
{ 4, 1.85,  0.0007757401, 1.8304451327, 0.5503967160, 1.85422},
{ 4, 1.90,  0.0006740398, 1.8484487383, 0.5502376937, 1.78503},
{ 4, 1.95, 0.00058655391, 1.8742215688, 0.5489738941, 1.72203},
{ 4, 2.00, 0.00051911189, 1.9069436300, 0.5468009434, 1.66473},
{ 4, 2.05, 0.00047127936, 1.9287029587, 0.5459425678, 1.61648},
{ 4, 2.10, 0.00042991098, 1.9468344976, 0.5456431243, 1.57393},
{ 4, 2.15, 0.00039952939, 1.9598362025, 0.5457007164, 1.53683},
{ 4, 2.20, 0.00036728958, 1.9734445042, 0.5460368586, 1.50285},
{ 4, 2.25, 0.00034355459, 1.9833876672, 0.5464865366, 1.47281},
{ 4, 2.30, 0.00032238422, 1.9922532404, 0.5470505385, 1.44556},
{ 4, 2.35, 0.00030354772, 2.0001857065, 0.5476974129, 1.42076},
{ 4, 2.40, 0.00029003195, 2.0059275365, 0.5482583426, 1.39847},
{ 4, 2.45, 0.00027493243, 2.0124393824, 0.5489923951, 1.37768},
{ 4, 2.50, 0.00026418063, 2.0171876964, 0.5495887377, 1.35887},
{ 5, 1.15,  0.0088036926, 1.4211620799, 0.5484370222, 20.4779},
{ 5, 1.20,  0.0045432118, 1.4604589193, 0.5502520137, 13.6414},
{ 5, 1.25,  0.0025659469, 1.5114537479, 0.5482371505, 9.84508},
{ 5, 1.30,  0.0014949902, 1.5662004850, 0.5453959646, 7.53251},
{ 5, 1.35, 0.00092874124, 1.5940645314, 0.5464869375, 6.14913},
{ 5, 1.40,  0.0005820084, 1.6193311874, 0.5477484983, 5.18009},
{ 5, 1.45, 0.00041837131, 1.6702721179, 0.5446584432, 4.4105},
{ 5, 1.50, 0.00032139657, 1.7106607912, 0.5430584562, 3.85869},
{ 5, 1.55, 0.00025831183, 1.7411526213, 0.5424476190, 3.44991},
{ 5, 1.60, 0.00021156623, 1.7694517444, 0.5419943230, 3.12801},
{ 5, 1.65, 0.00018112326, 1.8069777863, 0.5396536287, 2.86115},
{ 5, 1.70, 0.00015177086, 1.8378820613, 0.5382171164, 2.64774},
{ 5, 1.75, 0.00012345178, 1.8819830388, 0.5352849778, 2.46038},
{ 5, 1.80, 0.00010093043, 1.9225188886, 0.5327117935, 2.30621},
{ 5, 1.85, 8.5743423e-05, 1.9766846627, 0.5286610959, 2.16884},
{ 5, 1.90, 7.5167678e-05, 2.0116590189, 0.5267194291, 2.06069},
{ 5, 1.95, 6.5915521e-05, 2.0401734131, 0.5256331063, 1.9692},
{ 5, 2.00, 5.7747201e-05, 2.0640495669, 0.5251476146, 1.89055},
{ 5, 2.05, 5.1546264e-05, 2.0815539113, 0.5250991831, 1.82348},
{ 5, 2.10, 4.6026099e-05, 2.0967549727, 0.5253193374, 1.76438},
{ 5, 2.15, 4.1922392e-05, 2.1078775421, 0.5256707657, 1.71291},
{ 5, 2.20,  3.755449e-05, 2.1195059411, 0.5262644187, 1.66605},
{ 5, 2.25, 3.4356546e-05, 2.1278007625, 0.5268876862, 1.62478},
{ 5, 2.30, 3.1539667e-05, 2.1348794082, 0.5276099184, 1.58754},
{ 5, 2.35, 2.9093026e-05, 2.1406731238, 0.5284272705, 1.55386},
{ 5, 2.40, 2.7399686e-05, 2.1442622022, 0.5291533269, 1.5237},
{ 5, 2.45, 2.5638396e-05, 2.1486533599, 0.5300021540, 1.49569},
{ 5, 2.50, 2.4438353e-05, 2.1554191248, 0.5302677770, 1.46983},
{ 6, 1.15,  0.0018919684, 1.4284593523, 0.5456388809, 44.0159},
{ 6, 1.20, 0.00087379161, 1.4871949080, 0.5434184013, 25.5241},
{ 6, 1.25, 0.00052387586, 1.5596009923, 0.5384733141, 16.4885},
{ 6, 1.30, 0.00030833805, 1.6293990176, 0.5339545697, 11.6668},
{ 6, 1.35, 0.00018595126, 1.6818294794, 0.5315541173, 8.92588},
{ 6, 1.40, 0.00010913759, 1.7181576557, 0.5313982413, 7.20067},
{ 6, 1.45,  7.446073e-05, 1.7747725210, 0.5283789353, 5.91462},
{ 6, 1.50, 5.3826324e-05, 1.8169800206, 0.5271499621, 5.03369},
{ 6, 1.55, 4.0746477e-05, 1.8477830463, 0.5269687744, 4.40025},
{ 6, 1.60, 3.1441179e-05, 1.8721012077, 0.5274488807, 3.91897},
{ 6, 1.65, 2.6100718e-05, 1.8868745350, 0.5282375527, 3.55423},
{ 6, 1.70, 2.1528068e-05, 1.9342762857, 0.5255216043, 3.21515},
{ 6, 1.75, 1.7177115e-05, 1.9721895688, 0.5238810155, 2.94959},
{ 6, 1.80, 1.3650115e-05, 1.9987007558, 0.5232375255, 2.74023},
{ 6, 1.85, 1.0598995e-05, 2.0219705218, 0.5229660080, 2.56562},
{ 6, 1.90, 8.8157904e-06, 2.0671787180, 0.5203913566, 2.40316},
{ 6, 1.95, 7.6286922e-06, 2.1003673879, 0.5190199946, 2.27148},
{ 6, 2.00, 6.5649967e-06, 2.1272513974, 0.5182973590, 2.16058},
{ 6, 2.05, 5.7476558e-06, 2.1458295051, 0.5181563583, 2.06773},
{ 6, 2.10, 5.0756513e-06, 2.1511783437, 0.5195284847, 1.99018},
{ 6, 2.15, 4.4661935e-06, 2.1743106617, 0.5184286324, 1.91645},
{ 6, 2.20, 3.8877561e-06, 2.1867440456, 0.5188571997, 1.8531},
{ 6, 2.25, 3.4672484e-06, 2.1957177257, 0.5193179074, 1.79767},
{ 6, 2.30, 3.1012426e-06, 2.2033754263, 0.5198649665, 1.74797},
{ 6, 2.35, 2.7894219e-06, 2.2096913080, 0.5204820251, 1.70327},
{ 6, 2.40, 2.5794626e-06, 2.2135312754, 0.5210295408, 1.6635},
{ 6, 2.45, 2.3571404e-06, 2.2317879841, 0.5204756852, 1.62378},
{ 6, 2.50, 2.1615297e-06, 2.2535598472, 0.5194274130, 1.58705},
{ 7, 1.15, 0.00078476028, 1.5248706519, 0.5288306317, 77.3008},
{ 7, 1.20, 0.00027127166, 1.5739348793, 0.5287992619, 42.3005},
{ 7, 1.25, 0.00012594628, 1.6245240723, 0.5279217770, 26.3653},
{ 7, 1.30, 7.0214545e-05, 1.6835745981, 0.5257484101, 17.7989},
{ 7, 1.35, 4.1972457e-05, 1.7343424414, 0.5239793844, 13.0169},
{ 7, 1.40,  2.378019e-05, 1.7845017738, 0.5224266045, 9.99361},
{ 7, 1.45, 1.3863408e-05, 1.8180597789, 0.5221834768, 8.09152},
{ 7, 1.50, 9.1605353e-06, 1.8680822720, 0.5206277502, 6.64755},
{ 7, 1.55,  6.479159e-06, 1.9188980015, 0.5183134674, 5.61353},
{ 7, 1.60, 4.6544571e-06, 1.9536166143, 0.5178695891, 4.87825},
{ 7, 1.65, 3.5489761e-06, 1.9786267068, 0.5178430252, 4.33279},
{ 7, 1.70, 2.7030348e-06, 2.0027666534, 0.5178577604, 3.89791},
{ 7, 1.75, 2.0533894e-06, 2.0289949199, 0.5176300336, 3.54136},
{ 7, 1.80, 1.6069122e-06, 2.0596412946, 0.5167551932, 3.24397},
{ 7, 1.85, 1.2936794e-06, 2.0720606842, 0.5178747891, 3.01314},
{ 7, 1.90, 1.0768664e-06, 2.0908981740, 0.5181009847, 2.81272},
{ 7, 1.95, 9.0890421e-07, 2.1086185697, 0.5184537843, 2.64203},
{ 7, 2.00, 7.7488775e-07, 2.1278284187, 0.5186377792, 2.49387},
{ 7, 2.05, 6.8025539e-07, 2.1300505355, 0.5201567726, 2.37603},
{ 7, 2.10, 6.0222531e-07, 2.1361214247, 0.5212397206, 2.27023},
{ 7, 2.15, 5.0130101e-07, 2.2231545475, 0.5137738214, 2.14207},
{ 7, 2.20, 4.2248762e-07, 2.2408449906, 0.5136504150, 2.05745},
{ 7, 2.25, 3.6494171e-07, 2.2556844458, 0.5135180892, 1.98365},
{ 7, 2.30, 3.1538194e-07, 2.2684881766, 0.5135795865, 1.91808},
{ 7, 2.35, 2.7282924e-07, 2.2815394648, 0.5136137508, 1.8589},
{ 7, 2.40, 2.4350524e-07, 2.2939223534, 0.5134777347, 1.80593},
{ 7, 2.45, 2.1263032e-07, 2.3041489588, 0.5138203874, 1.75785},
{ 7, 2.50, 1.9134836e-07, 2.3076482212, 0.5145035417, 1.716},
{ 8, 1.20, 7.8028732e-05, 1.6209261450, 0.5219287175, 72.156},
{ 8, 1.25, 2.7460918e-05, 1.6851585171, 0.5199250590, 41.0455},
{ 8, 1.30, 1.3421658e-05, 1.7442373315, 0.5182155619, 26.3042},
{ 8, 1.35, 7.5158217e-06, 1.7876782642, 0.5176319503, 18.6037},
{ 8, 1.40, 4.2472384e-06, 1.8294321912, 0.5171860211, 13.8907},
{ 8, 1.45, 2.5794802e-06, 1.8716918210, 0.5161733611, 10.8298},
{ 8, 1.50, 1.6131994e-06, 1.9213040541, 0.5145350888, 8.6747},
{ 8, 1.55, 1.0974814e-06, 1.9637229131, 0.5134005827, 7.19139},
{ 8, 1.60,  7.531955e-07, 2.0002761373, 0.5128849282, 6.11792},
{ 8, 1.65, 5.5097346e-07, 2.0275645736, 0.5127082324, 5.33605},
{ 8, 1.70, 4.0136726e-07, 2.0498410409, 0.5130237662, 4.73241},
{ 8, 1.75,  2.906467e-07, 2.0731585170, 0.5131757153, 4.2462},
{ 8, 1.80, 2.1834922e-07, 2.0907418726, 0.5136046561, 3.86143},
{ 8, 1.85, 1.6329905e-07, 2.1164552354, 0.5133333878, 3.53086},
{ 8, 1.90, 1.2828598e-07, 2.1261570160, 0.5143004427, 3.27482},
{ 8, 1.95, 1.0171134e-07, 2.1363206613, 0.5152354910, 3.05604},
{ 8, 2.00, 8.1881369e-08, 2.1397013368, 0.5166895497, 2.87454},
{ 8, 2.05, 6.9121193e-08, 2.1466071700, 0.5176145380, 2.71496},
{ 8, 2.10, 5.9525932e-08, 2.1510526592, 0.5186914118, 2.57749},
{ 8, 2.15, 5.2942463e-08, 2.2365737543, 0.5125850104, 2.40676},
{ 8, 2.20, 4.3612361e-08, 2.2635555483, 0.5116114910, 2.29267},
{ 8, 2.25, 3.6764793e-08, 2.2808513000, 0.5112823144, 2.19713},
{ 8, 2.30, 3.0899101e-08, 2.2961118291, 0.5111472899, 2.11277},
{ 8, 2.35, 2.5951523e-08, 2.3025419974, 0.5117804832, 2.04051},
{ 8, 2.40, 2.2598633e-08, 2.3146967576, 0.5117074796, 1.97372},
{ 8, 2.45, 1.9029665e-08, 2.3186279028, 0.5125332192, 1.91551},
{ 8, 2.50, 1.6752523e-08, 2.3321309669, 0.5124348616, 1.85985},
{ 9, 1.25, 8.0261034e-06, 1.7103888450, 0.5164129862, 66.625},
{ 9, 1.30, 3.2272675e-06, 1.7768638337, 0.5141821303, 39.7671},
{ 9, 1.35, 1.6398132e-06, 1.8259273732, 0.5131939428, 26.7133},
{ 9, 1.40, 8.5542435e-07, 1.8706775936, 0.5126332399, 19.1697},
{ 9, 1.45, 4.8998062e-07, 1.9079562176, 0.5122630978, 14.5697},
{ 9, 1.50, 2.8357238e-07, 1.9376265737, 0.5127460716, 11.5443},
{ 9, 1.55, 1.7632448e-07, 1.9831130904, 0.5114276594, 9.30856},
{ 9, 1.60, 1.1241387e-07, 2.0031047508, 0.5126932345, 7.83305},
{ 9, 1.65, 8.0252028e-08, 2.0285383278, 0.5127726006, 6.70894},
{ 9, 1.70,  5.741767e-08, 2.0574910347, 0.5124426828, 5.82994},
{ 9, 1.75, 4.0256578e-08, 2.0895174008, 0.5117693191, 5.13211},
{ 9, 1.80, 2.8882533e-08, 2.1256913951, 0.5104744301, 4.57008},
{ 9, 1.90, 1.5495159e-08, 2.1522089772, 0.5119205380, 3.81106},
{ 9, 1.95, 1.1746262e-08, 2.1630258913, 0.5127106930, 3.52484},
{ 9, 2.00, 9.1018891e-09, 2.1705676763, 0.5137317339, 3.28511},
{10, 1.30, 7.6037409e-07, 1.8037273215, 0.5111938042, 59.9053},
{10, 1.35, 3.5081762e-07, 1.8371602342, 0.5120315623, 39.3257},
{10, 1.40, 1.7691912e-07, 1.8539150977, 0.5144182834, 28.066},
{10, 1.45, 9.5898587e-08, 1.9083076984, 0.5123796414, 20.1661},
{10, 1.50, 5.1649488e-08, 1.9549482618, 0.5110787344, 15.2855},
{10, 1.55, 2.9344166e-08, 1.9935039498, 0.5103977246, 12.0965},
{10, 1.60, 1.6984065e-08, 2.0134235964, 0.5114751650, 9.98157},
{10, 1.65, 1.1201377e-08, 2.0278278839, 0.5124699585, 8.46133},
{10, 1.70, 7.7392472e-09, 2.0550888710, 0.5123896568, 7.2305},
{10, 1.75, 5.4226206e-09, 2.0984256136, 0.5108595057, 6.22071},
{10, 1.80, 3.8051062e-09, 2.1302650200, 0.5100239564, 5.47481},
{10, 1.85, 2.6039483e-09, 2.1585185874, 0.5095352903, 4.88572},
{10, 1.90, 1.8492238e-09, 2.1707445630, 0.5101429537, 4.44066},
{10, 1.95, 1.3147032e-09, 2.1817982009, 0.5108284721, 4.07098},
{10, 2.00, 9.6449676e-10, 2.2018582745, 0.5109329270, 3.74388},
{11, 1.30, 1.8775393e-07, 1.8205994670, 0.5094502380, 90.9541},
{11, 1.35, 7.7356019e-08, 1.8518214951, 0.5103702455, 57.475},
{11, 1.40, 3.2528373e-08, 1.8880512818, 0.5108048393, 38.5515},
{11, 1.45, 1.6597763e-08, 1.9165471477, 0.5113840356, 27.7079},
{11, 1.50, 8.8600675e-09, 1.9625071236, 0.5103156252, 20.383},
{11, 1.55, 4.8306858e-09, 2.0049198756, 0.5093020852, 15.6747},
{11, 1.60, 2.5766957e-09, 2.0295371707, 0.5098854450, 12.6302},
{11, 1.65,  1.559903e-09, 2.0444486400, 0.5107653615, 10.5256},
{11, 1.70, 9.9258899e-10, 2.0627527204, 0.5113810966, 8.91137},
{11, 1.75, 6.9033069e-10, 2.0881262613, 0.5114471968, 7.63199},
{11, 1.80, 4.9444651e-10, 2.1384158198, 0.5092846283, 6.54168},
{11, 1.85, 3.3291902e-10, 2.1714286439, 0.5084380285, 5.75283},
{11, 1.90, 2.3006797e-10, 2.1889103178, 0.5086174417, 5.16381},
{11, 1.95, 1.5880468e-10, 2.2013464966, 0.5091627247, 4.68976},
{11, 2.00, 1.1177184e-10, 2.2147511730, 0.5096421117, 4.29102},
{12, 1.35, 1.7038991e-08, 1.8619112597, 0.5093832337, 84.1509},
{12, 1.40, 6.5438748e-09, 1.9069147481, 0.5089479889, 53.6538},
{12, 1.45, 2.9874764e-09, 1.9318398074, 0.5098082325, 37.6304},
{12, 1.50, 1.4920459e-09, 1.9628483155, 0.5100985753, 27.3958},
{12, 1.55, 8.0989276e-10, 2.0129847811, 0.5085327805, 20.3394},
{12, 1.60, 4.1660575e-10, 2.0517921747, 0.5079102398, 15.8407},
{12, 1.65, 2.3539727e-10, 2.0698388400, 0.5085131064, 12.9551},
{12, 1.70, 1.3497289e-10, 2.0887365361, 0.5090417146, 10.8069},
{12, 1.75, 8.3256938e-11, 2.1069557330, 0.5095920671, 9.18319},
{12, 1.80, 5.8834619e-11, 2.1359415217, 0.5091887069, 7.87128},
{12, 1.90, 2.6412908e-11, 2.2006369514, 0.5075889699, 6.01663},
{12, 1.95, 1.7189689e-11, 2.2146741638, 0.5080017404, 5.41103},
{12, 2.00, 1.2174796e-11, 2.2431392199, 0.5075191177, 4.87223},
{13, 1.40, 1.3183319e-09, 1.9191093732, 0.5077468429, 75.1325},
{13, 1.45, 5.3619083e-10, 1.9470729711, 0.5082824666, 50.9341},
{13, 1.50, 2.2322801e-10, 1.9739868782, 0.5088959158, 36.3221},
{13, 1.55, 1.1526271e-10, 2.0033259275, 0.5090866602, 26.9396},
{13, 1.60, 6.1426379e-11, 2.0510989715, 0.5078563876, 20.273},
{13, 1.65, 3.4052024e-11, 2.0818527907, 0.5074775828, 16.0817},
{13, 1.70, 1.8584229e-11, 2.1045643925, 0.5077204961, 13.1667},
{13, 1.75, 1.0810235e-11, 2.1243271913, 0.5081633997, 11.0246},
{13, 1.80,  7.377578e-12, 2.1467756379, 0.5081996246, 9.37882},
{13, 1.85, 4.9615599e-12, 2.1851483818, 0.5071383353, 8.01062},
{13, 1.90, 3.2630681e-12, 2.2086376720, 0.5068763959, 7.01814},
{13, 1.95,   2.08439e-12, 2.2299240544, 0.5068023663, 6.22754},
{13, 2.00, 1.3665993e-12, 2.2491841265, 0.5068976036, 5.58775},
{14, 1.45, 1.0233715e-10, 1.9577546112, 0.5073185014, 69.119},
{14, 1.50, 4.0219827e-11, 1.9908513900, 0.5074206790, 47.5913},
{14, 1.55, 1.8311452e-11, 2.0199455000, 0.5076279818, 34.5365},
{14, 1.60, 9.3257786e-12, 2.0564687260, 0.5072732788, 25.7812},
{14, 1.65, 5.0810838e-12, 2.0910978691, 0.5066391006, 19.9888},
{14, 1.70, 2.6125118e-12, 2.1169092067, 0.5066457965, 16.0644},
{14, 1.75, 1.3838674e-12, 2.1429649758, 0.5066883335, 13.1916},
{14, 1.80, 8.4432728e-13, 2.1058949451, 0.5120095847, 11.5897},
{14, 1.85, 5.6112166e-13, 2.1848347409, 0.5069391840, 9.4964},
{14, 1.90, 3.6520481e-13, 2.2144100981, 0.5063398541, 8.19239},
{14, 1.95, 2.2835769e-13, 2.2387205994, 0.5061044697, 7.18565},
{14, 2.00, 1.5660252e-13, 2.2662033397, 0.5057027272, 6.36102},
{15, 1.50, 6.7738438e-12, 2.0031350908, 0.5063691270, 62.5687},
{15, 1.55,  2.664519e-12, 2.0276909655, 0.5068592217, 44.7014},
{15, 1.60, 1.2501824e-12, 2.0594275589, 0.5068953514, 32.856},
{15, 1.65, 7.1442756e-13, 2.0958989255, 0.5061994449, 24.9268},
{15, 1.70, 3.6768121e-13, 2.1256130851, 0.5059364487, 19.631},
{15, 1.75, 1.8650547e-13, 2.1497846837, 0.5060649594, 15.9255},
{15, 1.80, 1.1458006e-13, 2.1307035522, 0.5091470333, 13.7024},
{15, 1.85, 7.0782501e-14, 2.1464045447, 0.5101577708, 11.5736},
{15, 1.90,  4.407902e-14, 2.2222258168, 0.5056705159, 9.54393},
{15, 1.95, 2.6655793e-14, 2.2419149837, 0.5057005374, 8.31949},
{15, 2.00, 1.7108261e-14, 2.2390528894, 0.5072248580, 7.4432},
{16, 1.50, 1.2100308e-12, 2.0130787701, 0.5055587965, 82.3371},
{16, 1.55, 4.6082202e-13, 2.0438032614, 0.5056309683, 56.9825},
{16, 1.60, 1.7883238e-13, 2.0329561822, 0.5089045671, 43.5164},
{16, 1.65, 9.2853815e-14, 2.0494514743, 0.5103582604, 32.9027},
{16, 1.70, 5.6614567e-14, 2.0925119791, 0.5083767402, 25.0837},
{16, 1.75,  2.875391e-14, 2.1461524027, 0.5062037834, 19.4122},
{16, 1.80, 1.6578982e-14, 2.1490040175, 0.5082721830, 16.1309},
{16, 1.85, 1.1782751e-14, 2.1811826814, 0.5072570059, 13.3536},
{16, 1.90, 8.9196865e-15, 2.1981176583, 0.5075840871, 11.3725},
{16, 1.95, 6.6530006e-15, 2.2340011350, 0.5060133105, 9.71037},
{16, 2.00, 5.0563492e-15, 2.2621631913, 0.5056924675, 8.42695}
};

template<typename T> T esknew (T v, T beta, T e0)
  {
  auto tmp = (1-v)*(1+v);
  auto tmp2 = tmp>=0;
  return tmp2*exp(beta*(pow(tmp*tmp2, e0)-1));
  }

auto selectKernel(size_t idx)
  {
  MR_assert(idx<KernelDB.size(), "no appropriate kernel found");
  auto supp = KernelDB[idx].W;
  auto beta = KernelDB[idx].beta*supp;
  auto e0 = KernelDB[idx].e0;
  auto lam = [beta,e0](double v){return esknew(v, beta, e0);};
  return make_shared<HornerKernel>(supp, supp+3, lam, GLFullCorrection(supp, lam));
  }

bool acceptable(size_t i)
  { return KernelDB[i].corr_range < 10.; }

template<typename T> constexpr inline size_t Wmax()
  { return is_same<T,float>::value ? 8 : 16; }

/*! Returns the best matching 2-parameter ES kernel for the given oversampling
    factor and error. */
template<typename T> auto selectKernel(double ofactor, double epsilon)
  {
  size_t Wmin = Wmax<T>();
  size_t idx = KernelDB.size();
  for (size_t i=0; i<KernelDB.size(); ++i)
    if ((KernelDB[i].ofactor<=ofactor) && (KernelDB[i].epsilon<=epsilon)
      && (KernelDB[i].W<=Wmin) && acceptable(i))
      {
      idx = i;
      Wmin = KernelDB[i].W;
      }
  return selectKernel(idx);
  }
template<typename T> auto selectKernel(double ofactor, double epsilon, size_t idx)
  {
  return (idx<KernelDB.size()) ?
    selectKernel(idx) : selectKernel<T>(ofactor, epsilon);
  }

template<typename T> auto getAvailableKernels(double epsilon,
  double ofactor_min=1.1, double ofactor_max=2.6)
  {
  vector<double> ofc(20, ofactor_max);
  vector<size_t> idx(20, KernelDB.size());
  size_t Wlim = Wmax<T>();
  for (size_t i=0; i<KernelDB.size(); ++i)
    {
    auto ofactor = KernelDB[i].ofactor;
    size_t W = KernelDB[i].W;
    if ((W<=Wlim) && (KernelDB[i].epsilon<=epsilon)
      && (ofactor<ofc[W]) && (ofactor>=ofactor_min) && acceptable(i))
      {
      ofc[W] = ofactor;
      idx[W] = i;
      }
    }
  vector<size_t> res;
  for (auto v: idx)
    if (v<KernelDB.size()) res.push_back(v);
  MR_assert(!res.empty(), "no appropriate kernel found");
  return res;
  }

}

using detail_gridding_kernel::GriddingKernel;
using detail_gridding_kernel::selectKernel;
using detail_gridding_kernel::getAvailableKernels;
using detail_gridding_kernel::HornerKernel;
using detail_gridding_kernel::TemplateKernel;
using detail_gridding_kernel::KernelParams;
using detail_gridding_kernel::KernelDB;

}

#endif