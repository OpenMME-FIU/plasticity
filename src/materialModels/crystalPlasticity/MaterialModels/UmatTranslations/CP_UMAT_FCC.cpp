#include <fem.hpp> // Fortran EMulation library of fable module
#include "../../../include/crystalPlasticity.h"

namespace cp_fcc {

using namespace fem::major_types;

using fem::common;

arr_ref<fem::real_star_8, 2>
matmult(
        arr_cref<fem::real_star_8, 2> a,
        arr_cref<fem::real_star_8, 2> b)
{
    fem::real_star_8 return_value = fem::zero<fem::real_star_8>();
    a(dimension(3, 3));
    b(dimension(3, 3));
    int i = fem::int0;
    int j = fem::int0;
    int k = fem::int0;
    FEM_DO_SAFE(i, 1, 3) {
        FEM_DO_SAFE(j, 1, 3) {
            return_value(i, j) = 0.e0;
            FEM_DO_SAFE(k, 1, 3) {
                return_value(i, j) += a(i, k) * b(k, j);
            }
        }
    }
    return return_value;
}

void
slip1(
  arr_ref<fem::real_star_8, 2> slpdir1,
  arr_ref<fem::real_star_8, 2> slpnor1,
  arr_cref<fem::real_star_8, 2> slpdir0,
  arr_cref<fem::real_star_8, 2> slpnor0,
  arr_cref<fem::real_star_8, 2> fel,
  arr_cref<fem::real_star_8, 2> felinv)
{
  slpdir1(dimension(3, 12));
  slpnor1(dimension(3, 12));
  slpdir0(dimension(3, 12));
  slpnor0(dimension(3, 12));
  fel(dimension(3, 3));
  felinv(dimension(3, 3));
  //C
  int n = fem::int0;
  int i = fem::int0;
  int k = fem::int0;
  FEM_DO_SAFE(n, 1, 12) {
    FEM_DO_SAFE(i, 1, 3) {
      slpdir1(i, n) = 0.0f;
      slpnor1(i, n) = 0.0f;
      FEM_DO_SAFE(k, 1, 3) {
        slpdir1(i, n) += fel(i, k) * slpdir0(k, n);
        slpnor1(i, n) += slpnor0(k, n) * felinv(k, i);
      }
    }
  }
  //C
}

void
green(
  arr_ref<fem::real_star_8, 2> eel,
  arr_cref<fem::real_star_8, 2> fel)
{
  eel(dimension(3, 3));
  fel(dimension(3, 3));
  //C
  int i = fem::int0;
  int k = fem::int0;
  arr_2d<3, 3, fem::real_star_8> felt(fem::fill0);
  FEM_DO_SAFE(i, 1, 3) {
    FEM_DO_SAFE(k, 1, 3) {
      felt(i, k) = fel(k, i);
    }
  }
  //C
  int n = fem::int0;
  FEM_DO_SAFE(i, 1, 3) {
    FEM_DO_SAFE(k, 1, 3) {
      eel(i, k) = 0;
      FEM_DO_SAFE(n, 1, 3) {
        eel(i, k) += felt(i, n) * fel(n, k);
      }
    }
  }
  //C
  eel(1, 1) = eel(1, 1) - 1.0f;
  eel(2, 2) = eel(2, 2) - 1.0f;
  eel(3, 3) = eel(3, 3) - 1.0f;
  FEM_DO_SAFE(i, 1, 3) {
    FEM_DO_SAFE(k, 1, 3) {
      eel(i, k) = eel(i, k) / 2.0f;
    }
  }
  //C
}

void
piola(
  arr_ref<fem::real_star_8, 2> spk2,
  arr_cref<fem::real_star_8, 4> c0,
  arr_cref<fem::real_star_8, 2> eel)
{
  spk2(dimension(3, 3));
  c0(dimension(3, 3, 3, 3));
  eel(dimension(3, 3));
  //C
  int i = fem::int0;
  int k = fem::int0;
  int m = fem::int0;
  int n = fem::int0;
  FEM_DO_SAFE(i, 1, 3) {
    FEM_DO_SAFE(k, 1, 3) {
      spk2(i, k) = 0.0f;
      FEM_DO_SAFE(m, 1, 3) {
        FEM_DO_SAFE(n, 1, 3) {
          spk2(i, k) += c0(i, k, m, n) * eel(m, n);
        }
      }
    }
  }
  //C
}

void
cauchy(
  arr_ref<fem::real_star_8, 2> s,
  arr_cref<fem::real_star_8, 2> fel,
  arr_cref<fem::real_star_8, 2> felinv,
  arr_cref<fem::real_star_8, 2> spk2)
{
  s(dimension(3, 3));
  fel(dimension(3, 3));
  felinv(dimension(3, 3));
  spk2(dimension(3, 3));
  //C
  fem::real_star_8 det = fel(1, 1) * (fel(2, 2) * fel(3, 3) - fel(2,
    3) * fel(3, 2)) - fel(1, 2) * (fel(2, 1) * fel(3, 3) - fel(2,
    3) * fel(3, 1)) + fel(1, 3) * (fel(2, 1) * fel(3, 2) - fel(2,
    2) * fel(3, 1));
  //C
  int i = fem::int0;
  int k = fem::int0;
  arr_2d<3, 3, fem::real_star_8> a(fem::fill0);
  int n = fem::int0;
  FEM_DO_SAFE(i, 1, 3) {
    FEM_DO_SAFE(k, 1, 3) {
      a(i, k) = 0;
      FEM_DO_SAFE(n, 1, 3) {
        a(i, k) += fel(i, n) * spk2(n, k);
      }
    }
  }
  FEM_DO_SAFE(i, 1, 3) {
    FEM_DO_SAFE(k, 1, 3) {
      s(i, k) = 0;
      FEM_DO_SAFE(n, 1, 3) {
        s(i, k) += a(i, n) * felinv(n, k);
      }
    }
  }
  //C
  if (det == 0) {
    FEM_DO_SAFE(i, 1, 3) {
      FEM_DO_SAFE(k, 1, 3) {
        s(i, k) = 0.0f;
      }
    }
  }
  if (det != 0) {
    FEM_DO_SAFE(i, 1, 3) {
      FEM_DO_SAFE(k, 1, 3) {
        s(i, k) = s(i, k) / det;
      }
    }
  }
  //C
}

void
rshrs(
  arr_ref<fem::real_star_8> tau,
  arr_cref<fem::real_star_8, 2> s,
  arr_cref<fem::real_star_8, 2> slpdir1,
  arr_cref<fem::real_star_8, 2> slpnor1)
{
  tau(dimension(12));
  s(dimension(3, 3));
  slpdir1(dimension(3, 12));
  slpnor1(dimension(3, 12));
  //C
  int i = fem::int0;
  int k = fem::int0;
  int l = fem::int0;
  FEM_DO_SAFE(i, 1, 12) {
    tau(i) = 0.0f;
    FEM_DO_SAFE(k, 1, 3) {
      FEM_DO_SAFE(l, 1, 3) {
        tau(i) += slpdir1(l, i) * slpnor1(k, i) * s(l, k);
      }
    }
  }
  //C
}

void
gd0(
  arr_ref<fem::real_star_8> gdot0,
  fem::real_star_8 const& d,
  fem::real_star_8 const& etha,
  arr_cref<fem::real_star_8> ro,
  arr_cref<fem::real_star_8> cts)
{
  gdot0(dimension(12));
  ro(dimension(12));
  cts(dimension(48));
  //C
  int i = fem::int0;
  FEM_DO_SAFE(i, 1, 12) {
    gdot0(i) = ro(i) * etha * d * cts(9) * cts(28);
  }
  //C
}

void
vgradient(
  arr_ref<fem::real_star_8, 2> vg,
  arr_cref<fem::real_star_8> gdot,
  arr_cref<fem::real_star_8, 2> slpdir0,
  arr_cref<fem::real_star_8, 2> slpnor0)
{
  vg(dimension(3, 3));
  gdot(dimension(12));
  slpdir0(dimension(3, 12));
  slpnor0(dimension(3, 12));
  //C
  int i = fem::int0;
  int k = fem::int0;
  int n = fem::int0;
  FEM_DO_SAFE(i, 1, 3) {
    FEM_DO_SAFE(k, 1, 3) {
      vg(i, k) = 0.0f;
      FEM_DO_SAFE(n, 1, 12) {
        vg(i, k) += slpdir0(i, n) * slpnor0(k, n) * gdot(n);
      }
    }
  }
  //C
}

void
inverse3(
  arr_cref<fem::real_star_8, 2> a,
  arr_ref<fem::real_star_8, 2> b)
{
  a(dimension(3, 3));
  b(dimension(3, 3));
  //C
  b(1, 1) = a(2, 2) * a(3, 3) - a(3, 2) * a(2, 3);
  b(1, 2) = a(3, 2) * a(1, 3) - a(1, 2) * a(3, 3);
  b(1, 3) = a(1, 2) * a(2, 3) - a(2, 2) * a(1, 3);
  b(2, 1) = a(3, 1) * a(2, 3) - a(2, 1) * a(3, 3);
  b(2, 2) = a(1, 1) * a(3, 3) - a(3, 1) * a(1, 3);
  b(2, 3) = a(2, 1) * a(1, 3) - a(1, 1) * a(2, 3);
  b(3, 1) = a(2, 1) * a(3, 2) - a(3, 1) * a(2, 2);
  b(3, 2) = a(3, 1) * a(1, 2) - a(1, 1) * a(3, 2);
  b(3, 3) = a(1, 1) * a(2, 2) - a(2, 1) * a(1, 2);
  fem::real_star_8 det = a(1, 1) * b(1, 1) + a(1, 2) * b(2, 1) + a(1,
    3) * b(3, 1);
  //C
  arr_2d<3, 3, fem::real_star_8> unit(fem::fill0);
  //unit = 0.0f;
  int i = fem::int0;
  int k = fem::int0;
  FEM_DO_SAFE(i, 1, 3) {
    unit(i, i) = 2.0f;
  }
    FEM_DO_SAFE(i, 1, 3)
    {
        FEM_DO_SAFE(k, 1, 3) {
            b(i, k) = b(i,k)/det;
        }
    }
  b = matmult(b, (unit - matmult(a, b)));
  //C
}

void
plasticgradientinv(
  arr_ref<fem::real_star_8, 2> fpinv,
  arr_cref<fem::real_star_8, 2> fpinv0,
  arr_ref<fem::real_star_8, 2> vg,
  fem::real_star_8 const& dtincr)
{
  fpinv(dimension(3, 3));
  fpinv0(dimension(3, 3));
  vg(dimension(3, 3));
  //C
  int i = fem::int0;
  int k = fem::int0;
  FEM_DO_SAFE(i, 1, 3) {
    FEM_DO_SAFE(k, 1, 3) {
      vg(i, k) = vg(i, k) * dtincr;
    }
  }
  fem::real_star_8 a1 = 0.0f;
  FEM_DO_SAFE(i, 1, 3) {
    FEM_DO_SAFE(k, 1, 3) {
      a1 += vg(i, k) * vg(i, k);
    }
  }
  //C
  arr_2d<3, 3, fem::real_star_8> a(fem::fill0);
  int n = fem::int0;
  FEM_DO_SAFE(i, 1, 3) {
    FEM_DO_SAFE(k, 1, 3) {
      a(i, k) = 0.0f;
      FEM_DO_SAFE(n, 1, 3) {
        a(i, k) += vg(i, n) * vg(n, k);
      }
    }
  }
  //C
  a1 = fem::sqrt(0.5f * a1);
  FEM_DO_SAFE(i, 1, 3) {
    if (fem::abs(a1) > 1e-8f) {
      FEM_DO_SAFE(k, 1, 3) {
        vg(i, k) = vg(i, k) * (fem::sin(a1) / a1) + a(i, k) * ((
          1.0f - fem::cos(a1)) / fem::pow(a1, 2.0f));
      }
    }
  }
  vg(1, 1) += 1.0f;
  vg(2, 2) += 1.0f;
  vg(3, 3) += 1.0f;
  //C
  arr_2d<3, 3, fem::real_star_8> a3(fem::fill0);
  inverse3(vg, a3);
  //C
  FEM_DO_SAFE(i, 1, 3) {
    FEM_DO_SAFE(k, 1, 3) {
      fpinv(i, k) = 0;
      FEM_DO_SAFE(n, 1, 3) {
        fpinv(i, k) += fpinv0(i, n) * a3(n, k);
      }
    }
  }
  //C
}

void
elasticgradient(
  arr_ref<fem::real_star_8, 2> fel,
  arr_cref<fem::real_star_8, 2> f1,
  arr_cref<fem::real_star_8, 2> fpinv)
{
  fel(dimension(3, 3));
  f1(dimension(3, 3));
  fpinv(dimension(3, 3));
  //C
  int i = fem::int0;
  int k = fem::int0;
  int n = fem::int0;
  FEM_DO_SAFE(i, 1, 3) {
    FEM_DO_SAFE(k, 1, 3) {
      fel(i, k) = 0;
      FEM_DO_SAFE(n, 1, 3) {
        fel(i, k) += f1(i, n) * fpinv(n, k);
      }
    }
  }
  //C
}

void
backstress(
  arr_ref<fem::real_star_8> b,
  arr_cref<fem::real_star_8> b0,
  arr_cref<fem::real_star_8> gdot,
  arr_cref<fem::real_star_8> cts,
  fem::real_star_8 const& dtincr,
  fem::real_star_8 const& fw,
  arr_cref<fem::real_star_8> statev,
  fem::real_star_8 const& etha,
  int const& /* kstep */)
{
  b(dimension(12));
  b0(dimension(12));
  gdot(dimension(12));
  cts(dimension(48));
  statev(dimension(400));
  //C
  int i = fem::int0;
  arr_1d<12, fem::real_star_8> fshill(fem::fill0);
  FEM_DO_SAFE(i, 1, 12) {
    fshill(i) = statev(i + 328);
  }
  //C
  fem::real_star_8 s1212 = fem::zero<fem::real_star_8>();
  if (cts(44) > 1.0f) {
    if (etha == 1.0f) {
      s1212 = cts(16);
    }
    else {
      s1212 = (0.32f / 0.24f) * cts(16);
    }
  }
  else {
    s1212 = cts(16);
  }
  //C
  fem::real_star_8 a0 = fem::zero<fem::real_star_8>();
  fem::real_star_8 a4 = fem::zero<fem::real_star_8>();
  FEM_DO_SAFE(i, 1, 12) {
    a0 = 1.0f + 4.0f * s1212 * cts(8) * fshill(i);
    a4 = 2.0f * cts(8) * (1.0f - 2.0f * s1212) / a0;
    //C
    b(i) = b0(i) + dtincr * (fw / (1.0f - fw)) * a4 * gdot(i);
  }
  //C
}

void
dannunload(
  arr_cref<fem::real_star_8> tau,
  arr_cref<fem::real_star_8> tau0,
  arr_ref<fem::real_star_8> statev,
  arr_cref<fem::real_star_8> b,
  arr_cref<fem::real_star_8> b0,
  fem::real_star_8 const& d,
  arr_cref<fem::real_star_8> cts,
  fem::real_star_8 const& dtincr)
{
  tau(dimension(12));
  tau0(dimension(12));
  statev(dimension(400));
  b(dimension(12));
  b0(dimension(12));
  cts(dimension(48));
  //C
  int i = fem::int0;
  fem::real_star_8 tdot = fem::zero<fem::real_star_8>();
  fem::real_star_8 bdot = fem::zero<fem::real_star_8>();
  arr_1d<12, fem::real_star_8> uannro(fem::fill0);
  FEM_DO_SAFE(i, 1, 12) {
    tdot = (tau(i) - tau0(i)) / dtincr;
    bdot = (b(i) - b0(i)) / dtincr;
    //C
    if (tau(i) - b(i) == 0.0f) {
      uannro(i) = 0.0f;
    }
    else {
      //C
      if (fem::sign(1.e0, tau(i) * tau0(i)) < 0.0f) {
        uannro(i) = 2.0f * 3.14f * fem::abs(tdot - bdot) / (d * cts(
          8) * cts(9));
      }
      else {
        uannro(i) = 0.0f;
      }
      //C
    }
  }
  //C
  FEM_DO_SAFE(i, 1, 12) {
    statev(i + 85) = fem::abs(uannro(i));
  }
  //C
}

void
ddensity(
  arr_ref<fem::real_star_8> ro,
  arr_cref<fem::real_star_8> ro0,
  arr_cref<fem::real_star_8> statev,
  fem::real_star_8 const& d,
  fem::real_star_8 const& etha,
  arr_cref<fem::real_star_8> gdot,
  arr_cref<fem::real_star_8> tau,
  arr_cref<fem::real_star_8> b,
  arr_cref<fem::real_star_8> cts,
  fem::real_star_8 const& dtincr,
  arr_cref<fem::real_star_8> ncrss)
{
  ro(dimension(12));
  ro0(dimension(12));
  statev(dimension(400));
  gdot(dimension(12));
  tau(dimension(12));
  b(dimension(12));
  cts(dimension(48));
  ncrss(dimension(12));
  //C
  fem::real_star_8 dkm = fem::zero<fem::real_star_8>();
  int i = fem::int0;
  arr_1d<12, fem::real_star_8> uannro(fem::fill0);
  fem::real_star_8 ethai = fem::zero<fem::real_star_8>();
  fem::real_star_8 ethaii = fem::zero<fem::real_star_8>();
  fem::real_star_8 etha_vein = fem::zero<fem::real_star_8>();
  fem::real_star_8 etha_psb = fem::zero<fem::real_star_8>();
  fem::real_star_8 etha_cell = fem::zero<fem::real_star_8>();
  fem::real_star_8 etha_lab = fem::zero<fem::real_star_8>();
  fem::real_star_8 r = fem::zero<fem::real_star_8>();
  if (cts(44) > 1.0f) {
    if (etha == 1.0f) {
      dkm = cts(10);
    }
    else {
      dkm = 2.0f * cts(10);
    }
    //C
    FEM_DO_SAFE(i, 1, 12) {
      if (ro0(i) > 7.0e13f) {
        uannro(i) = statev(i + 85);
      }
      else {
        uannro(i) = 0.0f;
      }
    }
    //C
    ethai = cts(36);
    ethaii = cts(37);
    etha_vein = cts(40);
    etha_psb = cts(41);
    etha_cell = cts(42);
    etha_lab = cts(43);
    //C
    if (etha == etha_vein) {
      r = 200.0f;
    }
    else if (etha == etha_psb) {
      r = 100.0f;
    }
    else if (etha == etha_lab) {
      r = 20.0f;
    }
    else if (etha == ethai) {
      r = 20.0f;
    }
    else if (etha == ethaii) {
      r = 20.0f;
    }
    else if (etha == etha_cell) {
      r = 10.0f;
    }
    //C
  }
  else {
    dkm = cts(10);
    FEM_DO_SAFE(i, 1, 12) {
      uannro(i) = 0.0f;
    }
    r = 1.0f;
  }
  //C
  fem::real_star_8 a1 = 1.0f * d / 12.0e-6f;
  fem::real_star_8 taucs = cts(8) * cts(9) / (4.0f * 3.14f * cts(12));
  //C
  arr_1d<12, fem::real_star_8> a2(fem::fill0);
  arr_1d<12, fem::real_star_8> a4(fem::fill0);
  FEM_DO_SAFE(i, 1, 12) {
    if (fem::abs(tau(i)) > 1.0e6f) {
      a2(i) = -cts(30) / (fem::abs(tau(i)) * cts(22) * cts(23));
    }
    else {
      a2(i) = 0.0f;
    }
    a4(i) = a2(i) * (taucs - fem::abs(tau(i) - b(i)));
  }
  //C
  fem::real_star_8 a5 = fem::zero<fem::real_star_8>();
  fem::real_star_8 a6 = fem::zero<fem::real_star_8>();
  fem::real_star_8 term1 = fem::zero<fem::real_star_8>();
  fem::real_star_8 term2 = fem::zero<fem::real_star_8>();
  fem::real_star_8 term4 = fem::zero<fem::real_star_8>();
  fem::real_star_8 romax = fem::zero<fem::real_star_8>();
  FEM_DO_SAFE(i, 1, 12) {
    a5 = 0.0f;
    a6 = 0.0f;
    if (ro0(i) > 1.1f * cts(31)) {
      if (ro0(ncrss(i)) > 1.1f * cts(31)) {
        a5 = ro0(ncrss(i)) * fem::exp(a4(i));
        a6 = ro0(i) * fem::exp(a4(ncrss(i)));
      }
    }
    //C
    term1 = (dkm / (etha * d * cts(9))) * fem::abs(gdot(i));
    term2 = (2.0f * cts(11) / cts(9)) * ro0(i) * fem::abs(gdot(i));
    term4 = a1 * (cts(29) * a5 - (1.0f - cts(29)) * a6);
    //C
    romax = (dkm / (etha * 1.0e-6f * 2.0f * cts(11)));
    //C
    ro(i) = ro0(i) + r * (dtincr * (term1 - term2 + term4 - uannro(i)));
    //C
    if (ro(i) < cts(31)) {
      ro(i) = cts(31);
    }
    if (ro(i) > romax) {
      ro(i) = romax;
    }
    //C
  }
  //C
}

void
threshstress(
  arr_ref<fem::real_star_8> ts,
  fem::real_star_8 const& d,
  arr_cref<fem::real_star_8> ro,
  arr_cref<fem::real_star_8> cts)
{
  ts(dimension(12));
  ro(dimension(12));
  cts(dimension(48));
  //C
  fem::real_star_8 a6 = cts(13) * cts(8) * cts(9);
  int i = fem::int0;
  FEM_DO_SAFE(i, 1, 12) {
    ts(i) = cts(15) + (a6 / (2.0f * d)) + cts(8) * cts(9) * fem::sqrt(
      cts(14) * ro(i));
  }
  //C
}

void
nrfunction(
  arr_ref<fem::real_star_8> func,
  arr_cref<fem::real_star_8> gdot,
  arr_cref<fem::real_star_8> gdot0,
  arr_cref<fem::real_star_8> tau,
  arr_cref<fem::real_star_8> b,
  arr_cref<fem::real_star_8> ts,
  arr_cref<fem::real_star_8> /* ro */,
  fem::real_star_8 const& /* d */,
  fem::real_star_8& sse,
  arr_cref<fem::real_star_8> cts,
  int& l1)
{
  func(dimension(12));
  gdot(dimension(12));
  gdot0(dimension(12));
  tau(dimension(12));
  b(dimension(12));
  ts(dimension(12));
  cts(dimension(48));
  //C
  fem::real_star_8 a1 = -cts(21) / (cts(22) * cts(23));
  fem::real_star_8 a2 = cts(26) * (cts(8) / cts(27));
  l1 = 0.0f;
  int i = fem::int0;
  arr_1d<12, fem::real_star_8> taueff(fem::fill0);
  arr_1d<12, fem::real_star_8> gau(fem::fill0);
  FEM_DO_SAFE(i, 1, 12) {
    taueff(i) = fem::abs(tau(i) - b(i)) - ts(i);
    //C
    if (taueff(i) <= 0.0f) {
      func(i) = gdot(i) - gdot0(i) * fem::exp(a1) * fem::sign(1.e0,
        tau(i) - b(i));
    }
    else {
      gau(i) = 1.0f - fem::pow(((taueff(i)) / a2), cts(24));
      if (gau(i) <= 0.0f) {
        func(i) = gdot(i) - gdot0(i) * fem::sign(1.e0, tau(i) - b(i));
        l1 = 1.0f;
      }
      else {
        func(i) = gdot(i) - gdot0(i) * fem::exp(a1 * (fem::pow((gau(i)),
          cts(25)))) * fem::sign(1.e0, tau(i) - b(i));
      }
    }
  }
  //C
  sse = 0.0f;
  fem::real_star_8 gmax = fem::abs(gdot(1));
  int k = fem::int0;
  FEM_DO_SAFE(k, 1, 12) {
    gmax = fem::max(fem::abs(gdot(k)), gmax);
  }
  gmax = fem::pow(gmax, 2.0f);
  FEM_DO_SAFE(k, 1, 12) {
    sse += fem::pow((gdot(k) * func(k)), 2.0f);
  }
  if (gmax > 0) {
    sse = fem::sqrt(sse / gmax) / 12.0f;
  }
  //C
}

void
check(
  arr_ref<fem::real_star_8> gdot,
  arr_cref<fem::real_star_8, 2> slpdir0,
  arr_cref<fem::real_star_8, 2> slpnor0,
  fem::real_star_8 const& dtincr,
  arr_cref<fem::real_star_8, 2> f1,
  arr_cref<fem::real_star_8, 4> c,
  arr_ref<fem::real_star_8> statev,
  arr_cref<fem::real_star_8> cts,
  fem::real_star_8& sse,
  arr_cref<fem::real_star_8> grad,
  arr_cref<fem::real_star_8> func0,
  arr_cref<fem::real_star_8, 2> fpinv0,
  fem::real_star_8& sseold,
  fem::real_star_8& sseref,
  arr_ref<fem::real_star_8, 2> fel,
  arr_ref<fem::real_star_8, 2> spk2,
  arr_ref<fem::real_star_8> tau,
  arr_cref<fem::real_star_8> tau0,
  arr_ref<fem::real_star_8> b,
  arr_ref<fem::real_star_8> ts,
  arr_ref<fem::real_star_8, 2> fpinv,
  arr_ref<fem::real_star_8> ro,
  arr_ref<fem::real_star_8, 2> s,
  arr_cref<fem::real_star_8> ro0,
  arr_cref<fem::real_star_8> b0,
  fem::real_star_8 const& etha,
  fem::real_star_8& l1,
  fem::real_star_8 const& fw,
  int const& /* nincrt */,
  fem::real_star_8 const& /* dtime */,
  fem::real_star_8 const& d,
  arr_cref<fem::real_star_8> ncrss,
  int const& kstep,
  int const& /* m1 */)
{
  gdot(dimension(12));
  slpdir0(dimension(3, 12));
  slpnor0(dimension(3, 12));
  f1(dimension(3, 3));
  c(dimension(3, 3, 3, 3));
  statev(dimension(400));
  cts(dimension(48));
  grad(dimension(12));
  func0(dimension(12));
  fpinv0(dimension(3, 3));
  fel(dimension(3, 3));
  spk2(dimension(3, 3));
  tau(dimension(12));
  tau0(dimension(12));
  b(dimension(12));
  ts(dimension(12));
  fpinv(dimension(3, 3));
  ro(dimension(12));
  s(dimension(3, 3));
  ro0(dimension(12));
  b0(dimension(12));
  ncrss(dimension(12));
  int k = fem::int0;
  arr_1d<12, fem::real_star_8> dgdot(fem::fill0);
  fem::real_star_8 sum = fem::zero<fem::real_star_8>();
  int i = fem::int0;
  bool improved = fem::bool0;
  arr_1d<12, fem::real_star_8> gtray(fem::fill0);
  arr_2d<3, 3, fem::real_star_8> vg(fem::fill0);
  arr_2d<3, 3, fem::real_star_8> eel(fem::fill0);
  arr_2d<3, 3, fem::real_star_8> felinv(fem::fill0);
  arr_2d<3, 12, fem::real_star_8> slpdir1(fem::fill0);
  arr_2d<3, 12, fem::real_star_8> slpnor1(fem::fill0);
  arr_1d<12, fem::real_star_8> gdot0(fem::fill0);
  arr_1d<12, fem::real_star_8> func(fem::fill0);
  //CREAL*8 (A-H,O-Z)
  //C
  sseref = sse;
  FEM_DO_SAFE(k, 1, 12) {
    dgdot(k) = func0(k);
  }
  //C: CHECKING THE DOWNHILL OF NR
  sum = 0.0f;
  FEM_DO_SAFE(i, 1, 12) {
    sum = sum - grad(i) * dgdot(i);
  }
  if (sum > 0.0f) {
    FEM_DO_SAFE(i, 1, 12) {
      dgdot(i) = -dgdot(i);
    }
  }
  FEM_DO_SAFE(i, 1, 12) {
    dgdot(i) = dgdot(i) * 2.0f;
  }
  improved = false;
  FEM_DO_SAFE(i, 1, 10) {
    sseold = sse;
    FEM_DO_SAFE(k, 1, 12) {
      dgdot(k) = dgdot(k) / 2.0f;
      gtray(k) = gdot(k) + dgdot(k);
    }
    //C
    //C: L VELOCITY GRADIENT  SUM OF GDOT*S0*M0
    vgradient(vg, gtray, slpdir0, slpnor0);
    //C:  PLASTIC PART OF DEFORMATION GRADIENT  FP(t+dt) AND FINALLY INVERSE OF PLASTIC GRADIENT  FP-1(t+dt)
    plasticgradientinv(fpinv, fpinv0, vg, dtincr);
    //C: ELASTIC PART OF DEFORMATION GRADIENT  FE(t+dt)
    elasticgradient(fel, f1, fpinv);
    //C: GREEN STRAIN TENSOR
    green(eel, fel);
    //C: SECOND PIOLA-KIRCHHOF STRESS
    piola(spk2, c, eel);
    //C: INVERSE OF ELASTIC GRADIENT  FEL-1(t+dt)
    inverse3(fel, felinv);
    //C: CAUCHY STRESS
    cauchy(s, fel, felinv, spk2);
    //C: CURRENT SLIP SYSTEMS
    slip1(slpdir1, slpnor1, slpdir0, slpnor0, fel, felinv);
    //C: RESOLVED SHEAR STRESS
    rshrs(tau, s, slpdir0, slpnor0);
    //C: B BACK STRESS       +
    backstress(b, b0, gtray, cts, dtincr, fw, statev, etha, kstep);
    //C: UNLOAD DISLOCATION ANNIHILATION
    dannunload(tau, tau0, statev, b, b0, d, cts, dtincr);
    //C: DISLOCATION DENSITY
    ddensity(ro, ro0, statev, d, etha, gtray, tau, b, cts, dtincr, ncrss);
    //C: TS: THRESHOLD STRESS, OR SLIP RESISTANCE
    threshstress(ts, d, ro, cts);
    //C: GAMA-DOT0   RO*L*B*Vg
    gd0(gdot0, d, etha, ro, cts);
    //C:  FUNCTION FOR NEWTON RAPHSON
    nrfunction(func, gtray, gdot0, tau, b, ts, ro, d, sse, cts, l1);
    //C: ANALYTICAL JACOBIAN MATRIX FOR NEWTON RAPHSON
    if ((sseold <= sseref) && (sse >= sseold) && (i > 1.0f)) {
      improved = true;
    }
    if (improved) {
      goto statement_200;
    }
  }
  statement_200:
  FEM_DO_SAFE(k, 1, 12) {
    gdot(k) += dgdot(k) * 2.0f;
  }
  //C
}

void
structmon(
  arr_ref<fem::real_star_8> statev,
  fem::real_star_8& etha,
  fem::real_star_8& fw,
  arr_cref<fem::real_star_8> cts)
{
  statev(dimension(400));
  cts(dimension(48));
  //CREAL*8 (A-H,O-Z)
  //C
  int i = fem::int0;
  arr_1d<12, fem::real_star_8> eslip(fem::fill0);
  FEM_DO_SAFE(i, 1, 12) {
    eslip(i) = fem::abs(statev(i + 136));
  }
  //C
  fem::real_star_8 dgammatx = 0.0f;
  //C
  dgammatx = fem::max(eslip(1), eslip(2), eslip(3), eslip(4), eslip(5),
    eslip(6), eslip(7), eslip(8), eslip(9), eslip(10), eslip(11), eslip(
    12));
  //C
  //C SPECIFY THE MULTI, DOUBLE OR SINGLE SLIP SYSTEMS
  //C
  arr_1d<12, fem::real_star_8> t(fem::fill0);
  FEM_DO_SAFE(i, 1, 12) {
    t(i) = eslip(i) / dgammatx;
  }
  //C
  int k = 0;
  FEM_DO_SAFE(i, 1, 12) {
    if (t(i) > 0.9f) {
      k++;
    }
  }
  dgammatx = dgammatx * k;
  //C
  fem::real_star_8 ecross = 0.0f;
  fem::real_star_8 ehirth = 0.0f;
  fem::real_star_8 ecot = 0.0f;
  //C
  if (t(1) > 0.9f) {
    ecross += eslip(4);
  }
  //C
  if (t(2) > 0.9f) {
    ecross += eslip(11);
  }
  //C
  if (t(3) > 0.9f) {
    ecross += eslip(9);
  }
  //C
  if (t(4) > 0.9f) {
    ecross += eslip(1);
  }
  //C
  if (t(5) > 0.9f) {
    ecross += eslip(8);
  }
  //C
  if (t(6) > 0.9f) {
    ecross += eslip(12);
  }
  //C
  if (t(7) > 0.9f) {
    ecross += eslip(10);
  }
  //C
  if (t(8) > 0.9f) {
    ecross += eslip(5);
  }
  //C
  if (t(9) > 0.9f) {
    ecross += eslip(3);
  }
  //C
  if (t(10) > 0.9f) {
    ecross += eslip(7);
  }
  //C
  if (t(11) > 0.9f) {
    ecross += eslip(2);
  }
  //C
  if (t(12) > 0.9f) {
    ecross += eslip(6);
  }
  //C
  if (ecross < statev(81)) {
    ecross = statev(81);
  }
  //C
  fem::real_star_8 gp = cts(18);
  fem::real_star_8 bf0 = cts(19);
  fem::real_star_8 bfi = cts(20);
  fem::real_star_8 fwi = cts(33);
  fem::real_star_8 fwii = cts(34);
  fem::real_star_8 fwiii = cts(35);
  fem::real_star_8 ethai = cts(36);
  fem::real_star_8 ethaii = cts(37);
  fem::real_star_8 gammacr = cts(15);
  //C
  if (dgammatx == 0.0f) {
    dgammatx = 100.0f;
  }
  //C
  if (ecross > gammacr) {
    etha = 1.0f;
  }
  else {
    if (ecot > ehirth) {
      etha = ethai;
    }
    else {
      etha = ethaii;
    }
  }
  //C
  if (etha == 1.0f) {
    fw = cts(35);
  }
  else {
    fw = bfi + (bf0 - bfi) * fem::exp(-dgammatx / gp);
  }
  //C
  statev(81) = ecross;
  statev(82) = ecot;
  statev(83) = ehirth;
  statev(84) = etha;
  statev(85) = fw;
  //C
}

void
structcyc(
  arr_ref<fem::real_star_8> statev,
  fem::real_star_8& etha,
  fem::real_star_8& fw,
  arr_cref<fem::real_star_8> cts,
  int const& kstep,
  arr_cref<fem::real_star_8> props)
{
  statev(dimension(400));
  cts(dimension(48));
  props(dimension(5));
  //C
  arr_1d<12, fem::real_star_8> nhirth1(fem::fill0);
  nhirth1(1) = 7.0f;
  arr_1d<12, fem::real_star_8> nhirth2(fem::fill0);
  nhirth2(1) = 10.0f;
  nhirth1(2) = 5.0f;
  nhirth2(2) = 8.0f;
  nhirth1(3) = 6.0f;
  nhirth2(3) = 12.0f;
  nhirth1(4) = 7.0f;
  nhirth2(4) = 10.0f;
  nhirth1(5) = 2.0f;
  nhirth2(5) = 11.0f;
  nhirth1(6) = 3.0f;
  nhirth2(6) = 9.0f;
  nhirth1(7) = 1.0f;
  nhirth2(7) = 4.0f;
  nhirth1(8) = 2.0f;
  nhirth2(8) = 11.0f;
  nhirth1(9) = 6.0f;
  nhirth2(9) = 12.0f;
  nhirth1(10) = 1.0f;
  nhirth2(10) = 4.0f;
  nhirth1(11) = 5.0f;
  nhirth2(11) = 8.0f;
  nhirth1(12) = 3.0f;
  nhirth2(12) = 9.0f;
  arr_1d<12, fem::real_star_8> ncrss(fem::fill0);
  ncrss(1) = 4.0f;
  ncrss(2) = 11.0f;
  ncrss(3) = 9.0f;
  ncrss(4) = 1.0f;
  ncrss(5) = 8.0f;
  ncrss(6) = 12.0f;
  ncrss(7) = 10.0f;
  ncrss(8) = 5.0f;
  ncrss(9) = 3.0f;
  ncrss(10) = 7.0f;
  ncrss(11) = 2.0f;
  ncrss(12) = 6.0f;
  //C
  int i = fem::int0;
  arr_1d<12, fem::real_star_8> erange(fem::fill0);
  FEM_DO_SAFE(i, 1, 12) {
    erange(i) = statev(i + 136);
  }
  //C
  fem::real_star_8 dgmax = 0.0f;
  FEM_DO_SAFE(i, 1, 12) {
    dgmax = fem::max(erange(i), dgmax);
  }
  //C
  fem::real_star_8 decross = 0.0f;
  fem::real_star_8 dehirth = 0.0f;
  fem::real_star_8 decot = 0.0f;
  int nmax2 = fem::int0;
  FEM_DO_SAFE(i, 1, 12) {
    if (dgmax == erange(i)) {
      decross += erange(ncrss(i));
      dehirth += erange(nhirth1(i)) + erange(nhirth2(i));
      nmax2 = i;
    }
  }
  //C
  fem::real_star_8 dgmax2 = 0.0f;
  FEM_DO_SAFE(i, 1, 12) {
    if (i != nmax2) {
      dgmax2 = fem::max(erange(i), dgmax2);
    }
  }
  //C
  int nmax3 = fem::int0;
  FEM_DO_SAFE(i, 1, 12) {
    if (dgmax2 == erange(i)) {
      decross += erange(ncrss(i));
      dehirth += erange(nhirth1(i)) + erange(nhirth2(i));
      nmax3 = i;
    }
  }
  //C
  fw = statev(85);
  fem::real_star_8 gp = cts(18);
  fem::real_star_8 bf0 = cts(19);
  fem::real_star_8 bfi = cts(20);
  fem::real_star_8 ethai = cts(36);
  fem::real_star_8 ethaii = cts(37);
  fem::real_star_8 dgm = cts(38);
  fem::real_star_8 dgmpsb = cts(39);
  fem::real_star_8 etha_vein = cts(40);
  fem::real_star_8 etha_psb = cts(41);
  fem::real_star_8 etha_cell = cts(42);
  fem::real_star_8 etha_lab = cts(43);
  //C
  fem::real_star_8 fw0 = statev(85);
  fem::real_star_8 etha0 = statev(84);
  //C
  if (dgmax < dgm) {
    etha = etha_vein;
  }
  else {
    if (decross < dgm && dgmax < dgmpsb) {
      etha = etha_psb;
    }
    else {
      if (dehirth <= 0.1f * dgm) {
        etha = etha_cell;
      }
      else {
        etha = etha_lab;
      }
    }
  }
  //C
  if (etha > etha0 && fw0 > 0.05f) {
    etha = etha0;
  }
  //C
  //C   OVERLOAD EFFECT
  if (kstep >= cts(46) && cts(46) != 0.0f) {
    if (props(4) == 1.0f) {
      etha = props(4);
    }
    else if (props(4) < (ethai + 0.001f)) {
      etha = props(4);
    }
    else if (props(4) < (ethaii + 0.001f)) {
      etha = props(4);
    }
  }
  //C
  if (etha == 1.0f) {
    fw = cts(35);
  }
  else {
    fw = bfi + (bf0 - bfi) * fem::exp(-dgmax / gp);
  }
  //C
  statev(81) = decross;
  statev(82) = decot;
  statev(83) = dehirth;
  statev(84) = etha;
  statev(85) = fw;
  //C
}

void
rotation(
  arr_ref<fem::real_star_8, 2> rotate,
  arr_cref<fem::real_star_8> cts)
{
  rotate(dimension(3, 3));
  cts(dimension(48));
  //C
  fem::real_star_8 sum1 = fem::sqrt(fem::pow2(cts(1)) + fem::pow2(cts(
    2)) + fem::pow2(cts(3)));
  fem::real_star_8 sum2 = fem::sqrt(fem::pow2(cts(2)) + fem::pow2(cts(1)));
  arr_2d<3, 3, fem::real_star_8> term1(fem::fill0);
  term1(1, 1) = cts(1) / sum1;
  term1(2, 1) = cts(2) / sum1;
  term1(3, 1) = cts(3) / sum1;
  //C
  term1(1, 2) = -cts(2) / sum2;
  term1(2, 2) = cts(1) / sum2;
  term1(3, 2) = 0.0f;
  //C
  term1(1, 3) = term1(2, 1) * term1(3, 2) - term1(3, 1) * term1(2, 2);
  term1(2, 3) = term1(3, 1) * term1(1, 2) - term1(1, 1) * term1(3, 2);
  term1(3, 3) = term1(1, 1) * term1(2, 2) - term1(2, 1) * term1(1, 2);
  int i = fem::int0;
  int j = fem::int0;
  FEM_DO_SAFE(i, 1, 3) {
    FEM_DO_SAFE(j, 1, 3) {
      rotate(i, j) = term1(j, i);
    }
  }
  //C
}

void
rotationpoly(
  arr_ref<fem::real_star_8, 2> rotate,
  arr_cref<fem::real_star_8> cts,
  arr_ref<fem::real_star_8> statev)
{
  rotate(dimension(3, 3));
  cts(dimension(48));
  statev(dimension(400));
  //C
  fem::real_star_8 phi_1 = cts(1);
  fem::real_star_8 phi = cts(2);
  fem::real_star_8 phi_2 = cts(3);
  //C
  fem::real_star_8 s1 = fem::sin(phi_1);
  fem::real_star_8 c1 = fem::cos(phi_1);
  fem::real_star_8 s2 = fem::sin(phi);
  fem::real_star_8 c2 = fem::cos(phi);
  fem::real_star_8 s3 = fem::sin(phi_2);
  fem::real_star_8 c3 = fem::cos(phi_2);
  //C
  rotate(1, 1) = c1 * c3 - s1 * s3 * c2;
  rotate(2, 1) = s1 * c3 + c1 * s3 * c2;
  rotate(3, 1) = s3 * s2;
  rotate(1, 2) = -c1 * s3 - s1 * c3 * c2;
  rotate(2, 2) = -s1 * s3 + c1 * c3 * c2;
  rotate(3, 2) = c3 * s2;
  rotate(1, 3) = s1 * s2;
  rotate(2, 3) = -c1 * s2;
  rotate(3, 3) = c2;
  //C
  statev(350) = s1;
  statev(351) = c1;
  statev(352) = phi_1;
}

void
initial(
  arr_ref<fem::real_star_8> b0,
  arr_ref<fem::real_star_8> ro0,
  arr_ref<fem::real_star_8> tau0,
  arr_ref<fem::real_star_8> statev,
  arr_ref<fem::real_star_8, 2> fpinv0,
  arr_cref<fem::real_star_8> cts)
{
  b0(dimension(12));
  ro0(dimension(12));
  tau0(dimension(12));
  statev(dimension(400));
  fpinv0(dimension(3, 3));
  cts(dimension(48));
  //C
  int i = fem::int0;
  int k = fem::int0;
  FEM_DO_SAFE(i, 1, 3) {
    FEM_DO_SAFE(k, 1, 3) {
      fpinv0(i, k) = 0.0f;
    }
  }
  fpinv0(1, 1) = 1.0f;
  fpinv0(2, 2) = 1.0f;
  fpinv0(3, 3) = 1.0f;
  //C
  FEM_DO_SAFE(i, 1, 100) {
    statev(i) = 0.0f;
  }
  FEM_DO_SAFE(i, 1, 12) {
    b0(i) = 0.0f;
    tau0(i) = 0.0f;
  }
  FEM_DO_SAFE(i, 1, 12) {
    if (cts(44) == 1.0f) {
      ro0(i) = cts(31);
      statev(32 + i) = cts(31);
    }
    else {
      ro0(i) = cts(31);
      statev(32 + i) = cts(31);
    }
  }
  //C
}

void
firstread(
  arr_cref<fem::real_star_8> statev,
  arr_ref<fem::real_star_8, 2> fpinv0,
  arr_ref<fem::real_star_8> ro0,
  arr_ref<fem::real_star_8> b0,
  arr_ref<fem::real_star_8> ts,
  arr_ref<fem::real_star_8> tau0)
{
  statev(dimension(400));
  fpinv0(dimension(3, 3));
  ro0(dimension(12));
  b0(dimension(12));
  ts(dimension(12));
  tau0(dimension(12));
  //C
  int n = 0.0f;
  int i = fem::int0;
  int k = fem::int0;
  FEM_DO_SAFE(i, 1, 3) {
    FEM_DO_SAFE(k, 1, 3) {
      n++;
      fpinv0(i, k) = statev(n);
    }
  }
  n = 32.0f;
  FEM_DO_SAFE(i, 1, 12) {
    n++;
    ro0(i) = statev(n);
  }
  n = 44.0f;
  FEM_DO_SAFE(i, 1, 12) {
    n++;
    b0(i) = statev(n);
  }
  FEM_DO_SAFE(i, 1, 12) {
    tau0(i) = statev(56 + i);
    ts(i) = statev(68 + i);
  }
  //C
}

void
elastic0(
  arr_ref<fem::real_star_8, 4> c0,
  arr_cref<fem::real_star_8> cts)
{
  c0(dimension(3, 3, 3, 3));
  cts(dimension(48));
  //C
  int k = fem::int0;
  int i = fem::int0;
  arr_2d<3, 3, fem::real_star_8> a1(fem::fill0);
  FEM_DO_SAFE(k, 1, 3) {
    FEM_DO_SAFE(i, 1, 3) {
      a1(k, i) = 0.0f;
    }
  }
  a1(1, 1) = 1.0f;
  a1(2, 2) = 1.0f;
  a1(3, 3) = 1.0f;
  //C
  int m = fem::int0;
  int n = fem::int0;
  FEM_DO_SAFE(i, 1, 3) {
    FEM_DO_SAFE(k, 1, 3) {
      FEM_DO_SAFE(m, 1, 3) {
        FEM_DO_SAFE(n, 1, 3) {
          c0(i, k, m, n) = cts(5) * a1(i, k) * a1(m, n) + cts(6) * (a1(i,
            m) * a1(k, n) + a1(i, n) * a1(m, k));
        }
      }
    }
  }
  c0(1, 1, 1, 1) = cts(4);
  c0(2, 2, 2, 2) = cts(4);
  c0(3, 3, 3, 3) = cts(4);
  //C
}

void
elastic(
  arr_ref<fem::real_star_8, 4> c,
  arr_cref<fem::real_star_8, 4> c0,
  arr_cref<fem::real_star_8, 2> r)
{
  c(dimension(3, 3, 3, 3));
  c0(dimension(3, 3, 3, 3));
  r(dimension(3, 3));
  //C
  int m = fem::int0;
  int n = fem::int0;
  int k = fem::int0;
  int l = fem::int0;
  arr<fem::real_star_8, 4> a(dimension(3, 3, 3, 3), fem::fill0);
  FEM_DO_SAFE(m, 1, 3) {
    FEM_DO_SAFE(n, 1, 3) {
      FEM_DO_SAFE(k, 1, 3) {
        FEM_DO_SAFE(l, 1, 3) {
          a(m, n, k, l) = r(k, 1) * (r(l, 1) * c0(m, n, 1, 1) + r(l,
            2) * c0(m, n, 1, 2) + r(l, 3) * c0(m, n, 1, 3)) + r(k,
            2) * (r(l, 1) * c0(m, n, 2, 1) + r(l, 2) * c0(m, n, 2,
            2) + r(l, 3) * c0(m, n, 2, 3)) + r(k, 3) * (r(l, 1) * c0(m,
            n, 3, 1) + r(l, 2) * c0(m, n, 3, 2) + r(l, 3) * c0(m, n,
            3, 3));
        }
      }
    }
  }
  int i = fem::int0;
  int j = fem::int0;
  FEM_DO_SAFE(i, 1, 3) {
    FEM_DO_SAFE(j, 1, 3) {
      FEM_DO_SAFE(k, 1, 3) {
        FEM_DO_SAFE(l, 1, 3) {
          c(i, j, k, l) = r(i, 1) * (r(j, 1) * a(1, 1, k, l) + r(j,
            2) * a(1, 2, k, l) + r(j, 3) * a(1, 3, k, l)) + r(i, 2) * (r(j,
            1) * a(2, 1, k, l) + r(j, 2) * a(2, 2, k, l) + r(j, 3) * a(2,
            3, k, l)) + r(i, 3) * (r(j, 1) * a(3, 1, k, l) + r(j,
            2) * a(3, 2, k, l) + r(j, 3) * a(3, 3, k, l));
        }
      }
    }
  }
  //C
}

void
slip0(
  arr_ref<fem::real_star_8, 2> slpdir0,
  arr_ref<fem::real_star_8, 2> slpnor0,
  arr_cref<fem::real_star_8, 2> rotate,
  arr_cref<fem::real_star_8> /* cts */,
  arr_ref<fem::real_star_8> ncrss)
{
  slpdir0(dimension(3, 12));
  slpnor0(dimension(3, 12));
  rotate(dimension(3, 3));
  ncrss(dimension(12));
  //C  THE SYSTEM OF SLIP SYSTEM NOTATION IS:  (111)[01-1],(111)[-101],(111)[1-10],(1-1-1)[0-11], (1-1-1)[-10-1],(1-1-1)[110],(-1-11)[0-1-1],(-1-11)[101],(-1-11)[-110],(-11-1)[011],(10-1)[101],(-11-1)[-1-10]
  //C
  int i = fem::int0;
  int k = fem::int0;
  FEM_DO_SAFE(i, 1, 12) {
    FEM_DO_SAFE(k, 1, 3) {
      slpdir0(k, i) = 1.0f / fem::sqrt(2.0f);
      slpnor0(k, i) = 1.0f / fem::sqrt(3.0f);
      if (i == k || i - k == 3) {
        slpdir0(k, i) = 0.0f;
      }
      if (i - k == 6) {
        slpdir0(k, i) = 0.0f;
      }
      if (i - k == 9) {
        slpdir0(k, i) = 0.0f;
      }
    }
  }
  slpdir0(1, 2) = -1.0f / fem::sqrt(2.0f);
  slpdir0(1, 5) = -1.0f / fem::sqrt(2.0f);
  slpdir0(1, 9) = -1.0f / fem::sqrt(2.0f);
  slpdir0(1, 12) = -1.0f / fem::sqrt(2.0f);
  slpdir0(2, 3) = -1.0f / fem::sqrt(2.0f);
  slpdir0(2, 4) = -1.0f / fem::sqrt(2.0f);
  slpdir0(2, 7) = -1.0f / fem::sqrt(2.0f);
  slpdir0(2, 12) = -1.0f / fem::sqrt(2.0f);
  slpdir0(3, 1) = -1.0f / fem::sqrt(2.0f);
  slpdir0(3, 5) = -1.0f / fem::sqrt(2.0f);
  slpdir0(3, 7) = -1.0f / fem::sqrt(2.0f);
  slpdir0(3, 11) = -1.0f / fem::sqrt(2.0f);
  //C
  slpnor0(1, 7) = -1.0f / fem::sqrt(3.0f);
  slpnor0(1, 8) = -1.0f / fem::sqrt(3.0f);
  slpnor0(1, 9) = -1.0f / fem::sqrt(3.0f);
  slpnor0(1, 10) = -1.0f / fem::sqrt(3.0f);
  slpnor0(1, 11) = -1.0f / fem::sqrt(3.0f);
  slpnor0(1, 12) = -1.0f / fem::sqrt(3.0f);
  slpnor0(2, 4) = -1.0f / fem::sqrt(3.0f);
  slpnor0(2, 5) = -1.0f / fem::sqrt(3.0f);
  slpnor0(2, 6) = -1.0f / fem::sqrt(3.0f);
  slpnor0(2, 7) = -1.0f / fem::sqrt(3.0f);
  slpnor0(2, 8) = -1.0f / fem::sqrt(3.0f);
  slpnor0(2, 9) = -1.0f / fem::sqrt(3.0f);
  slpnor0(3, 4) = -1.0f / fem::sqrt(3.0f);
  slpnor0(3, 5) = -1.0f / fem::sqrt(3.0f);
  slpnor0(3, 6) = -1.0f / fem::sqrt(3.0f);
  slpnor0(3, 10) = -1.0f / fem::sqrt(3.0f);
  slpnor0(3, 11) = -1.0f / fem::sqrt(3.0f);
  slpnor0(3, 12) = -1.0f / fem::sqrt(3.0f);
  //C
  arr_2d<3, 12, fem::real_star_8> nss(fem::fill0);
  FEM_DO_SAFE(i, 1, 12) {
    FEM_DO_SAFE(k, 1, 3) {
      nss(k, i) = 1.0f;
      if (i == k || i - k == 3) {
        nss(k, i) = 0.0f;
      }
      if (i - k == 6) {
        nss(k, i) = 0.0f;
      }
      if (i - k == 9) {
        nss(k, i) = 0.0f;
      }
    }
  }
  nss(1, 2) = -1.0f;
  nss(1, 5) = -1.0f;
  nss(1, 9) = -1.0f;
  nss(1, 12) = -1.0f;
  nss(2, 3) = -1.0f;
  nss(2, 4) = -1.0f;
  nss(2, 7) = -1.0f;
  nss(2, 12) = -1.0f;
  nss(3, 1) = -1.0f;
  nss(3, 5) = -1.0f;
  nss(3, 7) = -1.0f;
  nss(3, 11) = -1.0f;
  //C
  int j = fem::int0;
  FEM_DO_SAFE(i, 1, 12) {
    FEM_DO_SAFE(j, 1, 12) {
      if ((nss(1, i) * nss(1, j) + nss(2, i) * nss(2, j) + nss(3,
          i) * nss(3, j)) ==  - 2.0f) {
        //C
        ncrss(i) = j;
      }
    }
  }
  //C
  //C  SLIP SYSTEMS IN GLOBAL LOADING DIRECTION
  int l = fem::int0;
  arr_1d<3, fem::real_star_8> a1(fem::fill0);
  FEM_DO_SAFE(l, 1, 12) {
    FEM_DO_SAFE(i, 1, 3) {
      a1(i) = 0.f;
      FEM_DO_SAFE(k, 1, 3) {
        a1(i) += rotate(i, k) * slpdir0(k, l);
      }
    }
    FEM_DO_SAFE(i, 1, 3) {
      slpdir0(i, l) = a1(i);
    }
    FEM_DO_SAFE(i, 1, 3) {
      a1(i) = 0.f;
      FEM_DO_SAFE(k, 1, 3) {
        a1(i) += rotate(i, k) * slpnor0(k, l);
      }
    }
    FEM_DO_SAFE(i, 1, 3) {
      slpnor0(i, l) = a1(i);
    }
  }
  //C
}

void
dgrad(
  arr_ref<fem::real_star_8, 2> fel,
  arr_ref<fem::real_star_8, 2> f,
  arr_ref<fem::real_star_8, 2> f1,
  arr_cref<fem::real_star_8, 2> dfgrd0,
  arr_cref<fem::real_star_8, 2> dfgrd1,
  arr_cref<fem::real_star_8, 2> fpinv0,
  int const& nincr,
  int const& nincrt)
{
  fel(dimension(3, 3));
  f(dimension(3, 3));
  f1(dimension(3, 3));
  dfgrd0(dimension(3, 3));
  dfgrd1(dimension(3, 3));
  fpinv0(dimension(3, 3));
  //C
  int i = fem::int0;
  int k = fem::int0;
  FEM_DO_SAFE(i, 1, 3) {
    FEM_DO_SAFE(k, 1, 3) {
      f(i, k) = dfgrd0(i, k) + (dfgrd1(i, k) - dfgrd0(i, k)) *
        fem::real(fem::real(nincr - 1) / fem::real(nincrt));
      f1(i, k) = dfgrd0(i, k) + (dfgrd1(i, k) - dfgrd0(i, k)) *
        fem::real(fem::real(nincr) / fem::real(nincrt));
    }
  }
  //C
  int n = fem::int0;
  FEM_DO_SAFE(i, 1, 3) {
    FEM_DO_SAFE(k, 1, 3) {
      fel(i, k) = 0;
      FEM_DO_SAFE(n, 1, 3) {
        fel(i, k) += f(i, n) * fpinv0(n, k);
      }
    }
  }
  //C
}

void
cellsizemon(
  fem::real_star_8& d,
  arr_cref<fem::real_star_8> tau,
  arr_cref<fem::real_star_8> cts,
  arr_ref<fem::real_star_8> statev,
  int const& /* kstep */,
  arr_cref<fem::real_star_8> /* time */)
{
  tau(dimension(12));
  cts(dimension(46));
  statev(dimension(400));
  //C
  fem::real_star_8 taurmax = fem::max(fem::abs(tau(1)), fem::abs(tau(2)),
    fem::abs(tau(3)), fem::abs(tau(4)), fem::abs(tau(5)), fem::abs(tau(6)),
    fem::abs(tau(7)), fem::abs(tau(8)), fem::abs(tau(9)), fem::abs(tau(10)),
    fem::abs(tau(11)), fem::abs(tau(12)));
  //C
  taurmax = taurmax / 2.0f;
  //C
  fem::real_star_8 d0 = 12.0e-6f;
  fem::real_star_8 a5 = cts(45) * cts(8) * cts(9);
  d = a5 / taurmax;
  d = fem::min(d0, d);
  //C
  statev(100) = d;
  //C
}

void
cellsize(
  fem::real_star_8& d,
  arr_cref<fem::real_star_8> tau,
  arr_cref<fem::real_star_8> cts,
  arr_ref<fem::real_star_8> statev,
  int const& kstep,
  arr_cref<fem::real_star_8> time,
  arr_cref<fem::real_star_8> props)
{
  tau(dimension(12));
  cts(dimension(48));
  statev(dimension(400));
  time(dimension(2));
  props(dimension(5));
  //C
  int i = fem::int0;
  arr_1d<12, fem::real_star_8> ptau(fem::fill0);
  FEM_DO_SAFE(i, 1, 12) {
    ptau(i) = statev(i + 280);
  }
  //C
  arr_1d<12, fem::real_star_8> taurange(fem::fill0);
  arr_1d<12, fem::real_star_8> taumax(fem::fill0);
  arr_1d<12, fem::real_star_8> taumin(fem::fill0);
  FEM_DO_SAFE(i, 1, 12) {
    if (kstep < 3.0f && time(2) < 1.98f) {
      if (fem::abs(tau(i)) > fem::abs(ptau(i)) && kstep == 1.0f) {
        if (fem::sign(1.e0, tau(i) * ptau(i)) >= 0.0f) {
          taurange(i) = fem::abs(tau(i));
          statev(280 + i) = tau(i);
          taumax(i) = fem::max(0.0f, tau(i));
          taumin(i) = fem::min(0.0f, tau(i));
        }
        else {
          taumax(i) = fem::max(0.0f, tau(i));
          taumin(i) = fem::min(0.0f, tau(i));
          taurange(i) = fem::abs(taumax(i)) + fem::abs(taumin(i));
          statev(280 + i) = tau(i);
        }
      }
      else {
        //C
        if (fem::sign(1.e0, tau(i) * ptau(i)) < 0.0f) {
          taurange(i) = fem::abs(tau(i)) + fem::abs(ptau(i));
          statev(280 + i) = ptau(i);
          taumax(i) = fem::max(ptau(i), tau(i));
          taumin(i) = fem::min(ptau(i), tau(i));
        }
        else {
          taurange(i) = fem::abs(ptau(i));
          statev(280 + i) = ptau(i);
          taumax(i) = fem::max(0.0f, ptau(i));
          taumin(i) = fem::min(0.0f, ptau(i));
        }
      }
      statev(304 + i) = taumax(i);
      statev(316 + i) = taumin(i);
      statev(292 + i) = taurange(i);
    }
    else {
      taumax(i) = statev(304 + i);
      taumin(i) = statev(316 + i);
      if (tau(i) > taumax(i)) {
        taumax(i) = tau(i);
      }
      if (tau(i) < taumin(i)) {
        taumin(i) = tau(i);
      }
      statev(304 + i) = taumax(i);
      statev(316 + i) = taumin(i);
      taurange(i) = fem::abs(taumax(i)) + fem::abs(taumin(i));
      statev(292 + i) = taurange(i);
    }
  }
  //C
  fem::real_star_8 taurmax = fem::max(taurange(1), taurange(2),
    taurange(3), taurange(4), taurange(5), taurange(6), taurange(7),
    taurange(8), taurange(9), taurange(10), taurange(11), taurange(
    12));
  //C
  taurmax = taurmax / 2.0f;
  //C
  fem::real_star_8 d0 = 12.0e-6f;
  fem::real_star_8 a5 = cts(7) * cts(8) * cts(9);
  d = a5 / taurmax;
  d = fem::min(d0, d);
  //C
  //C   OVERLOAD EFFECT
  fem::real_star_8 ethai = cts(36);
  fem::real_star_8 ethaii = cts(37);
  if (kstep >= cts(46) && cts(46) != 0.0f) {
    if (props(4) == 1.0f) {
      d = fem::min(d, props(5));
    }
    else if (props(4) < (ethai + 0.001f)) {
      d = fem::min(d, props(5));
    }
    else if (props(4) < (ethaii + 0.001f)) {
      d = fem::min(d, props(5));
    }
  }
  //C
  statev(100) = d;
}

void
dshrsdgd(
  arr_ref<fem::real_star_8, 2> dtau,
  arr_cref<fem::real_star_8, 4> c0,
  arr_cref<fem::real_star_8, 2> slpdir0,
  arr_cref<fem::real_star_8, 2> slpnor0,
  fem::real_star_8 const& dtincr)
{
  dtau(dimension(12, 12));
  c0(dimension(3, 3, 3, 3));
  slpdir0(dimension(3, 12));
  slpnor0(dimension(3, 12));
  //C
  int i = fem::int0;
  int k = fem::int0;
  int l = fem::int0;
  int m = fem::int0;
  arr_2d<3, 3, fem::real_star_8> a1(fem::fill0);
  int ii = fem::int0;
  int kk = fem::int0;
  arr_2d<3, 3, fem::real_star_8> a2(fem::fill0);
  int ll = fem::int0;
  int mm = fem::int0;
  fem::real_star_8 sum = fem::zero<fem::real_star_8>();
  FEM_DO_SAFE(i, 1, 12) {
    FEM_DO_SAFE(k, 1, 12) {
      dtau(i, k) = 0;
      FEM_DO_SAFE(l, 1, 3) {
        FEM_DO_SAFE(m, 1, 3) {
          a1(l, m) = slpdir0(l, k) * slpnor0(m, k);
        }
      }
      //C
      FEM_DO_SAFE(ii, 1, 3) {
        FEM_DO_SAFE(kk, 1, 3) {
          a2(ii, kk) = 0.0f;
          FEM_DO_SAFE(ll, 1, 3) {
            FEM_DO_SAFE(mm, 1, 3) {
              a2(ii, kk) += c0(ii, kk, ll, mm) * a1(ll, mm);
            }
          }
        }
      }
      //C
      FEM_DO_SAFE(ii, 1, 3) {
        FEM_DO_SAFE(kk, 1, 3) {
          a1(ii, kk) = slpdir0(ii, i) * slpnor0(kk, i);
        }
      }
      //C
      sum = 0;
      FEM_DO_SAFE(ii, 1, 3) {
        FEM_DO_SAFE(kk, 1, 3) {
          sum += a1(ii, kk) * a2(ii, kk);
        }
      }
      //C
      dtau(i, k) = -1.0f * sum * dtincr;
    }
  }
  //C
}

void
inigd(
  arr_ref<fem::real_star_8> gdot,
  arr_cref<fem::real_star_8> gdot0,
  arr_cref<fem::real_star_8> tau,
  arr_cref<fem::real_star_8> b,
  arr_cref<fem::real_star_8> ts,
  arr_cref<fem::real_star_8> cts,
  int& l)
{
  gdot(dimension(12));
  gdot0(dimension(12));
  tau(dimension(12));
  b(dimension(12));
  ts(dimension(12));
  cts(dimension(46));
  //C
  fem::real_star_8 a1 = -cts(21) / (cts(22) * cts(23));
  fem::real_star_8 a2 = cts(26) * (cts(8) / cts(27));
  l = 0.0f;
  int i = fem::int0;
  arr_1d<12, fem::real_star_8> taueff(fem::fill0);
  arr_1d<12, fem::real_star_8> gau(fem::fill0);
  FEM_DO_SAFE(i, 1, 12) {
    taueff(i) = fem::abs(tau(i) - b(i)) - ts(i);
    if (taueff(i) <= 0.0f) {
      gdot(i) = gdot0(i) * fem::exp(a1) * fem::sign(1.e0, tau(i) - b(i));
    }
    else {
      gau(i) = 1.0f - fem::pow(((taueff(i)) / a2), cts(24));
      if (gau(i) <= 0.0f) {
        gdot(i) = gdot0(i) * fem::sign(1.e0, tau(i) - b(i));
        l = 1.0f;
      }
      else {
        gdot(i) = gdot0(i) * fem::exp(a1 * (fem::pow((gau(i)), cts(
          25)))) * fem::sign(1.e0, tau(i) - b(i));
      }
    }
  }
  //C
}

void
hill(
  arr_ref<fem::real_star_8> fshill,
  arr_cref<fem::real_star_8> tau,
  arr_cref<fem::real_star_8> tau0,
  arr_ref<fem::real_star_8> statev,
  arr_cref<fem::real_star_8> gdot,
  fem::real_star_8 const& dtime,
  fem::real_star_8 const& /* d */,
  int const& /* nincrt */)
{
  fshill(dimension(12));
  tau(dimension(12));
  tau0(dimension(12));
  statev(dimension(400));
  gdot(dimension(12));
  //C
  int i = fem::int0;
  arr_1d<12, fem::real_star_8> dgamma(fem::fill0);
  FEM_DO_SAFE(i, 1, 12) {
    dgamma(i) = dtime * gdot(i);
  }
  //C
  FEM_DO_SAFE(i, 1, 12) {
    //C      FSHILL(I)=(tanh(0.005*0.002/D))*0.5*abs(DGAMMA(I))/
    //C     1 (abs(TAU(I)-TAU0(I))/real(NINCRT))
    if (tau(i) - tau0(i) != 0.0f) {
      //C      FSHILL(I)=0.5*ABS( DGAMMA(I)/(TAU(I)-TAU0(I)) )
      //C          FSHILL(I)=(tanh(0.005*0.002/D))*0.5*abs(DGAMMA(I))/
      //C     1 (abs(TAU(I)-TAU0(I))/REAL(NINCRT))
      //C
      fshill(i) = 0.5f * (fem::abs(dgamma(i)) / fem::abs((fem::abs(
        tau(i)) - fem::abs(tau0(i)))));
    }
    else {
      fshill(i) = 0.0f;
    }
  }
  FEM_DO_SAFE(i, 1, 12) {
    statev(i + 328) = fshill(i);
  }
  //C
}

void
ajacobian(
  arr_ref<fem::real_star_8, 2> djac,
  fem::real_star_8 const& d,
  fem::real_star_8 const& etha,
  arr_cref<fem::real_star_8> ro,
  arr_cref<fem::real_star_8> gdot,
  arr_cref<fem::real_star_8> gdot0,
  arr_cref<fem::real_star_8> tau,
  arr_cref<fem::real_star_8, 2> dtau,
  arr_cref<fem::real_star_8> b,
  arr_cref<fem::real_star_8> ts,
  arr_cref<fem::real_star_8> statev,
  arr_cref<fem::real_star_8> func,
  arr_ref<fem::real_star_8> grad,
  arr_cref<fem::real_star_8> cts,
  fem::real_star_8 const& dtincr,
  fem::real_star_8 const& fw,
  arr_cref<fem::real_star_8> ro0)
{
  djac(dimension(12, 12));
  ro(dimension(12));
  gdot(dimension(12));
  gdot0(dimension(12));
  tau(dimension(12));
  dtau(dimension(12, 12));
  b(dimension(12));
  ts(dimension(12));
  statev(dimension(400));
  func(dimension(12));
  grad(dimension(12));
  cts(dimension(48));
  ro0(dimension(12));
  //C
  int i = fem::int0;
  int k = fem::int0;
  arr_2d<12, 12, fem::real_star_8> db(fem::fill0);
  arr_2d<12, 12, fem::real_star_8> ds(fem::fill0);
  FEM_DO_SAFE(i, 1, 12) {
    FEM_DO_SAFE(k, 1, 12) {
      db(i, k) = 0.0f;
      ds(i, k) = 0.0f;
      djac(i, k) = 0.0f;
    }
  }
  //C
  //C   DERIVITIVE OF BACK STRESS WRT GDOT
  //C
  fem::real_star_8 a0 = fem::zero<fem::real_star_8>();
  arr_1d<12, fem::real_star_8> a01(fem::fill0);
  FEM_DO_SAFE(i, 1, 12) {
    a0 = 1.0f + 4.0f * cts(16) * cts(8) * statev(i + 73);
    a01(i) = 2.0f * cts(8) * (1.0f - 2.0f * cts(16)) / a0;
    //C
    db(i, i) = dtincr * (fw / (1.0f - fw)) * a01(i);
  }
  //C
  //C       DERIVITIVE OF S WRT GDOT
  //C
  fem::real_star_8 a20 = cts(8) * cts(9) * cts(14);
  fem::real_star_8 term1 = fem::zero<fem::real_star_8>();
  fem::real_star_8 term2 = fem::zero<fem::real_star_8>();
  FEM_DO_SAFE(i, 1, 12) {
    //C
    term1 = (cts(10) / (etha * cts(9) * d)) * fem::sign(1.e0, gdot(i));
    term2 = (2.0f * cts(11) / cts(9)) * ro0(i) * fem::sign(1.e0, gdot(i));
    //C
    ds(i, i) = dtincr * (term1 - term2) * a20 / fem::sqrt(cts(14) * ro(i));
    //C
  }
  //C
  fem::real_star_8 a7 = -cts(21) / (cts(22) * cts(23));
  fem::real_star_8 a8 = cts(26) * (cts(8) / cts(27));
  fem::real_star_8 a10 = -a7 * cts(24) * cts(25) / a8;
  arr_1d<12, fem::real_star_8> taueff(fem::fill0);
  fem::real_star_8 term3 = fem::zero<fem::real_star_8>();
  fem::real_star_8 term4 = fem::zero<fem::real_star_8>();
  fem::real_star_8 term5 = fem::zero<fem::real_star_8>();
  arr_1d<12, fem::real_star_8> gau(fem::fill0);
  fem::real_star_8 term = fem::zero<fem::real_star_8>();
  FEM_DO_SAFE(i, 1, 12) {
    taueff(i) = fem::abs(tau(i) - b(i)) - ts(i);
    term1 = 0.00f;
    term2 = 0.00f;
    term3 = 0.00f;
    term4 = 0.00f;
    term5 = 0.00f;
    //C
    djac(i, i) = 1.0f;
    //C
    if (taueff(i) <= 0.0f) {
      djac(i, i) = 1.0f;
    }
    else {
      gau(i) = 1.0f - fem::pow((taueff(i) / a8), cts(24));
      if (gau(i) < 0.0f) {
        djac(i, i) = 1.0f;
      }
      else {
        //C
        FEM_DO_SAFE(k, 1, 12) {
          term1 = gdot0(i) * fem::exp(a7 * fem::pow((gau(i)), cts(25)));
          term2 = fem::pow(gau(i), (cts(25) - 1.0f));
          term3 = fem::pow((taueff(i) / a8), (cts(24) - 1.0f));
          term = a10 * term1 * term2 * term3;
          //C
          term4 = dtau(i, k) - db(i, k);
          term5 = ds(i, k);
          //C
          djac(i, k) = djac(i, k) - term * (term4 * fem::sign(1.e0,
            tau(i) - b(i)) - term5) * fem::sign(1.e0, tau(i) - b(i));
        }
      }
    }
  }
  //C
  int j = fem::int0;
  FEM_DO_SAFE(j, 1, 12) {
    grad(j) = 0.0f;
    FEM_DO_SAFE(i, 1, 12) {
      grad(j) = grad(j) - func(i) * djac(i, j);
    }
  }
  grad(1) = 2.0f * grad(1);
  grad(2) = 2.0f * grad(2);
  grad(3) = 2.0f * grad(3);
  //C
}

void
njacobian(
  arr_ref<fem::real_star_8, 2> djac,
  arr_cref<fem::real_star_8> func,
  arr_cref<fem::real_star_8> gdot,
  arr_ref<fem::real_star_8> pgdot,
  arr_ref<fem::real_star_8> dg0,
  arr_ref<fem::real_star_8> func0,
  arr_ref<fem::real_star_8> grad)
{
  djac(dimension(12, 12));
  func(dimension(12));
  gdot(dimension(12));
  pgdot(dimension(12));
  dg0(dimension(12));
  func0(dimension(12));
  grad(dimension(12));
  //C
  int i = fem::int0;
  arr_1d<12, fem::real_star_8> dg(fem::fill0);
  arr_1d<12, fem::real_star_8> ddg(fem::fill0);
  arr_1d<12, fem::real_star_8> dfunc(fem::fill0);
  FEM_DO_SAFE(i, 1, 12) {
    dg(i) = gdot(i) - pgdot(i);
    ddg(i) = dg(i) - dg0(i);
    dfunc(i) = func(i) - func0(i);
  }
  //C!  CHECK FOR IF INITIAL CALL
  fem::real_star_8 a4 = 0;
  int k = fem::int0;
  FEM_DO_SAFE(i, 1, 12) {
    FEM_DO_SAFE(k, 1, 12) {
      a4 += fem::abs(djac(i, k));
    }
  }
  //C!  INITIAL J
  fem::real_star_8 gdotgt = fem::zero<fem::real_star_8>();
  arr_1d<12, fem::real_star_8> a2(fem::fill0);
  arr_2d<12, 12, fem::real_star_8> a3(fem::fill0);
  arr_1d<12, fem::real_star_8> a1(fem::fill0);
  if (a4 == 0.0f) {
    FEM_DO_SAFE(i, 1, 12) {
      djac(i, i) = 1.0f;
    }
  }
  else {
    //C
    gdotgt = 0.0f;
    FEM_DO_SAFE(i, 1, 12) {
      gdotgt += fem::pow(ddg(i), 2.0f);
    }
    a2 = 0.0f;
    a3 = 0.0f;
    FEM_DO_SAFE(i, 1, 12) {
      a1(i) = 0.0f;
      FEM_DO_SAFE(k, 1, 12) {
        a1(i) += djac(i, k) * ddg(k);
      }
    }
    FEM_DO_SAFE(i, 1, 12) {
      a2(i) = dfunc(i) - a1(i);
    }
    FEM_DO_SAFE(i, 1, 12) {
      FEM_DO_SAFE(k, 1, 12) {
        a3(i, k) = a2(i) * ddg(k);
      }
    }
    FEM_DO_SAFE(i, 1, 12) {
      FEM_DO_SAFE(k, 1, 12) {
        djac(i, k) += (1 / gdotgt) * a3(i, k);
      }
    }
    //C
  }
  //C
  int j = fem::int0;
  FEM_DO_SAFE(j, 1, 12) {
    grad(j) = 0.0f;
    FEM_DO_SAFE(i, 1, 12) {
      grad(j) = grad(j) - func(i) * djac(i, j);
    }
  }
  grad(1) = 2.0f * grad(1);
  grad(2) = 2.0f * grad(2);
  grad(3) = 2.0f * grad(3);
  //C
  FEM_DO_SAFE(i, 1, 12) {
    pgdot(i) = gdot(i);
    dg0(i) = dg(i);
    func0(i) = func(i);
  }
  //C
}

void
stressupdate(
  arr_ref<fem::real_star_8> stress,
  arr_cref<fem::real_star_8, 2> s,
  int const& ndi,
  int const& nshr)
{
  stress(dimension(6));
  s(dimension(3, 3));
  //C
  int i = fem::int0;
  FEM_DO_SAFE(i, 1, ndi) {
    stress(i) = s(i, i);
  }
  if (nshr == 1.0f) {
    stress(ndi + 1) = s(1, 2);
  }
  if (nshr == 3.0f) {
    stress(4) = s(1, 2);
    stress(5) = s(1, 3);
    stress(6) = s(2, 3);
  }
  //C
}

void
ludcmp(
  int const& n,
  arr_ref<fem::real_star_8, 2> a,
  arr_ref<fem::real_star_8> indx)
{
  a(dimension(n, n));
  indx(dimension(n));
  //C
  fem::real_star_8 tiny = 1.0e-20f;
  //C
  int i = fem::int0;
  fem::real_star_8 aamax = fem::zero<fem::real_star_8>();
  int j = fem::int0;
  arr<fem::real_star_8> vv(dimension(n), fem::fill0);
  FEM_DO_SAFE(i, 1, n) {
    aamax = 0.0f;
    //C
    FEM_DO_SAFE(j, 1, n) {
      if (fem::abs(a(i, j)) > aamax) {
        aamax = fem::abs(a(i, j));
      }
    }
    //C
    vv(i) = 1.0f / aamax;
  }
  //C
  fem::real_star_8 sum = fem::zero<fem::real_star_8>();
  int k = fem::int0;
  fem::real_star_8 dum = fem::zero<fem::real_star_8>();
  int imax = fem::int0;
  FEM_DO_SAFE(j, 1, n) {
    FEM_DO_SAFE(i, 1, j - 1) {
      sum = a(i, j);
      FEM_DO_SAFE(k, 1, i - 1) {
        sum = sum - a(i, k) * a(k, j);
      }
    }
    aamax = 0.0f;
    FEM_DO_SAFE(i, j, n) {
      sum = a(i, j);
      FEM_DO_SAFE(k, 1, j - 1) {
        sum = sum - a(i, k) * a(k, j);
      }
      a(i, j) = sum;
      dum = vv(i) * fem::abs(sum);
      if (dum >= aamax) {
        imax = i;
        aamax = dum;
      }
    }
    //C
    if (j != imax) {
      FEM_DO_SAFE(k, 1, n) {
        dum = a(imax, k);
        a(imax, k) = a(j, k);
        a(j, k) = dum;
      }
      vv(imax) = vv(j);
    }
    indx(j) = imax;
    if (a(j, j) == 0.f) {
      a(j, j) = tiny;
    }
    if (j != n) {
      dum = 1.f / a(j, j);
      FEM_DO_SAFE(i, j + 1, n) {
        a(i, j) = a(i, j) * dum;
      }
    }
  }
  //C
}

//C
//C----------------------------------------------------------------------
//C
void
lubksb(
  arr_cref<fem::real_star_8, 2> a,
  int const& n,
  arr_cref<fem::real_star_8> indx,
  arr_ref<fem::real_star_8> b)
{
  a(dimension(n, n));
  indx(dimension(n));
  b(dimension(n));
  //C
  int ii = 0;
  int i = fem::int0;
  int ll = fem::int0;
  fem::real_star_8 sum = fem::zero<fem::real_star_8>();
  int j = fem::int0;
  FEM_DO_SAFE(i, 1, n) {
    ll = indx(i);
    sum = b(ll);
    b(ll) = b(i);
    //C
    if (ii != 0) {
      FEM_DO_SAFE(j, ii, i - 1) {
        sum = sum - a(i, j) * b(j);
      }
    }
    else if (sum != 0.f) {
      ii = i;
    }
    //C
    b(i) = sum;
  }
  //C
  FEM_DOSTEP(i, n, 1, -1) {
    sum = b(i);
    //C
    if (i < n) {
      FEM_DO_SAFE(j, i + 1, n) {
        sum = sum - a(i, j) * b(j);
      }
    }
    //C
    b(i) = sum / a(i, i);
  }
  //C
}

void
stff(
  arr_ref<fem::real_star_8, 2> ddsdde,
  arr_cref<fem::real_star_8, 4> c,
  arr_cref<fem::real_star_8> gdot,
  arr_cref<fem::real_star_8> tau,
  arr_cref<fem::real_star_8> b,
  arr_cref<fem::real_star_8> ts,
  arr_cref<fem::real_star_8, 2> slpdir0,
  arr_cref<fem::real_star_8, 2> slpnor0,
  arr_cref<fem::real_star_8, 2> fel,
  arr_cref<fem::real_star_8, 2> felinv,
  arr_cref<fem::real_star_8> cts,
  fem::real_star_8 const& dtime)
{
  ddsdde(dimension(6, 6));
  c(dimension(3, 3, 3, 3));
  gdot(dimension(12));
  tau(dimension(12));
  b(dimension(12));
  ts(dimension(12));
  slpdir0(dimension(3, 12));
  slpnor0(dimension(3, 12));
  fel(dimension(3, 3));
  felinv(dimension(3, 3));
  cts(dimension(48));
  //C
  int n = fem::int0;
  int i = fem::int0;
  arr_2d<3, 12, fem::real_star_8> slpdir1(fem::fill0);
  arr_2d<3, 12, fem::real_star_8> slpnor1(fem::fill0);
  int k = fem::int0;
  FEM_DO_SAFE(n, 1, 12) {
    FEM_DO_SAFE(i, 1, 3) {
      slpdir1(i, n) = 0.0f;
      slpnor1(i, n) = 0.0f;
      FEM_DO_SAFE(k, 1, 3) {
        slpdir1(i, n) += fel(i, k) * slpdir0(k, n);
        slpnor1(i, n) += slpnor0(k, n) * felinv(k, i);
      }
    }
  }
  //C
  fem::real_star_8 a8 = cts(26) * (cts(8) / cts(27));
  fem::real_star_8 a9 = -cts(21) / (cts(22) * cts(23));
  fem::real_star_8 a10 = -a9 * cts(24) * cts(25) / a8;
  //C
  int j = fem::int0;
  int l = fem::int0;
  arr<fem::real_star_8, 4> ddpds(dimension(3, 3, 3, 3), fem::fill0);
  int m = fem::int0;
  arr_1d<12, fem::real_star_8> taueff(fem::fill0);
  arr_1d<12, fem::real_star_8> gau(fem::fill0);
  FEM_DO_SAFE(i, 1, 3) {
    FEM_DO_SAFE(j, 1, 3) {
      FEM_DO_SAFE(k, 1, 3) {
        FEM_DO_SAFE(l, 1, 3) {
          ddpds(i, j, k, l) = 0.0f;
          FEM_DO_SAFE(m, 1, 12) {
            //C
            taueff(m) = fem::abs(tau(m) - b(m)) - ts(m);
            //C
            if (taueff(m) > 0.0f) {
              gau(m) = 1.0f - fem::pow(((taueff(m)) / a8), cts(24));
            }
            else {
              gau(m) = 1.0f;
              taueff(m) = 0.0f;
            }
            //C
            if (gau(m) < 0.0f) {
              gau(m) = 0.0f;
            }
            if (taueff(m) > 0.0f) {
              ddpds(i, j, k, l) += (slpdir1(i, m) * slpnor1(j, m) + slpnor1(i,
                m) * slpdir1(j, m)) * (slpdir1(k, m) * slpnor1(l,
                m) + slpnor1(k, m) * slpdir1(l, m)) * (gdot(m)) * (
                fem::pow(gau(m), (cts(25) - 1.0f))) * (fem::pow((
                taueff(m) / a8), (cts(24) - 1.0f))) * a10;
            }
          }
        }
      }
    }
  }
  //C
  FEM_DO_SAFE(i, 1, 3) {
    FEM_DO_SAFE(j, 1, 3) {
      FEM_DO_SAFE(k, 1, 3) {
        FEM_DO_SAFE(l, 1, 3) {
          ddpds(i, j, k, l) = ddpds(i, j, k, l) * dtime / 4.0f;
          ddpds(i, j, k, l) = 0.0f;
        }
      }
    }
  }
  //C
  arr<fem::real_star_8, 4> e(dimension(3, 3, 3, 3), fem::fill0);
  int m1 = fem::int0;
  int m2 = fem::int0;
  FEM_DO_SAFE(i, 1, 3) {
    FEM_DO_SAFE(j, 1, 3) {
      FEM_DO_SAFE(k, 1, 3) {
        FEM_DO_SAFE(l, 1, 3) {
          e(i, j, k, l) = 0.0f;
          FEM_DO_SAFE(m1, 1, 3) {
            FEM_DO_SAFE(m2, 1, 3) {
              e(i, j, k, l) += c(i, j, m1, m2) * ddpds(m1, m2, k, l);
            }
          }
        }
      }
    }
  }
  //C
  arr_2d<3, 3, fem::real_star_8> del(fem::fill0);
  FEM_DO_SAFE(i, 1, 3) {
    FEM_DO_SAFE(k, 1, 3) {
      del(i, k) = 0.0f;
    }
  }
  del(1, 1) = 1.0f;
  del(2, 2) = 1.0f;
  del(3, 3) = 1.0f;
  //C
  FEM_DO_SAFE(i, 1, 3) {
    FEM_DO_SAFE(j, 1, 3) {
      FEM_DO_SAFE(k, 1, 3) {
        FEM_DO_SAFE(l, 1, 3) {
          e(i, j, k, l) += 0.5f * (del(i, k) * del(j, l) + del(i,
            l) * del(j, k));
        }
      }
    }
  }
  //C
  int ia = fem::int0;
  int ib = fem::int0;
  arr_2d<6, 6, fem::real_star_8> a1(fem::fill0);
  FEM_DO_SAFE(i, 1, 3) {
    FEM_DO_SAFE(j, i, 3) {
      ia = i;
      if (i != j) {
        ia = 9 - i - j;
      }
      FEM_DO_SAFE(k, 1, 3) {
        FEM_DO_SAFE(l, k, 3) {
          ib = k;
          if (k != l) {
            ib = 9 - k - l;
          }
          a1(ia, ib) = e(i, j, k, l);
        }
      }
    }
  }
  //C
  arr_2d<6, 6, fem::real_star_8> a0(fem::fill0);
  FEM_DO_SAFE(i, 1, 6) {
    FEM_DO_SAFE(k, 1, 6) {
      a0(i, k) = a1(i, k);
    }
  }
  //C
  arr_2d<6, 6, fem::real_star_8> a2(fem::fill0);
  FEM_DO_SAFE(i, 1, 6) {
    FEM_DO_SAFE(k, 1, 6) {
      a2(i, k) = 0.0f;
    }
    a2(i, i) = 1.0f;
  }
  //C
  arr_1d<6, fem::real_star_8> indx(fem::fill0);
  ludcmp(6, a0, indx);
  FEM_DO_SAFE(i, 1, 6) {
    lubksb(a0, 6, indx, a2(1, i));
  }
  //C
  arr<fem::real_star_8, 4> a3(dimension(3, 3, 3, 3), fem::fill0);
  FEM_DO_SAFE(i, 1, 3) {
    FEM_DO_SAFE(j, 1, 3) {
      ia = i;
      if (i != j) {
        ia = 9 - i - j;
      }
      FEM_DO_SAFE(k, 1, 3) {
        FEM_DO_SAFE(l, 1, 3) {
          ib = k;
          if (k != l) {
            ib = 9 - k - l;
          }
          a3(i, j, k, l) = a2(ia, ib);
          if (ia > 3) {
            a3(i, j, k, l) = a3(i, j, k, l) / 2.0f;
          }
          if (ib > 3) {
            a3(i, j, k, l) = a3(i, j, k, l) / 2.0f;
          }
        }
      }
    }
  }
  //C
  arr<fem::real_star_8, 4> a4(dimension(3, 3, 3, 3), fem::fill0);
  FEM_DO_SAFE(i, 1, 3) {
    FEM_DO_SAFE(j, 1, 3) {
      FEM_DO_SAFE(k, 1, 3) {
        FEM_DO_SAFE(l, 1, 3) {
          a4(i, j, k, l) = 0.0f;
          FEM_DO_SAFE(m1, 1, 3) {
            FEM_DO_SAFE(m2, 1, 3) {
              a4(i, j, k, l) += a3(i, j, m1, m2) * c(m1, m2, k, l);
            }
          }
        }
      }
    }
  }
  //C
  arr_2d<6, 6, fem::real_star_8> a5(fem::fill0);
  FEM_DO_SAFE(i, 1, 3) {
    FEM_DO_SAFE(j, i, 3) {
      ia = i;
      if (i != j) {
        ia = i + j + 1;
      }
      FEM_DO_SAFE(k, 1, 3) {
        FEM_DO_SAFE(l, k, 3) {
          ib = k;
          if (k != l) {
            ib = k + l + 1;
          }
          a5(ia, ib) = a4(i, j, k, l);
          if (ib >= 4) {
            a5(ia, ib) = a5(ia, ib) * 2.0f;
          }
        }
      }
    }
  }
  //C
  FEM_DO_SAFE(i, 1, 6) {
    FEM_DO_SAFE(j, 1, 6) {
      ddsdde(i, j) = 0.0f;
    }
  }
  //C       DO I=1,6
  //C          DO J=1,6
  //C               DDSDDE(I,J)=A5(I,J)
  //C          END DO
  //C      END DO
  FEM_DO_SAFE(i, 1, 3) {
    FEM_DO_SAFE(j, i, 3) {
      ia = i;
      if (i != j) {
        ia = i + j + 1;
      }
      FEM_DO_SAFE(k, 1, 3) {
        FEM_DO_SAFE(l, k, 3) {
          ib = k;
          if (k != l) {
            ib = k + l + 1;
          }
          ddsdde(ia, ib) = c(i, j, k, l);
        }
      }
    }
  }
  //C
}

void
endstore(
  arr_cref<fem::real_star_8> gdot,
  arr_ref<fem::real_star_8> statev,
  arr_cref<fem::real_star_8, 2> fpinv,
  arr_cref<fem::real_star_8> ro,
  arr_cref<fem::real_star_8> b,
  arr_cref<fem::real_star_8> tau,
  arr_cref<fem::real_star_8> ts,
  arr_cref<fem::real_star_8, 2> slpdir0,
  arr_cref<fem::real_star_8, 2> slpnor0,
  fem::real_star_8 const& dtincr,
  int const& kstep,
  arr_cref<fem::real_star_8, 2> rotate)
{
  gdot(dimension(12));
  statev(dimension(400));
  fpinv(dimension(3, 3));
  ro(dimension(12));
  b(dimension(12));
  tau(dimension(12));
  ts(dimension(12));
  slpdir0(dimension(3, 12));
  slpnor0(dimension(3, 12));
  rotate(dimension(3, 3));
  //C
  int i = fem::int0;
  arr_1d<12, fem::real_star_8> dgamma(fem::fill0);
  arr_1d<12, fem::real_star_8> eslip(fem::fill0);
  FEM_DO_SAFE(i, 1, 12) {
    dgamma(i) = dtincr * gdot(i);
    statev(i + 20) += dgamma(i);
    eslip(i) = statev(i + 20);
  }
  //C
  arr_1d<12, fem::real_star_8> erange(fem::fill0);
  if (kstep == 1.0f) {
    FEM_DO_SAFE(i, 1, 12) {
      statev(i + 100) = eslip(i);
      erange(i) = eslip(i);
    }
  }
  else if (kstep == 2.0f) {
    FEM_DO_SAFE(i, 1, 12) {
      statev(i + 112) = eslip(i);
      erange(i) = statev(i + 100);
    }
  }
  else if (kstep == 3.0f) {
    FEM_DO_SAFE(i, 1, 12) {
      statev(i + 124) = eslip(i);
      erange(i) = (statev(i + 100) - statev(i + 112));
    }
  }
  else if (kstep == 4.0f) {
    FEM_DO_SAFE(i, 1, 12) {
      statev(i + 100) = eslip(i);
      erange(i) = (statev(i + 112) - statev(i + 124));
    }
  }
  else if (kstep == 5.0f) {
    FEM_DO_SAFE(i, 1, 12) {
      statev(i + 112) = eslip(i);
      erange(i) = (statev(i + 124) - statev(i + 100));
    }
  }
  else if (kstep == 6.0f) {
    FEM_DO_SAFE(i, 1, 12) {
      statev(i + 124) = eslip(i);
      erange(i) = (statev(i + 100) - statev(i + 112));
    }
  }
  else if (kstep == 7.0f) {
    FEM_DO_SAFE(i, 1, 12) {
      statev(i + 100) = eslip(i);
      erange(i) = (statev(i + 112) - statev(i + 124));
    }
  }
  else if (kstep == 8.0f) {
    FEM_DO_SAFE(i, 1, 12) {
      statev(i + 112) = eslip(i);
      erange(i) = (statev(i + 124) - statev(i + 100));
    }
  }
  else if (kstep == 9.0f) {
    FEM_DO_SAFE(i, 1, 12) {
      statev(i + 124) = eslip(i);
      erange(i) = (statev(i + 100) - statev(i + 112));
    }
  }
  else if (kstep == 10.0f) {
    FEM_DO_SAFE(i, 1, 12) {
      statev(i + 100) = eslip(i);
      erange(i) = (statev(i + 112) - statev(i + 124));
    }
  }
  else if (kstep == 11.0f) {
    FEM_DO_SAFE(i, 1, 12) {
      statev(i + 112) = eslip(i);
      erange(i) = (statev(i + 124) - statev(i + 100));
    }
  }
  else if (kstep == 12.0f) {
    FEM_DO_SAFE(i, 1, 12) {
      statev(i + 124) = eslip(i);
      erange(i) = (statev(i + 100) - statev(i + 112));
    }
  }
  else if (kstep == 13.0f) {
    FEM_DO_SAFE(i, 1, 12) {
      statev(i + 100) = eslip(i);
      erange(i) = (statev(i + 112) - statev(i + 124));
    }
  }
  else if (kstep == 14.0f) {
    FEM_DO_SAFE(i, 1, 12) {
      statev(i + 112) = eslip(i);
      erange(i) = (statev(i + 124) - statev(i + 100));
    }
  }
  else if (kstep == 15.0f) {
    FEM_DO_SAFE(i, 1, 12) {
      statev(i + 124) = eslip(i);
      erange(i) = (statev(i + 100) - statev(i + 112));
    }
  }
  else if (kstep == 16.0f) {
    FEM_DO_SAFE(i, 1, 12) {
      statev(i + 100) = eslip(i);
      erange(i) = (statev(i + 112) - statev(i + 124));
    }
  }
  else if (kstep == 17.0f) {
    FEM_DO_SAFE(i, 1, 12) {
      statev(i + 112) = eslip(i);
      erange(i) = (statev(i + 124) - statev(i + 100));
    }
  }
  else if (kstep == 18.0f) {
    FEM_DO_SAFE(i, 1, 12) {
      statev(i + 124) = eslip(i);
      erange(i) = (statev(i + 100) - statev(i + 112));
    }
  }
  else if (kstep == 19.0f) {
    FEM_DO_SAFE(i, 1, 12) {
      statev(i + 100) = eslip(i);
      erange(i) = (statev(i + 112) - statev(i + 124));
    }
  }
  else if (kstep == 20.0f) {
    FEM_DO_SAFE(i, 1, 12) {
      statev(i + 112) = eslip(i);
      erange(i) = (statev(i + 124) - statev(i + 100));
    }
  }
  else if (kstep == 21.0f) {
    FEM_DO_SAFE(i, 1, 12) {
      statev(i + 124) = eslip(i);
      erange(i) = (statev(i + 100) - statev(i + 112));
    }
  }
  else if (kstep == 22.0f) {
    FEM_DO_SAFE(i, 1, 12) {
      statev(i + 100) = eslip(i);
      erange(i) = (statev(i + 112) - statev(i + 124));
    }
  }
  else if (kstep == 23.0f) {
    FEM_DO_SAFE(i, 1, 12) {
      statev(i + 112) = eslip(i);
      erange(i) = (statev(i + 124) - statev(i + 100));
    }
  }
  else if (kstep == 24.0f) {
    FEM_DO_SAFE(i, 1, 12) {
      statev(i + 124) = eslip(i);
      erange(i) = (statev(i + 100) - statev(i + 112));
    }
  }
  else if (kstep == 25.0f) {
    FEM_DO_SAFE(i, 1, 12) {
      statev(i + 100) = eslip(i);
      erange(i) = (statev(i + 112) - statev(i + 124));
    }
  }
  else if (kstep == 26.0f) {
    FEM_DO_SAFE(i, 1, 12) {
      statev(i + 112) = eslip(i);
      erange(i) = (statev(i + 124) - statev(i + 100));
    }
  }
  else if (kstep == 27.0f) {
    FEM_DO_SAFE(i, 1, 12) {
      statev(i + 124) = eslip(i);
      erange(i) = (statev(i + 100) - statev(i + 112));
    }
  }
  else if (kstep == 28.0f) {
    FEM_DO_SAFE(i, 1, 12) {
      statev(i + 100) = eslip(i);
      erange(i) = (statev(i + 112) - statev(i + 124));
    }
  }
  else if (kstep == 29.0f) {
    FEM_DO_SAFE(i, 1, 12) {
      statev(i + 112) = eslip(i);
      erange(i) = (statev(i + 124) - statev(i + 100));
    }
  }
  else if (kstep == 30.0f) {
    FEM_DO_SAFE(i, 1, 12) {
      statev(i + 124) = eslip(i);
      erange(i) = (statev(i + 100) - statev(i + 112));
    }
  }
  else if (kstep == 31.0f) {
    FEM_DO_SAFE(i, 1, 12) {
      statev(i + 100) = eslip(i);
      erange(i) = (statev(i + 112) - statev(i + 124));
    }
  }
  else if (kstep == 32.0f) {
    FEM_DO_SAFE(i, 1, 12) {
      statev(i + 112) = eslip(i);
      erange(i) = (statev(i + 124) - statev(i + 100));
    }
  }
  else if (kstep == 33.0f) {
    FEM_DO_SAFE(i, 1, 12) {
      statev(i + 124) = eslip(i);
      erange(i) = (statev(i + 100) - statev(i + 112));
    }
  }
  else if (kstep == 34.0f) {
    FEM_DO_SAFE(i, 1, 12) {
      statev(i + 100) = eslip(i);
      erange(i) = (statev(i + 112) - statev(i + 124));
    }
  }
  else if (kstep == 35.0f) {
    FEM_DO_SAFE(i, 1, 12) {
      statev(i + 112) = eslip(i);
      erange(i) = (statev(i + 124) - statev(i + 100));
    }
  }
  else if (kstep == 36.0f) {
    FEM_DO_SAFE(i, 1, 12) {
      statev(i + 124) = eslip(i);
      erange(i) = (statev(i + 100) - statev(i + 112));
    }
  }
  else if (kstep == 37.0f) {
    FEM_DO_SAFE(i, 1, 12) {
      statev(i + 100) = eslip(i);
      erange(i) = (statev(i + 112) - statev(i + 124));
    }
  }
  else if (kstep == 38.0f) {
    FEM_DO_SAFE(i, 1, 12) {
      statev(i + 112) = eslip(i);
      erange(i) = (statev(i + 124) - statev(i + 100));
    }
  }
  else if (kstep == 39.0f) {
    FEM_DO_SAFE(i, 1, 12) {
      statev(i + 124) = eslip(i);
      erange(i) = (statev(i + 100) - statev(i + 112));
    }
  }
  else if (kstep == 40.0f) {
    FEM_DO_SAFE(i, 1, 12) {
      statev(i + 100) = eslip(i);
      erange(i) = (statev(i + 112) - statev(i + 124));
    }
  }
  else if (kstep == 41.0f) {
    FEM_DO_SAFE(i, 1, 12) {
      statev(i + 112) = eslip(i);
      erange(i) = (statev(i + 124) - statev(i + 100));
    }
  }
  else if (kstep == 42.0f) {
    FEM_DO_SAFE(i, 1, 12) {
      statev(i + 124) = eslip(i);
      erange(i) = (statev(i + 100) - statev(i + 112));
    }
  }
  else if (kstep == 43.0f) {
    FEM_DO_SAFE(i, 1, 12) {
      statev(i + 100) = eslip(i);
      erange(i) = (statev(i + 112) - statev(i + 124));
    }
  }
  else if (kstep == 44.0f) {
    FEM_DO_SAFE(i, 1, 12) {
      statev(i + 112) = eslip(i);
      erange(i) = (statev(i + 124) - statev(i + 100));
    }
  }
  else if (kstep == 45.0f) {
    FEM_DO_SAFE(i, 1, 12) {
      statev(i + 124) = eslip(i);
      erange(i) = (statev(i + 100) - statev(i + 112));
    }
  }
  else if (kstep == 46.0f) {
    FEM_DO_SAFE(i, 1, 12) {
      statev(i + 100) = eslip(i);
      erange(i) = (statev(i + 112) - statev(i + 124));
    }
  }
  else if (kstep == 47.0f) {
    FEM_DO_SAFE(i, 1, 12) {
      statev(i + 112) = eslip(i);
      erange(i) = (statev(i + 124) - statev(i + 100));
    }
  }
  else if (kstep == 48.0f) {
    FEM_DO_SAFE(i, 1, 12) {
      statev(i + 124) = eslip(i);
      erange(i) = (statev(i + 100) - statev(i + 112));
    }
  }
  else if (kstep == 49.0f) {
    FEM_DO_SAFE(i, 1, 12) {
      statev(i + 100) = eslip(i);
      erange(i) = (statev(i + 112) - statev(i + 124));
    }
  }
  else if (kstep == 50.0f) {
    FEM_DO_SAFE(i, 1, 12) {
      statev(i + 112) = eslip(i);
      erange(i) = (statev(i + 124) - statev(i + 100));
    }
  }
  else if (kstep == 51.0f) {
    FEM_DO_SAFE(i, 1, 12) {
      statev(i + 124) = eslip(i);
      erange(i) = (statev(i + 100) - statev(i + 112));
    }
  }
  //C
  FEM_DO_SAFE(i, 1, 12) {
    erange(i) = fem::abs(erange(i));
    statev(i + 136) = erange(i);
  }
  //C
  int k = fem::int0;
  arr_2d<3, 3, fem::real_star_8> depl(fem::fill0);
  FEM_DO_SAFE(i, 1, 3) {
    FEM_DO_SAFE(k, 1, 3) {
      depl(i, k) = 0.0f;
    }
  }
  int n = fem::int0;
  FEM_DO_SAFE(i, 1, 3) {
    FEM_DO_SAFE(k, 1, 3) {
      FEM_DO_SAFE(n, 1, 12) {
        depl(i, k) += 0.5f * (slpdir0(i, n) * slpnor0(k, n) + slpdir0(k,
          n) * slpnor0(i, n)) * dgamma(n);
      }
    }
  }
  //C
  n = 9.0f;
  arr_2d<3, 3, fem::real_star_8> epl(fem::fill0);
  FEM_DO_SAFE(i, 1, 3) {
    FEM_DO_SAFE(k, 1, 3) {
      n++;
      epl(i, k) = statev(n) + depl(i, k);
    }
  }
  //C
  fem::real_star_8 sum = 0.0f;
  fem::real_star_8 sum1 = 0.0f;
  FEM_DO_SAFE(i, 1, 3) {
    FEM_DO_SAFE(k, 1, 3) {
      sum1 += fem::pow2(depl(i, k));
      sum += fem::pow2(epl(i, k));
    }
  }
  fem::real_star_8 eplacc = fem::sqrt((2.0f / 3.0f) * sum1);
  fem::real_star_8 epleff = fem::sqrt((2.0f / 3.0f) * sum);
  //C
  statev(19) = epleff;
  statev(20) = eplacc;
  //C
  n = 0.0f;
  int j = fem::int0;
  FEM_DO_SAFE(i, 1, 3) {
    FEM_DO_SAFE(j, 1, 3) {
      n++;
      statev(n) = fpinv(i, j);
    }
  }
  //C
  n = 9.0f;
  FEM_DO_SAFE(i, 1, 3) {
    FEM_DO_SAFE(k, 1, 3) {
      n++;
      statev(n) = epl(i, k);
    }
  }
  //C
  arr_1d<12, fem::real_star_8> gdot0(fem::fill0);
  FEM_DO_SAFE(i, 1, 12) {
    statev(i + 32) = ro(i);
    statev(i + 44) = b(i);
    statev(i + 56) = tau(i);
    statev(i + 68) = ts(i);
    statev(i + 388) = gdot0(i);
  }
  //C
  n = 340.0f;
  FEM_DO_SAFE(i, 1, 3) {
    FEM_DO_SAFE(k, 1, 3) {
      n++;
      statev(n) = rotate(i, k);
    }
  }
  //C
}

//C*********************************************************************
//C
//C   This UMAT corresponds to the work by Castelluccio and Co-workers for
//C   predicting FCC single crystals under monotonic loading.
//C
//C   Please cite as follows:
//C    "S Dindarlou and G M Castelluccio. Substructure-Sensitive
//C    Crystal Plasticity with Material-Invariant Parameters,
//C    International Journal of Plasticity 155 (2022): 103306"
//C
//C*********************************************************************
//C
//C  Constants are predetermined on a separate file read by the pre-compiler
//C
//C  Crystal orientation can be defined with Miller indexes or Euler angles.
//C  See Readme.txt file.
//C
//C*********************************************************************
//C
void
umat(
  arr_ref<fem::real_star_8> stress,
  arr_ref<fem::real_star_8> statev,
  arr_ref<fem::real_star_8, 2> ddsdde,
  fem::real_star_8& sse,
  fem::real_star_8 const& /* spd */,
  fem::real_star_8 const& /* scd */,
  fem::real_star_8 const& /* rpl */,
  arr_cref<fem::real_star_8> /* ddsddt */,
  arr_cref<fem::real_star_8> /* drplde */,
  fem::real_star_8 const& /* drpldt */,
  arr_cref<fem::real_star_8> /* stran */,
  arr_cref<fem::real_star_8> /* dstran */,
  arr_cref<fem::real_star_8> time,
  fem::real_star_8 const& dtime,
  fem::real_star_8 const& /* temp */,
  fem::real_star_8 const& /* dtemp */,
  arr_cref<fem::real_star_8> /* predef */,
  arr_cref<fem::real_star_8> /* dpred */,
  str_cref /* cmname */,
  int const& ndi,
  int const& nshr,
  int const& ntens,
  int const& nstatv,
  arr_cref<fem::real_star_8> props,
  int const& nprops,
  arr_cref<fem::real_star_8> /* coords */,
  arr_cref<fem::real_star_8, 2> /* drot */,
  fem::real_star_8& pnewdt,
  fem::real_star_8 const& /* celent */,
  arr_cref<fem::real_star_8, 2> dfgrd0,
  arr_cref<fem::real_star_8, 2> dfgrd1,
  int const& noel,
  int const& npt,
  fem::real_star_8 const& /* layer */,
  int const& /* kspt */,
  int const& kstep,
  int const& /* kinc */)
{
  stress(dimension(ntens));
  statev(dimension(nstatv));
  ddsdde(dimension(ntens, ntens));
  time(dimension(2));
  props(dimension(nprops));
  dfgrd0(dimension(3, 3));
  dfgrd1(dimension(3, 3));
  const fem::real_star_8 mi_ea = 1;
  arr_1d<48, fem::real_star_8> cts(fem::fill0);
  const fem::real_star_8 h1 = 1;
  arr_2d<3, 3, fem::real_star_8> rotate(fem::fill0);
  const fem::real_star_8 h2 = 0;
  const fem::real_star_8 h3 = 0;
  const fem::real_star_8 c11 = 249000000000.0f;
  const fem::real_star_8 c12 = 155000000000.0f;
  const fem::real_star_8 c44 = 114600000000.0f;
  const fem::real_star_8 skcy = 4.5f;
  const fem::real_star_8 shm = 80600000000.0f;
  const fem::real_star_8 bv = 2.56e-10f;
  const fem::real_star_8 dkm = 1;
  const fem::real_star_8 ysg = 2e-09f;
  const fem::real_star_8 ysc = 5e-08f;
  const fem::real_star_8 ale = 1;
  const fem::real_star_8 aii = 0.1f;
  const fem::real_star_8 gammacr = 0.03f;
  const fem::real_star_8 s1212 = 0.24f;
  const fem::real_star_8 fh = 0;
  const fem::real_star_8 gp = 0.002f;
  const fem::real_star_8 bf0 = 0.45f;
  const fem::real_star_8 bfi = 0.15f;
  const fem::real_star_8 f0 = 250000;
  const fem::real_star_8 bk = 8.314f;
  const fem::real_star_8 t = 298;
  const fem::real_star_8 pp = 0.67f;
  const fem::real_star_8 pq = 1.5f;
  const fem::real_star_8 s0t = 120000000;
  const fem::real_star_8 shm0 = 92600000000.0f;
  const fem::real_star_8 af = 1000000000000.0f;
  const fem::real_star_8 phi = 0.49f;
  const fem::real_star_8 vcs0 = 1.67e-26f;
  const fem::real_star_8 dd0 = 100000000.0f;
  const fem::real_star_8 ss = 0;
  const fem::real_star_8 fwi = 0.23f;
  const fem::real_star_8 fwii = 0.24f;
  const fem::real_star_8 fwiii = 0.175f;
  const fem::real_star_8 ethai = 3.7f;
  const fem::real_star_8 ethaii = 5.15f;
  const fem::real_star_8 dgm = 0.0001f;
  const fem::real_star_8 dgmpsb = 0.018f;
  const fem::real_star_8 etha_vein = 50;
  const fem::real_star_8 etha_psb = 15;
  const fem::real_star_8 etha_cell = 1;
  const fem::real_star_8 etha_lab = 5;
  const fem::real_star_8 load = 1;
  const fem::real_star_8 skmon = 4.5f;
  const fem::real_star_8 ovcyn = 0;
  const fem::real_star_8 crrot = 0;
  arr_1d<12, fem::real_star_8> b0(fem::fill0);
  arr_1d<12, fem::real_star_8> ro0(fem::fill0);
  arr_1d<12, fem::real_star_8> tau0(fem::fill0);
  arr_2d<3, 3, fem::real_star_8> fpinv0(fem::fill0);
  arr_1d<12, fem::real_star_8> ts(fem::fill0);
  arr<fem::real_star_8, 4> c0(dimension(3, 3, 3, 3), fem::fill0);
  arr<fem::real_star_8, 4> c(dimension(3, 3, 3, 3), fem::fill0);
  arr_2d<3, 12, fem::real_star_8> slpdir0(fem::fill0);
  arr_2d<3, 12, fem::real_star_8> slpnor0(fem::fill0);
  arr_1d<12, fem::real_star_8> ncrss(fem::fill0);
  int nincr = fem::int0;
  int nincrt = fem::int0;
  fem::real_star_8 dtincr = fem::zero<fem::real_star_8>();
  arr_2d<3, 3, fem::real_star_8> fel(fem::fill0);
  arr_2d<3, 3, fem::real_star_8> f(fem::fill0);
  arr_2d<3, 3, fem::real_star_8> f1(fem::fill0);
  arr_2d<3, 3, fem::real_star_8> eel(fem::fill0);
  arr_2d<3, 3, fem::real_star_8> spk2(fem::fill0);
  arr_2d<3, 3, fem::real_star_8> felinv(fem::fill0);
  arr_2d<3, 3, fem::real_star_8> s(fem::fill0);
  arr_2d<3, 12, fem::real_star_8> slpdir1(fem::fill0);
  arr_2d<3, 12, fem::real_star_8> slpnor1(fem::fill0);
  arr_1d<12, fem::real_star_8> tau(fem::fill0);
  fem::real_star_8 etha = fem::zero<fem::real_star_8>();
  fem::real_star_8 fw = fem::zero<fem::real_star_8>();
  fem::real_star_8 d = fem::zero<fem::real_star_8>();
  arr_1d<12, fem::real_star_8> gdot0(fem::fill0);
  arr_1d<12, fem::real_star_8> gdot(fem::fill0);
  int l = fem::int0;
  arr_2d<12, 12, fem::real_star_8> dtau(fem::fill0);
  int i = fem::int0;
  arr_1d<12, fem::real_star_8> pgdot(fem::fill0);
  arr_1d<12, fem::real_star_8> dg0(fem::fill0);
  arr_1d<12, fem::real_star_8> func0(fem::fill0);
  int j = fem::int0;
  arr_2d<12, 12, fem::real_star_8> djac(fem::fill0);
  bool converged = fem::bool0;
  arr_2d<3, 3, fem::real_star_8> vg(fem::fill0);
  arr_2d<3, 3, fem::real_star_8> fpinv(fem::fill0);
  arr_1d<12, fem::real_star_8> b(fem::fill0);
  arr_1d<12, fem::real_star_8> ro(fem::fill0);
  arr_1d<12, fem::real_star_8> func(fem::fill0);
  fem::real_star_8 l1 = fem::zero<fem::real_star_8>();
  arr_1d<12, fem::real_star_8> grad(fem::fill0);
  arr_2d<12, 12, fem::real_star_8> djac0(fem::fill0);
  arr_1d<12, fem::real_star_8> indx(fem::fill0);
  fem::real_star_8 sseold = fem::zero<fem::real_star_8>();
  fem::real_star_8 sseref = fem::zero<fem::real_star_8>();
  int m1 = fem::int0;
  arr_1d<12, fem::real_star_8> fshill(fem::fill0);
  //C
  //CREAL*8 (A-H,O-Z)
  //C
  if (noel+npt==0)
  {
      std::cout << "inside umat\n";
  }
  cts(48) = mi_ea;
  //C
  if (h1 == 999) {
    cts(1) = props(1);
    cts(2) = props(2);
    cts(3) = props(3);
    rotationpoly(rotate, cts, statev);
  }
  else {
    cts(1) = h1;
    cts(2) = h2;
    cts(3) = h3;
    if (cts(48) == 1) {
      //C: ROTATION MATRIX FROM CRYSTAL SYSTEM TO GOLBAL
      rotation(rotate, cts);
    }
    else {
      //C: ROTATION MATRIX FROM CRYSTAL SYSTEM TO GOLBAL
      rotationpoly(rotate, cts, statev);
      //C                    CALL ROTATION(ROTATE,CTS)
    }
  }
  cts(4) = c11;
  cts(5) = c12;
  cts(6) = c44;
  cts(7) = skcy;
  cts(8) = shm;
  cts(9) = bv;
  cts(10) = dkm;
  cts(11) = ysg;
  cts(12) = ysc;
  cts(13) = ale;
  cts(14) = aii;
  cts(15) = gammacr;
  cts(16) = s1212;
  cts(17) = fh;
  cts(18) = gp;
  cts(19) = bf0;
  cts(20) = bfi;
  cts(21) = f0;
  cts(22) = bk;
  cts(23) = t;
  cts(24) = pp;
  cts(25) = pq;
  cts(26) = s0t;
  cts(27) = shm0;
  cts(28) = af;
  cts(29) = phi;
  cts(30) = vcs0;
  cts(31) = dd0;
  cts(32) = ss;
  cts(33) = fwi;
  cts(34) = fwii;
  cts(35) = fwiii;
  cts(36) = ethai;
  cts(37) = ethaii;
  cts(38) = dgm;
  cts(39) = dgmpsb;
  cts(40) = etha_vein;
  cts(41) = etha_psb;
  cts(42) = etha_cell;
  cts(43) = etha_lab;
  cts(44) = load;
  cts(45) = skmon;
  cts(46) = ovcyn;
  cts(47) = crrot;
  //C
  //C     STATEV:
  //C               STATEV(1)-STATEV(9)      FPINV(3,3)
  //C               STATEV(10)-STATEV(18)    EPL(3,3) PLASTIC STRAIN
  //C               STATEV(19)               E PLEFF , THE EFFECTIVE SHEAR STRAIN OF EPL(3,3) (IN THIS INCEREMENT)
  //C               STATEV(20)               EPLT , THE ACCUMULATED EFFECTIVE SHEAR STRAIN OF EPL(3,3)  (TOTAL)
  //C               STATEV(21)-STATEV(32)    DGAMMAT, TOTAL SHEAR IN EACH SLIP SYSTEM
  //C               STATEV(33)-STATEV(44)    RO    DISLOCATION DENSITY
  //C               STATEV(45)-STATEV(56)    B     BACK STRESS
  //C               STATEV(57)-STATEV(68)    TAU   SHEAR STRESS
  //C               STATEV(69)-STATEV(80)    TS    ATHERMAL STRESS
  //C               STATEV(81)    CROSS SLIP STRAIN
  //C               STATEV(82)    COTTRELL STRAIN
  //C               STATEV(83)    HIRTH STRAIN
  //C               STATEV(84)    ETHA
  //C               STATEV(85)    FW
  //C               STATEV(86)-STATEV(97)    UNLOAD DISLOCATION  ANNIHILATION
  //C               STATEV(100)              CELL SIZE
  //C               STATEV(101)-STATEV(136)  THE SHEAR VALUES IN EACH SLIP SYSTEM TO COMPUTE SHEAR AMPLITUDE
  //C               STATEV(137)-STATEV(148)  SHEAR AMPLITUDE IN EACH SLIP SYSTEM
  //C               STATEV(149)-STATEV(160)  SLIP SYSTEM SHEAR IN OVERLOADS
  //C               STATEV(161)-STATEV(162)  ETHA AND FW IN OVERLOADS
  //C               STATEV(281)-STATEV(328)  SHEAR STRESS RANGES
  //C               STATEV(329)-STATEV(340)  FSHILL
  //C
  //C         IF (TIME(1) .GT. 0.1  .AND. TIME(1) .LT. 0.15) THEN
  //C      OPEN (unit=24,file ='R:\Main Folder\Base Code\A1.txt')
  //C        DO I=1,3
  //C      WRITE (24, 1000) (F(I,J),J=1,3),TIME(2)
  //C 1000 FORMAT ( 12(1x, e13.4) )
  //C        END DO
  //C      END IF
  //C
  //C: INITIAL CONDITION
  if (time(2) == 0.0f) {
    initial(b0, ro0, tau0, statev, fpinv0, cts);
  }
  else {
    //C: READING STATE VARIABLES
    firstread(statev, fpinv0, ro0, b0, ts, tau0);
  }
  //C: INITIAL FORTH RANK ELASTIC TENSOR
  elastic0(c0, cts);
  //C: ROTATE LOCAL FORTH RANK ELASTIC TENSOR INTO THE GLOBAL LOADING SYSTEM
  elastic(c, c0, rotate);
  //C: SLIP SYSTEMS IN GLOBAL LOADING DIRECTION (THE EFFECT OF CRYSTAL ORIENTATION)
  slip0(slpdir0, slpnor0, rotate, cts, ncrss);
  //C: DEFORMATION GRADIENT F AND FEL
  nincr = 1;
  nincrt = 1;
  statement_100:
  dtincr = fem::real(dtime / fem::real(nincrt));
  dgrad(fel, f, f1, dfgrd0, dfgrd1, fpinv0, nincr, nincrt);
  //C: GREEN STRAIN TENSOR
  green(eel, fel);
  //C: SECOND PIOLA-KIRCHHOF STRESS
  piola(spk2, c, eel);
  //C: INVERSE OF ELASTIC GRADIENT  FEL-1(t+dt)
  inverse3(fel, felinv);
  //C: CAUCHY STRESS
  cauchy(s, fel, felinv, spk2);
  //C: CURRENT SLIP SYSTEMS
  slip1(slpdir1, slpnor1, slpdir0, slpnor0, fel, felinv);
  //C: RESOLVED SHEAR STRESS
  rshrs(tau, s, slpdir0, slpnor0);
  //C
  //C: STRUCTURE
  if (cts(44) == 1.0f) {
    structmon(statev, etha, fw, cts);
    cellsizemon(d, tau, cts, statev, kstep, time);
  }
  else {
    structcyc(statev, etha, fw, cts, kstep, props);
    cellsize(d, tau, cts, statev, kstep, time, props);
  }
  //C
  //C: GAMA-DOT0   RO*L*B*Vg
  gd0(gdot0, d, etha, ro0, cts);
  //C: TS: THRESHOLD STRESS, OR SLIP RESISTANCE
  threshstress(ts, d, ro0, cts);
  //C: INITIAL GAMA DOT FOR SLIP SYSTEMS
  inigd(gdot, gdot0, tau, b0, ts, cts, l);
  //C: DERIVITIVE OF SHEAR STRESS WRT GAMADOT
  dshrsdgd(dtau, c0, slpdir0, slpnor0, dtincr);
  //C-------------------------------------------------------------------------------------------------------
  //C: INITIATE THE NEWTON-RAPHSON ITERATION
  FEM_DO_SAFE(i, 1, 12) {
    pgdot(i) = 0.0f;
    dg0(i) = 0.0f;
    func0(i) = 0.0f;
    FEM_DO_SAFE(j, 1, 12) {
      djac(i, j) = 0.0f;
    }
  }
  converged = false;
  while (!converged) {
    converged = true;
    //C: L VELOCITY GRADIENT  SUM OF GDOT*S0*M0
    vgradient(vg, gdot, slpdir0, slpnor0);
    //C:  PLASTIC PART OF DEFORMATION GRADIENT  FP(t+dt) AND FINALLY INVERSE OF PLASTIC GRADIENT  FP-1(t+dt)
    plasticgradientinv(fpinv, fpinv0, vg, dtincr);
    //C: ELASTIC PART OF DEFORMATION GRADIENT  FE(t+dt)
    elasticgradient(fel, f1, fpinv);
    //C: GREEN STRAIN TENSOR
    green(eel, fel);
    //C: SECOND PIOLA-KIRCHHOF STRESS
    piola(spk2, c, eel);
    //C: INVERSE OF ELASTIC GRADIENT  FEL-1(t+dt)
    inverse3(fel, felinv);
    //C: CAUCHY STRESS
    cauchy(s, fel, felinv, spk2);
    //C: CURRENT SLIP SYSTEMS
    slip1(slpdir1, slpnor1, slpdir0, slpnor0, fel, felinv);
    //C: RESOLVED SHEAR STRESS
    rshrs(tau, s, slpdir0, slpnor0);
    //C: B BACK STRESS       +
    backstress(b, b0, gdot, cts, dtincr, fw, statev, etha, kstep);
    //C: UNLOAD DISLOCATION ANNIHILATION
    dannunload(tau, tau0, statev, b, b0, d, cts, dtincr);
    //C: DISLOCATION DENSITY
    ddensity(ro, ro0, statev, d, etha, gdot, tau, b, cts, dtincr, ncrss);
    //C: TS: THRESHOLD STRESS, OR SLIP RESISTANCE
    threshstress(ts, d, ro, cts);
    //C: GAMA-DOT0   RO*L*B*Vg
    gd0(gdot0, d, etha, ro, cts);
    //C:  FUNCTION FOR NEWTON RAPHSON
    nrfunction(func, gdot, gdot0, tau, b, ts, ro, d, sse, cts, l1);
    //C: ANALYTICAL AND NUMERICAL JACOBIAN MATRIX FOR NEWTON RAPHSON
    if (cts(32) == 0.0f) {
      ajacobian(djac, d, etha, ro, gdot, gdot0, tau, dtau, b, ts,
        statev, func, grad, cts, dtincr, fw, ro0);
    }
    else {
      njacobian(djac, func, gdot, pgdot, dg0, func0, grad);
    }
    //C
    FEM_DO_SAFE(i, 1, 12) {
      func0(i) = func(i);
      pgdot(i) = gdot(i);
      FEM_DO_SAFE(j, 1, 12) {
        djac0(i, j) = djac(i, j);
      }
    }
    //C
    //C:  SOLVE J X DELG = FUNC AND RESULT ANALYSIS
    ludcmp(12, djac0, indx);
    lubksb(djac0, 12, indx, func0);
    //C
    //C: NEW FUNC VALUE AND RESULTS ANALYSIS
    check(gdot, slpdir0, slpnor0, dtincr, f1, c, statev, cts, sse,
      grad, func0, fpinv0, sseold, sseref, fel, spk2, tau, tau0, b,
      ts, fpinv, ro, s, ro0, b0, etha, l1, fw, nincrt, dtime, d,
      ncrss, kstep, m1);
    //C
    //C: CHECKING CONVERGENCE
    if (sseold > 1.0e-5f) {
      converged = false;
    }
    if ((sseold > sseref / 2.0f) && (!converged)) {
      nincr = 2 * nincr - 1;
      nincrt = 2 * nincrt;
      //C
      if (nincrt > 1.0e5f) {
        pnewdt = 0.75f;
        return;
      }
      //C
      goto statement_100;
    }
    //C
    //C  NEWTON RAPHSON END
  }
  //C
  //C: SUBINCREMENT
  if (nincr < nincrt) {
    if (nincr == (nincr / 2) * 2) {
      nincr = (nincr / 2) + 1;
      nincrt = nincrt / 2;
    }
    else {
      nincr++;
    }
    FEM_DO_SAFE(i, 1, 3) {
      FEM_DO_SAFE(j, 1, 3) {
        fpinv0(i, j) = fpinv(i, j);
      }
    }
    FEM_DO_SAFE(i, 1, 12) {
      b0(i) = b(i);
      ro0(i) = ro(i);
      tau0(i) = tau(i);
      statev(i + 20) += dtincr * gdot(i);
    }
    goto statement_100;
  }
  //C
  //C: HILL FACTOR
  hill(fshill, tau, tau0, statev, gdot, dtincr, d, nincrt);
  //C: UPDATE STRESS
  stressupdate(stress, s, ndi, nshr);
  //C: UPDATE MATERIAL TANGENT STIFFNESS MATRIX AND STRESS UPDATE
  stff(ddsdde, c, gdot, tau, b, ts, slpdir0, slpnor0, fel, felinv, cts, dtime);
  //C: STORE STATE VARIABLES EPLEFF=EFFECTIVE PLASTIC SHEAR STRAIN,
  endstore(gdot, statev, fpinv, ro, b, tau, ts, slpdir0, slpnor0,
    dtincr, kstep, rotate);
  //C
  //C:  END OF UMAT
  //C
  //C 1ST MI/ PHI_1/ FOR POLY, H1
  //C 2ND MI/ PHI
  //C 3RD MI/ PHI_2
  //C ELASTIC CONSTAN
  //C CYCLIC SIMILITUDE CONSTAN
  //C SHEAR MODULUS
  //C BERGERS VECTOR
  //C DISLOCATION K MULTI
  //C EDGE DISLOCATION ANNIHIL
  //C SCREW DISLOCATION ANNIHI
  //C ALPHA LE, LINE TENSION COEFF
  //C AII, DISLOCATION-DISLOCATI
  //C CRITICAL CROSS SLIP STRAIN
  //C BACK STRESS COEFICIENT
  //C FSHILL
  //C WALL STRUCTURE EXP PARA
  //C WALL VOLUME FRACTION AT T
  //C SLIP ACTIVATION ENERGY
  //C BOLTZMAN CONSTANT
  //C TEMPERATURE
  //C THE EXPONENT OF CONSTITU
  //C THE EXPONENT OF CONSTITUT
  //C THERMAL SLIP RESISTA
  //C SHEAR MODULUS AT
  //C ATTEMPT FREQU
  //C CROSS SLIP EFFICIENCY PAR
  //C INITIAL CROSS SLIP ACT
  //C INITIAL DISLOCATIO
  //C SOLUTION STRATEGY  SS=1.0 FOR BR
  //C MONOTONIC STRUCTUR
  //C MONOTONIC STRUCTU
  //C MONOTONIC STRUCT
  //C MONOTONIC STRUCTUR
  //C MONOTONIC STRUCTU
  //C STRUCTURE PARAMETE
  //C PSB STRUCTURE PARAM
  //C VEIN STRUCTURE ETHA
  //C PSB STRUCTURE ETHA
  //C CELL STRUCTURE ETHA
  //C LABYRINTHIC STRUCTU
  //C loading condition L
  //C MONOTONIC SIMILIT
  //C OVERLOAD CYCLE NUMB
  //C CRYSTAL ROTATION
  //C MILLER INDICES=1, EULER ANGLES=0
}

fem::real_star_8
log2(
  fem::real_star_8 const& x)
{
  fem::real_star_8 return_value = fem::zero<fem::real_star_8>();
  //C
  return_value = fem::log(x) / fem::log(2.0e0);
  return return_value;
}

fem::real_star_8
pade_coefficient(
  fem::integer_star_4 const& j)
{
  fem::real_star_8 return_value = fem::zero<fem::real_star_8>();
  arr_1d<7, fem::real_star_8> pt_table(fem::fill0);
  pt_table = ( / 1.0e0, 1.0e0 / 2.0e0, 5.0e0 / 44.0e0, 1.0e0 / 66.0e0,
    1.0e0 / 792.0e0, 1.0e0 / 15840.0e0, 1.0e0 / 665280.0e0 / );
  //C
  return_value = pt_table(j + 1);
  return return_value;
}

void
pade_factor(
  arr_cref<fem::real_star_8, 2> a,
  arr_ref<fem::real_star_8, 2> s)
{
  a(dimension(3, 3));
  s(dimension(3, 3));
  fem::integer_star_4 ii = fem::zero<fem::integer_star_4>();
  fem::integer_star_4 jj = fem::zero<fem::integer_star_4>();
  arr_2d<3, 3, fem::real_star_8> unit(fem::fill0);
  FEM_DO_SAFE(ii, 1, 3) {
    FEM_DO_SAFE(jj, 1, 3) {
      unit(ii, jj) = 1.e0;
    }
  }
  //C
  s = unit;
  fem::integer_star_4 j = fem::zero<fem::integer_star_4>();
  arr_2d<3, 3, fem::real_star_8> matmult(fem::fill0);
  FEM_DOSTEP(j, 6, 1, -1) {
    s = unit + (pade_coefficient(j) / pade_coefficient(j - 1)) * matmult(a, s);
  }
  //C
}

void
approximate_exp(
  arr_cref<fem::real_star_8, 2> a,
  arr_ref<fem::real_star_8, 2> exp_a)
{
  a(dimension(3, 3));
  exp_a(dimension(3, 3));
  //C
  bool scale = false;
  fem::real_star_8 maxavalue = fem::max(fem::abs(a));
  int matrixpower = fem::int0;
  fem::real_star_8 modmaxavalue = fem::zero<fem::real_star_8>();
  arr_2d<3, 3, fem::real_star_8> a_mod(fem::fill0);
  if (maxavalue > 0.5f) {
    scale = true;
    matrixpower = fem::nint(log2(2.0e0 * maxavalue) + 0.e5);
    modmaxavalue = fem::pow(2.0e0, matrixpower);
    a_mod = a / modmaxavalue;
  }
  else {
    a_mod = a;
  }
  arr_2d<3, 3, fem::real_star_8> d(fem::fill0);
  pade_factor(-a_mod, d);
  arr_2d<3, 3, fem::real_star_8> n(fem::fill0);
  pade_factor(a_mod, n);
  arr_2d<3, 3, fem::real_star_8> d_inv(fem::fill0);
  inverse3(d, d_inv);
  arr_2d<3, 3, fem::real_star_8> matmult(fem::fill0);
  exp_a = matmult(d_inv, n);
  //C
  int loopvar = fem::int0;
  if (scale) {
    FEM_DO_SAFE(loopvar, 1, matrixpower) {
      exp_a = matmult(exp_a, exp_a);
    }
  }
  //C
}

void
plasticgradientinv1(
  arr_ref<fem::real_star_8, 2> fpinv,
  arr_cref<fem::real_star_8, 2> fpinv0,
  arr_cref<fem::real_star_8, 2> vg,
  fem::real_star_8 const& dtincr)
{
  fpinv(dimension(3, 3));
  fpinv0(dimension(3, 3));
  vg(dimension(3, 3));
  //C
  fem::real_star_8 sum = 0.0f;
  int i = fem::int0;
  int k = fem::int0;
  arr_2d<3, 3, fem::real_star_8> a(fem::fill0);
  FEM_DO_SAFE(i, 1, 3) {
    FEM_DO_SAFE(k, 1, 3) {
      a(i, k) = vg(i, k) * dtincr;
      sum += fem::abs(vg(i, k) * dtincr);
    }
  }
  a(1, 1) += 1.0f;
  a(2, 2) += 1.0f;
  a(3, 3) += 1.0f;
  //C
  arr_2d<3, 3, fem::real_star_8> a1(fem::fill0);
  int n = fem::int0;
  arr_2d<3, 3, fem::real_star_8> a2(fem::fill0);
  arr_2d<3, 3, fem::real_star_8> a3(fem::fill0);
  arr_2d<3, 3, fem::real_star_8> matmult(fem::fill0);
  arr_2d<3, 3, fem::real_star_8> a4(fem::fill0);
  if (sum < 1.0e-5f) {
    //C
    inverse3(a, a1);
    //C
    FEM_DO_SAFE(i, 1, 3) {
      FEM_DO_SAFE(k, 1, 3) {
        fpinv(i, k) = 0;
        FEM_DO_SAFE(n, 1, 3) {
          fpinv(i, k) += fpinv0(i, n) * a1(n, k);
        }
      }
    }
    //C
  }
  else {
    //C
    inverse3(fpinv0, a2);
    approximate_exp(vg * dtincr, a3);
    //C
    a4 = matmult(a3, a2);
    inverse3(a4, fpinv);
  }
  //C
}



} // namespace cp_fcc

//g++ -o "fable_cout" -Wall -Wno-sign-compare -Winvalid-pch -Wno-deprecated-declarations -g -O0 -I"/home/atallman/PycharmProjects/fable/modules/cctbx_project/fable" -I"/home/atallman/PycharmProjects/fable/modules/cctbx_project" "fable_cout.cpp"
//g++ -o "fable_cout" -Wall -Wno-sign-compare -Winvalid-pch -Wno-deprecated-declarations -g -O0 -I"/home/atallman/PycharmProjects/fable/modules/cctbx_project/fable" -I"/home/atallman/PycharmProjects/fable/modules/cctbx_project" "fable_cout.cpp"
