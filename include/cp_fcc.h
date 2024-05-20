//
// Created by atallman on 5/14/24.
//

#ifndef CP_FCC_H
#define CP_FCC_H

#include "fem.hpp"


namespace cp_fcc {

    using namespace fem::major_types;

    using fem::common;

    void umat(arr_ref<fem::real_star_8> stress,
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
              int const& /* noel */,
              int const& /* npt */,
              fem::real_star_8 const& /* layer */,
              int const& /* kspt */,
              int const& kstep,
              int const& /* kinc */);
}

#endif //CP_FCC_H
