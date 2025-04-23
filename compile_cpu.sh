#!/bin/bash


#ml prgenv/nvidia nvidia/22.11 netcdf4/4.9.1 nco/4.9.7 fcm/2021.05.0 openblas/0.3.21
ml prgenv/nvidia nvidia/24.1 netcdf4 nco fcm openblas

set -x

mkdir -p lib mod

make PROFILE=pgi BITIDENTITY_TESTING=1 BLASLIB=openblas -j32 ifsdriver

