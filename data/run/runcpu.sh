#!/bin/bash
#SBATCH -N1
#SBATCH --time 00:10:00
#SBATCH --export="NONE"
#SBATCH --partition=par
#SBATCH -o runcpu.%j.out

#ml prgenv/nvidia nvidia/22.11 netcdf4/4.9.1 nco/4.9.7 fcm/2021.05.0 openblas/0.3.21
ml prgenv/nvidia nvidia/24.1 netcdf4 nco fcm openblas

set -x
set -e

pwd

bin=../../bin/ecrad_ifs_blocked

ls -l $bin
ldd $bin

export DR_HOOK=0
export DR_HOOK_OPT=prof

export ARCH=CPU
export LACC=F
arch="${ARCH,,}"

dir=run${arch}.$SLURM_JOBID.dir
mkdir -p $dir

for test in small large
do

$bin \
  configCY47R3_${arch}.nam \
  $test.in.nc $dir/$test.out.nc

ncdiff \
  $dir/$test.out.nc \
  ref${arch}/$test.ref.nc \
  $dir/$test.diff.nc

ncdump $dir/$test.diff.nc  | ./filter.pl

done
