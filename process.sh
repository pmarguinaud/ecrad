#!/bin/bash

set -e

dir=$1

if [ "x$dir" = "x" ]
then
  echo "Usage: $0 dir"
  exit 1
fi

all=$2

if [ "x$all" = "x" ]
then
  all=0
fi

for f in \
   driver/ecrad_driver_read_input.F90                        \
   driver/ecrad_ifs_driver.F90                               \
   driver/ecrad_ifs_driver_blocked.F90                       \
   driver/ifs_blocking.F90                                   \
   ifs/cloud_overlap_decorr_len.F90                          \
   ifs/ice_effective_radius.F90                              \
   ifs/liquid_effective_radius.F90                           \
   ifs/radiation_scheme.F90                                  \
   ifs/radiation_setup.F90                                   \
   ifs/satur.F90                                             \
   ifsrrtm/rrtm_cmbgb1.F90                                   \
   ifsrrtm/rrtm_cmbgb10.F90                                  \
   ifsrrtm/rrtm_cmbgb11.F90                                  \
   ifsrrtm/rrtm_cmbgb12.F90                                  \
   ifsrrtm/rrtm_cmbgb13.F90                                  \
   ifsrrtm/rrtm_cmbgb14.F90                                  \
   ifsrrtm/rrtm_cmbgb15.F90                                  \
   ifsrrtm/rrtm_cmbgb16.F90                                  \
   ifsrrtm/rrtm_cmbgb2.F90                                   \
   ifsrrtm/rrtm_cmbgb3.F90                                   \
   ifsrrtm/rrtm_cmbgb4.F90                                   \
   ifsrrtm/rrtm_cmbgb5.F90                                   \
   ifsrrtm/rrtm_cmbgb6.F90                                   \
   ifsrrtm/rrtm_cmbgb7.F90                                   \
   ifsrrtm/rrtm_cmbgb8.F90                                   \
   ifsrrtm/rrtm_cmbgb9.F90                                   \
   ifsrrtm/rrtm_gas_optical_depth.F90                        \
   ifsrrtm/rrtm_init_140gp.F90                               \
   ifsrrtm/rrtm_prepare_gases.F90                            \
   ifsrrtm/rrtm_setcoef_140gp.F90                            \
   ifsrrtm/rrtm_taumol1.F90                                  \
   ifsrrtm/rrtm_taumol10.F90                                 \
   ifsrrtm/rrtm_taumol11.F90                                 \
   ifsrrtm/rrtm_taumol12.F90                                 \
   ifsrrtm/rrtm_taumol13.F90                                 \
   ifsrrtm/rrtm_taumol14.F90                                 \
   ifsrrtm/rrtm_taumol15.F90                                 \
   ifsrrtm/rrtm_taumol16.F90                                 \
   ifsrrtm/rrtm_taumol2.F90                                  \
   ifsrrtm/rrtm_taumol3.F90                                  \
   ifsrrtm/rrtm_taumol4.F90                                  \
   ifsrrtm/rrtm_taumol5.F90                                  \
   ifsrrtm/rrtm_taumol6.F90                                  \
   ifsrrtm/rrtm_taumol7.F90                                  \
   ifsrrtm/rrtm_taumol8.F90                                  \
   ifsrrtm/rrtm_taumol9.F90                                  \
   ifsrrtm/srtm_cmbgb16.F90                                  \
   ifsrrtm/srtm_cmbgb17.F90                                  \
   ifsrrtm/srtm_cmbgb18.F90                                  \
   ifsrrtm/srtm_cmbgb19.F90                                  \
   ifsrrtm/srtm_cmbgb20.F90                                  \
   ifsrrtm/srtm_cmbgb21.F90                                  \
   ifsrrtm/srtm_cmbgb22.F90                                  \
   ifsrrtm/srtm_cmbgb23.F90                                  \
   ifsrrtm/srtm_cmbgb24.F90                                  \
   ifsrrtm/srtm_cmbgb25.F90                                  \
   ifsrrtm/srtm_cmbgb26.F90                                  \
   ifsrrtm/srtm_cmbgb27.F90                                  \
   ifsrrtm/srtm_cmbgb28.F90                                  \
   ifsrrtm/srtm_cmbgb29.F90                                  \
   ifsrrtm/srtm_gas_optical_depth.F90                        \
   ifsrrtm/srtm_init.F90                                     \
   ifsrrtm/srtm_setcoef.F90                                  \
   ifsrrtm/srtm_taumol16.F90                                 \
   ifsrrtm/srtm_taumol17.F90                                 \
   ifsrrtm/srtm_taumol18.F90                                 \
   ifsrrtm/srtm_taumol19.F90                                 \
   ifsrrtm/srtm_taumol20.F90                                 \
   ifsrrtm/srtm_taumol21.F90                                 \
   ifsrrtm/srtm_taumol22.F90                                 \
   ifsrrtm/srtm_taumol23.F90                                 \
   ifsrrtm/srtm_taumol24.F90                                 \
   ifsrrtm/srtm_taumol25.F90                                 \
   ifsrrtm/srtm_taumol26.F90                                 \
   ifsrrtm/srtm_taumol27.F90                                 \
   ifsrrtm/srtm_taumol28.F90                                 \
   ifsrrtm/srtm_taumol29.F90                                 \
   ifsrrtm/surrtpk.F90                                       \
   ifsrrtm/surrtrf.F90                                       \
   ifsrrtm/susrtm.F90                                        \
   radiation/radiation_adding_ica_lw.F90                     \
   radiation/radiation_adding_ica_sw.F90                     \
   radiation/radiation_aerosol.F90                           \
   radiation/radiation_aerosol_optics.F90                    \
   radiation/radiation_aerosol_optics_data.F90               \
   radiation/radiation_cloud.F90                             \
   radiation/radiation_cloud_cover.F90                       \
   radiation/radiation_cloud_generator.F90                   \
   radiation/radiation_cloud_generator_acc.F90               \
   radiation/radiation_cloud_optics.F90                      \
   radiation/radiation_cloud_optics_data.F90                 \
   radiation/radiation_cloudless_lw.F90                      \
   radiation/radiation_cloudless_sw.F90                      \
   radiation/radiation_config.F90                            \
   radiation/radiation_ecckd.F90                             \
   radiation/radiation_ecckd_gas.F90                         \
   radiation/radiation_ecckd_interface.F90                   \
   radiation/radiation_flux.F90                              \
   radiation/radiation_gas.F90                               \
   radiation/radiation_general_cloud_optics.F90              \
   radiation/radiation_general_cloud_optics_data.F90         \
   radiation/radiation_homogeneous_lw.F90                    \
   radiation/radiation_homogeneous_sw.F90                    \
   radiation/radiation_ice_optics_baran.F90                  \
   radiation/radiation_ice_optics_baran2016.F90              \
   radiation/radiation_ice_optics_baran2017.F90              \
   radiation/radiation_ice_optics_fu.F90                     \
   radiation/radiation_ice_optics_yi.F90                     \
   radiation/radiation_ifs_rrtm.F90                          \
   radiation/radiation_ifs_rrtmgp.F90                        \
   radiation/radiation_interface.F90                         \
   radiation/radiation_liquid_optics_slingo.F90              \
   radiation/radiation_liquid_optics_socrates.F90            \
   radiation/radiation_lw_derivatives.F90                    \
   radiation/radiation_mcica_acc_lw.F90                      \
   radiation/radiation_mcica_acc_sw.F90                      \
   radiation/radiation_mcica_lw.F90                          \
   radiation/radiation_mcica_sw.F90                          \
   radiation/radiation_monochromatic.F90                     \
   radiation/radiation_overlap.F90                           \
   radiation/radiation_pdf_sampler.F90                       \
   radiation/radiation_random_numbers.F90                    \
   radiation/radiation_regions.F90                           \
   radiation/radiation_save.F90                              \
   radiation/radiation_single_level.F90                      \
   radiation/radiation_spartacus_lw.F90                      \
   radiation/radiation_spartacus_sw.F90                      \
   radiation/radiation_spectral_definition.F90               \
   radiation/radiation_thermodynamics.F90                    \
   radiation/radiation_tripleclouds_lw.F90                   \
   radiation/radiation_tripleclouds_sw.F90                   \
   radiation/radiation_two_stream.F90                        

do

  echo "==> $f <=="

  if [ /home/sor/ecrad/mix/$f -ot $f ] || [ $all -eq 1 ]
  then

    ./split.pl $f
    ./lacc.pl $f
    ./merge.pl $f $dir

    b=$(basename $f .F90)
    h=include/$b.intfb.h

    if [ -f $h ]
    then
      echo "==> $h <=="
      ./split.pl $h
      ./lacc.pl $h
      ./merge.pl $h $dir
    fi

  fi

done


