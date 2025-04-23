! radiation_cloud_optics_data.F90 - Type to store cloud optical properties
!
! (C) Copyright 2014- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
! Author:  Robin Hogan
! Email:   r.j.hogan@ecmwf.int
!

#include "ecrad_config.h"

module radiation_cloud_optics_data

  use parkind1

  implicit none
  

  !---------------------------------------------------------------------
  ! This type holds the configuration information to compute
  ! cloud optical properties
  type cloud_optics_type
     ! Band-specific coefficients are provided separately in the
     ! shortwave and longwave, and are dimensioned (nband,ncoeff),
     ! where ncoeff depends on the nature of the parameterization
     real(jprb), allocatable, dimension(:,:) :: &
          &  liq_coeff_lw, liq_coeff_sw, &
          &  ice_coeff_lw, ice_coeff_sw
     ! General coefficients are vectors of length ncoeffgen, which
     ! depends on the nature of the parameterization; note that most
     ! parameterizations use only band-specific coefficients
     real(jprb), allocatable, dimension(:) :: &
          &  liq_coeff_gen, ice_coeff_gen

   contains
     
     
     
     
     

  procedure :: setup_GPU => setup_cloud_optics_data_GPU

  procedure :: create_device_GPU

  procedure :: update_host_GPU

  procedure :: update_device_GPU

  procedure :: delete_device_GPU

  procedure :: setup_CPU => setup_cloud_optics_data_CPU

  end type cloud_optics_type

contains

  !---------------------------------------------------------------------
  ! Setup cloud optics coefficients by reading them from a file
  


  

  

  

  

  subroutine setup_cloud_optics_data_GPU(this, liq_file_name, ice_file_name, iverbose, lacc)

    use yomhook
    use easy_netcdf

    class(cloud_optics_type), intent(inout) :: this
    character(len=*), intent(in)            :: liq_file_name, ice_file_name
    integer, intent(in), optional           :: iverbose

    ! The NetCDF file containing the coefficients
    type(netcdf_file)  :: file
    integer            :: iverb
    real(jphook) :: hook_handle
    logical, intent (in) :: lacc

    if (lhook) call dr_hook('radiation_cloud_optics_data:setup',0,hook_handle)

    if (present(iverbose)) then
      iverb = iverbose
    else
      iverb = 2
    end if

    ! Open the droplet scattering file and configure the way it is
    ! read
    call file%open(trim(liq_file_name), iverbose=iverb)
    call file%transpose_matrices()

    ! Read the band-specific coefficients
    call file%get('coeff_lw',this%liq_coeff_lw)
    call file%get('coeff_sw',this%liq_coeff_sw)

    ! Read the general  coefficients
    if (file%exists('coeff_gen')) then
      call file%get('coeff_gen',this%liq_coeff_gen)
    end if

    ! Close droplet scattering file
    call file%close()

    ! Open the ice scattering file and configure the way it is read
    call file%open(trim(ice_file_name), iverbose=iverb)
    call file%transpose_matrices()

    ! Read the band-specific  coefficients
    call file%get('coeff_lw',this%ice_coeff_lw)
    call file%get('coeff_sw',this%ice_coeff_sw)

    ! Read the general  coefficients
    if (file%exists('coeff_gen')) then
      call file%get('coeff_gen',this%ice_coeff_gen)
    end if

    ! Close ice scattering file
    call file%close()

    if (lhook) call dr_hook('radiation_cloud_optics_data:setup',1,hook_handle)

  end subroutine setup_cloud_optics_data_GPU

  subroutine create_device_GPU(this, lacc)
    class(cloud_optics_type), intent(inout) :: this
    logical, intent (in) :: lacc

    !$ACC ENTER DATA COPYIN(this%liq_coeff_lw) IF((allocated(this%liq_coeff_lw)) .AND. lacc) ASYNC(1)
    !$ACC ENTER DATA COPYIN(this%liq_coeff_sw) IF((allocated(this%liq_coeff_sw)) .AND. lacc) ASYNC(1)
    !$ACC ENTER DATA COPYIN(this%ice_coeff_lw) IF((allocated(this%ice_coeff_lw)) .AND. lacc) ASYNC(1)
    !$ACC ENTER DATA COPYIN(this%ice_coeff_sw) IF((allocated(this%ice_coeff_sw)) .AND. lacc) ASYNC(1)
    !$ACC ENTER DATA COPYIN(this%liq_coeff_gen) IF((allocated(this%liq_coeff_gen)) .AND. lacc) ASYNC(1)
    !$ACC ENTER DATA COPYIN(this%ice_coeff_gen) IF((allocated(this%ice_coeff_gen)) .AND. lacc) ASYNC(1)
  end subroutine create_device_GPU

  subroutine update_host_GPU(this, lacc)
    class(cloud_optics_type), intent(inout) :: this
    logical, intent (in) :: lacc

    !$ACC UPDATE HOST(this%liq_coeff_lw) IF((allocated(this%liq_coeff_lw)) .AND. lacc) ASYNC(1)
    !$ACC UPDATE HOST(this%liq_coeff_sw) IF((allocated(this%liq_coeff_sw)) .AND. lacc) ASYNC(1)
    !$ACC UPDATE HOST(this%ice_coeff_lw) IF((allocated(this%ice_coeff_lw)) .AND. lacc) ASYNC(1)
    !$ACC UPDATE HOST(this%ice_coeff_sw) IF((allocated(this%ice_coeff_sw)) .AND. lacc) ASYNC(1)
    !$ACC UPDATE HOST(this%liq_coeff_gen) IF((allocated(this%liq_coeff_gen)) .AND. lacc) ASYNC(1)
    !$ACC UPDATE HOST(this%ice_coeff_gen) IF((allocated(this%ice_coeff_gen)) .AND. lacc) ASYNC(1)
  end subroutine update_host_GPU

  subroutine update_device_GPU(this, lacc)
    class(cloud_optics_type), intent(inout) :: this
    logical, intent (in) :: lacc

    !$ACC UPDATE DEVICE(this%liq_coeff_lw) IF((allocated(this%liq_coeff_lw)) .AND. lacc) ASYNC(1)
    !$ACC UPDATE DEVICE(this%liq_coeff_sw) IF((allocated(this%liq_coeff_sw)) .AND. lacc) ASYNC(1)
    !$ACC UPDATE DEVICE(this%ice_coeff_lw) IF((allocated(this%ice_coeff_lw)) .AND. lacc) ASYNC(1)
    !$ACC UPDATE DEVICE(this%ice_coeff_sw) IF((allocated(this%ice_coeff_sw)) .AND. lacc) ASYNC(1)
    !$ACC UPDATE DEVICE(this%liq_coeff_gen) IF((allocated(this%liq_coeff_gen)) .AND. lacc) ASYNC(1)
    !$ACC UPDATE DEVICE(this%ice_coeff_gen) IF((allocated(this%ice_coeff_gen)) .AND. lacc) ASYNC(1)
  end subroutine update_device_GPU

  subroutine delete_device_GPU(this, lacc)
    class(cloud_optics_type), intent(inout) :: this
    logical, intent (in) :: lacc

    !$ACC EXIT DATA DELETE(this%liq_coeff_lw) IF((allocated(this%liq_coeff_lw)) .AND. lacc) ASYNC(1)
    !$ACC EXIT DATA DELETE(this%liq_coeff_sw) IF((allocated(this%liq_coeff_sw)) .AND. lacc) ASYNC(1)
    !$ACC EXIT DATA DELETE(this%ice_coeff_lw) IF((allocated(this%ice_coeff_lw)) .AND. lacc) ASYNC(1)
    !$ACC EXIT DATA DELETE(this%ice_coeff_sw) IF((allocated(this%ice_coeff_sw)) .AND. lacc) ASYNC(1)
    !$ACC EXIT DATA DELETE(this%liq_coeff_gen) IF((allocated(this%liq_coeff_gen)) .AND. lacc) ASYNC(1)
    !$ACC EXIT DATA DELETE(this%ice_coeff_gen) IF((allocated(this%ice_coeff_gen)) .AND. lacc) ASYNC(1)
  end subroutine delete_device_GPU

  subroutine setup_cloud_optics_data_CPU(this, liq_file_name, ice_file_name, iverbose)

    use yomhook
    use easy_netcdf

    class(cloud_optics_type), intent(inout) :: this
    character(len=*), intent(in)            :: liq_file_name, ice_file_name
    integer, intent(in), optional           :: iverbose

    ! The NetCDF file containing the coefficients
    type(netcdf_file)  :: file
    integer            :: iverb
    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_cloud_optics_data:setup',0,hook_handle)

    if (present(iverbose)) then
      iverb = iverbose
    else
      iverb = 2
    end if

    ! Open the droplet scattering file and configure the way it is
    ! read
    call file%open(trim(liq_file_name), iverbose=iverb)
    call file%transpose_matrices()

    ! Read the band-specific coefficients
    call file%get('coeff_lw',this%liq_coeff_lw)
    call file%get('coeff_sw',this%liq_coeff_sw)

    ! Read the general  coefficients
    if (file%exists('coeff_gen')) then
      call file%get('coeff_gen',this%liq_coeff_gen)
    end if

    ! Close droplet scattering file
    call file%close()

    ! Open the ice scattering file and configure the way it is read
    call file%open(trim(ice_file_name), iverbose=iverb)
    call file%transpose_matrices()

    ! Read the band-specific  coefficients
    call file%get('coeff_lw',this%ice_coeff_lw)
    call file%get('coeff_sw',this%ice_coeff_sw)

    ! Read the general  coefficients
    if (file%exists('coeff_gen')) then
      call file%get('coeff_gen',this%ice_coeff_gen)
    end if

    ! Close ice scattering file
    call file%close()

    if (lhook) call dr_hook('radiation_cloud_optics_data:setup',1,hook_handle)

  end subroutine setup_cloud_optics_data_CPU

end module radiation_cloud_optics_data
