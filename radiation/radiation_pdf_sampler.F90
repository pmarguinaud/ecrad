! radiation_pdf_sampler.F90 - Get samples from a PDF for McICA
!
! (C) Copyright 2015- ECMWF.
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

module radiation_pdf_sampler

  use parkind1

  implicit none
  

  !---------------------------------------------------------------------
  ! Derived type for sampling from a lognormal or gamma distribution,
  ! or other PDF, used to generate water content or optical depth
  ! scalings for use in the Monte Carlo Independent Column
  ! Approximation (McICA)
  type pdf_sampler_type
    ! Number of points in look-up table for cumulative distribution
    ! function (CDF) and fractional standard deviation (FSD)
    ! dimensions
    integer :: ncdf, nfsd

    ! First value of FSD and the reciprocal of the interval between
    ! FSD values (which are assumed to be uniformly distributed)
    real(jprb) :: fsd1, inv_fsd_interval

    ! Value of the distribution for each CDF and FSD bin
    real(jprb), allocatable, dimension(:,:) :: val

  contains

    
    
    
    
    
    

    
    
    
    

  procedure :: setup_GPU => setup_pdf_sampler_GPU

  procedure :: sample_GPU => sample_from_pdf_GPU

  procedure :: masked_sample_GPU => sample_from_pdf_masked_GPU

  procedure :: block_sample_GPU => sample_from_pdf_block_GPU

  procedure :: masked_block_sample_GPU => sample_from_pdf_masked_block_GPU

  procedure :: deallocate_GPU => deallocate_pdf_sampler_GPU

  procedure :: create_device_GPU

  procedure :: update_host_GPU

  procedure :: update_device_GPU

  procedure :: delete_device_GPU

  procedure :: setup_CPU => setup_pdf_sampler_CPU

  procedure :: sample_CPU => sample_from_pdf_CPU

  procedure :: masked_sample_CPU => sample_from_pdf_masked_CPU

  procedure :: block_sample_CPU => sample_from_pdf_block_CPU

  procedure :: masked_block_sample_CPU => sample_from_pdf_masked_block_CPU

  procedure :: deallocate_CPU => deallocate_pdf_sampler_CPU

  end type pdf_sampler_type

contains

  !---------------------------------------------------------------------
  ! Load look-up table from a file
  

  !---------------------------------------------------------------------
  ! Deallocate data in pdf_sampler_type derived type
  


  !---------------------------------------------------------------------
  ! Extract the value from a PDF with fractional standard deviation
  ! "fsd" corresponding to the cumulative distribution function value
  ! "cdf", and return it in val. Since this is an elemental
  ! subroutine, fsd, cdf and val may be arrays.
  


  !---------------------------------------------------------------------
  ! For true elements of mask, extract the values of a PDF with
  ! fractional standard deviation "fsd" corresponding to the
  ! cumulative distribution function values "cdf", and return in
  ! val. For false elements of mask, return zero in val.
  

  !---------------------------------------------------------------------
  ! Extract the values of a PDF with fractional standard deviation
  ! "fsd" corresponding to the cumulative distribution function values
  ! "cdf", and return in val. This version works on 2D blocks of data.
  

  !---------------------------------------------------------------------
  ! Extract the values of a PDF with fractional standard deviation
  ! "fsd" corresponding to the cumulative distribution function values
  ! "cdf", and return in val. This version works on 2D blocks of data.
  


  

  

  

  

  subroutine setup_pdf_sampler_GPU(this, file_name, iverbose, lacc)

    use yomhook
    use easy_netcdf

    class(pdf_sampler_type), intent(inout) :: this
    character(len=*),        intent(in)    :: file_name
    integer, optional,       intent(in)    :: iverbose

    type(netcdf_file)  :: file
    integer            :: iverb
    real(jprb), allocatable :: fsd(:)

    real(jphook) :: hook_handle
    logical, intent (in) :: lacc

    if (lhook) call dr_hook('radiation_pdf_sampler:setup',0,hook_handle)

    if (present(iverbose)) then
      iverb = iverbose
    else
      iverb = 2
    end if

    if (allocated(this%val)) then
      deallocate(this%val)
    end if

    call file%open(trim(file_name), iverbose=iverb)

    call file%get('fsd',fsd)
    call file%get('x',  this%val)

    call file%close()

    this%ncdf = size(this%val,1)
    this%nfsd = size(this%val,2)
    this%fsd1 = fsd(1)
    this%inv_fsd_interval = 1.0_jprb / (fsd(2)-fsd(1))

    deallocate(fsd)

    if (lhook) call dr_hook('radiation_pdf_sampler:setup',1,hook_handle)

  end subroutine setup_pdf_sampler_GPU

  subroutine deallocate_pdf_sampler_GPU(this, lacc)

    use yomhook

    class(pdf_sampler_type), intent(inout) :: this
    real(jphook) :: hook_handle
    logical, intent (in) :: lacc

    if (lhook) call dr_hook('radiation_pdf_sampler:deallocate',0,hook_handle)

    if (allocated(this%val)) then
      deallocate(this%val)
    end if

    if (lhook) call dr_hook('radiation_pdf_sampler:deallocate',1,hook_handle)

  end subroutine deallocate_pdf_sampler_GPU

  elemental subroutine sample_from_pdf_GPU(this, fsd, cdf, val, lacc)

    class(pdf_sampler_type), intent(in)  :: this

    ! Fractional standard deviation (0 to 4) and cumulative
    ! distribution function (0 to 1)
    real(jprb),              intent(in)  :: fsd, cdf

    ! Sample from distribution
    real(jprb),              intent(out) :: val

    ! Index to look-up table
    integer    :: ifsd, icdf

    ! Weights in bilinear interpolation
    real(jprb) :: wfsd, wcdf
    logical, intent (in) :: lacc

    ! Bilinear interpolation with bounds
    wcdf = cdf * (this%ncdf-1) + 1.0_jprb
    icdf = max(1, min(int(wcdf), this%ncdf-1))
    wcdf = max(0.0_jprb, min(wcdf - icdf, 1.0_jprb))

    wfsd = (fsd-this%fsd1) * this%inv_fsd_interval + 1.0_jprb
    ifsd = max(1, min(int(wfsd), this%nfsd-1))
    wfsd = max(0.0_jprb, min(wfsd - ifsd, 1.0_jprb))

    val =    (1.0_jprb-wcdf)*(1.0_jprb-wfsd) * this%val(icdf  ,ifsd)   &
         & + (1.0_jprb-wcdf)*          wfsd  * this%val(icdf  ,ifsd+1) &
         & +           wcdf *(1.0_jprb-wfsd) * this%val(icdf+1,ifsd)   &
         & +           wcdf *          wfsd  * this%val(icdf+1,ifsd+1)

  end subroutine sample_from_pdf_GPU

  subroutine sample_from_pdf_masked_GPU(this, nsamp, fsd, cdf, val, mask, lacc)

    class(pdf_sampler_type), intent(in)  :: this

    ! Number of samples
    integer,    intent(in) :: nsamp

    ! Fractional standard deviation (0 to 4) and cumulative
    ! distribution function (0 to 1)
    real(jprb), intent(in)  :: fsd(nsamp), cdf(nsamp)

    ! Sample from distribution
    real(jprb), intent(out) :: val(:)

    ! Mask
    logical,    intent(in) :: mask(nsamp)

    ! Loop index
    integer    :: jsamp

    ! Index to look-up table
    integer    :: ifsd, icdf

    ! Weights in bilinear interpolation
    real(jprb) :: wfsd, wcdf
    logical, intent (in) :: lacc

    do jsamp = 1,nsamp
      if (mask(jsamp)) then
        ! Bilinear interpolation with bounds
        wcdf = cdf(jsamp) * (this%ncdf-1) + 1.0_jprb
        icdf = max(1, min(int(wcdf), this%ncdf-1))
        wcdf = max(0.0_jprb, min(wcdf - icdf, 1.0_jprb))

        wfsd = (fsd(jsamp)-this%fsd1) * this%inv_fsd_interval + 1.0_jprb
        ifsd = max(1, min(int(wfsd), this%nfsd-1))
        wfsd = max(0.0_jprb, min(wfsd - ifsd, 1.0_jprb))

        val(jsamp)=(1.0_jprb-wcdf)*(1.0_jprb-wfsd) * this%val(icdf  ,ifsd)   &
             &    +(1.0_jprb-wcdf)*          wfsd  * this%val(icdf  ,ifsd+1) &
             &    +          wcdf *(1.0_jprb-wfsd) * this%val(icdf+1,ifsd)   &
             &    +          wcdf *          wfsd  * this%val(icdf+1,ifsd+1)
      else
        val(jsamp) = 0.0_jprb
      end if
    end do
  end subroutine sample_from_pdf_masked_GPU

  subroutine sample_from_pdf_block_GPU(this, nz, ng, fsd, cdf, val, lacc)

    class(pdf_sampler_type), intent(in)  :: this

    ! Number of samples
    integer,    intent(in) :: nz, ng

    ! Fractional standard deviation (0 to 4) and cumulative
    ! distribution function (0 to 1)
    real(jprb), intent(in)  :: fsd(nz), cdf(ng, nz)

    ! Sample from distribution
    real(jprb), intent(out) :: val(:,:)

    ! Loop index
    integer    :: jz, jg

    ! Index to look-up table
    integer    :: ifsd, icdf

    ! Weights in bilinear interpolation
    real(jprb) :: wfsd, wcdf
    logical, intent (in) :: lacc

    do jz = 1,nz
      do jg = 1,ng
        if (cdf(jg, jz) > 0.0_jprb) then
          ! Bilinear interpolation with bounds
          wcdf = cdf(jg,jz) * (this%ncdf-1) + 1.0_jprb
          icdf = max(1, min(int(wcdf), this%ncdf-1))
          wcdf = max(0.0_jprb, min(wcdf - icdf, 1.0_jprb))

          wfsd = (fsd(jz)-this%fsd1) * this%inv_fsd_interval + 1.0_jprb
          ifsd = max(1, min(int(wfsd), this%nfsd-1))
          wfsd = max(0.0_jprb, min(wfsd - ifsd, 1.0_jprb))

          val(jg,jz)=(1.0_jprb-wcdf)*(1.0_jprb-wfsd) * this%val(icdf  ,ifsd)   &
               &    +(1.0_jprb-wcdf)*          wfsd  * this%val(icdf  ,ifsd+1) &
               &    +          wcdf *(1.0_jprb-wfsd) * this%val(icdf+1,ifsd)   &
               &    +          wcdf *          wfsd  * this%val(icdf+1,ifsd+1)
        else
          val(jg,jz) = 0.0_jprb
        end if
      end do
    end do

  end subroutine sample_from_pdf_block_GPU

  subroutine sample_from_pdf_masked_block_GPU(this, nz, ng, fsd, cdf, val, mask, lacc)

    class(pdf_sampler_type), intent(in)  :: this

    ! Number of samples
    integer,    intent(in) :: nz, ng

    ! Fractional standard deviation (0 to 4) and cumulative
    ! distribution function (0 to 1)
    real(jprb), intent(in)  :: fsd(nz), cdf(ng, nz)

    ! Sample from distribution
    real(jprb), intent(out) :: val(:,:)

    ! Mask
    logical,    intent(in), optional :: mask(nz)

    ! Loop index
    integer    :: jz, jg

    ! Index to look-up table
    integer    :: ifsd, icdf

    ! Weights in bilinear interpolation
    real(jprb) :: wfsd, wcdf
    logical, intent (in) :: lacc

    do jz = 1,nz

      if (mask(jz)) then

        do jg = 1,ng
          if (cdf(jg, jz) > 0.0_jprb) then
            ! Bilinear interpolation with bounds
            wcdf = cdf(jg,jz) * (this%ncdf-1) + 1.0_jprb
            icdf = max(1, min(int(wcdf), this%ncdf-1))
            wcdf = max(0.0_jprb, min(wcdf - icdf, 1.0_jprb))

            wfsd = (fsd(jz)-this%fsd1) * this%inv_fsd_interval + 1.0_jprb
            ifsd = max(1, min(int(wfsd), this%nfsd-1))
            wfsd = max(0.0_jprb, min(wfsd - ifsd, 1.0_jprb))

            val(jg,jz)=(1.0_jprb-wcdf)*(1.0_jprb-wfsd) * this%val(icdf  ,ifsd)   &
                 &    +(1.0_jprb-wcdf)*          wfsd  * this%val(icdf  ,ifsd+1) &
                 &    +          wcdf *(1.0_jprb-wfsd) * this%val(icdf+1,ifsd)   &
                 &    +          wcdf *          wfsd  * this%val(icdf+1,ifsd+1)
          else
            val(jg,jz) = 0.0_jprb
          end if
        end do

      end if

    end do

  end subroutine sample_from_pdf_masked_block_GPU

  subroutine create_device_GPU(this, lacc)
    class(pdf_sampler_type), intent(inout) :: this
    logical, intent (in) :: lacc

    !$ACC ENTER DATA COPYIN(this%val) IF((allocated(this%val)) .AND. lacc) ASYNC(1)
  end subroutine create_device_GPU

  subroutine update_host_GPU(this, lacc)
    class(pdf_sampler_type), intent(inout) :: this
    logical, intent (in) :: lacc

    !$ACC UPDATE HOST(this%val) IF((allocated(this%val)) .AND. lacc) ASYNC(1)
  end subroutine update_host_GPU

  subroutine update_device_GPU(this, lacc)
    class(pdf_sampler_type), intent(inout) :: this
    logical, intent (in) :: lacc

    !$ACC UPDATE DEVICE(this%val) IF((allocated(this%val)) .AND. lacc) ASYNC(1)
  end subroutine update_device_GPU

  subroutine delete_device_GPU(this, lacc)
    class(pdf_sampler_type), intent(inout) :: this
    logical, intent (in) :: lacc

    !$ACC EXIT DATA DELETE(this%val) IF((allocated(this%val)) .AND. lacc) ASYNC(1)
  end subroutine delete_device_GPU

  subroutine setup_pdf_sampler_CPU(this, file_name, iverbose)

    use yomhook
    use easy_netcdf

    class(pdf_sampler_type), intent(inout) :: this
    character(len=*),        intent(in)    :: file_name
    integer, optional,       intent(in)    :: iverbose

    type(netcdf_file)  :: file
    integer            :: iverb
    real(jprb), allocatable :: fsd(:)

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_pdf_sampler:setup',0,hook_handle)

    if (present(iverbose)) then
      iverb = iverbose
    else
      iverb = 2
    end if

    if (allocated(this%val)) then
      deallocate(this%val)
    end if

    call file%open(trim(file_name), iverbose=iverb)

    call file%get('fsd',fsd)
    call file%get('x',  this%val)

    call file%close()

    this%ncdf = size(this%val,1)
    this%nfsd = size(this%val,2)
    this%fsd1 = fsd(1)
    this%inv_fsd_interval = 1.0_jprb / (fsd(2)-fsd(1))

    deallocate(fsd)

    if (lhook) call dr_hook('radiation_pdf_sampler:setup',1,hook_handle)

  end subroutine setup_pdf_sampler_CPU

  subroutine deallocate_pdf_sampler_CPU(this)

    use yomhook

    class(pdf_sampler_type), intent(inout) :: this
    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_pdf_sampler:deallocate',0,hook_handle)

    if (allocated(this%val)) then
      deallocate(this%val)
    end if

    if (lhook) call dr_hook('radiation_pdf_sampler:deallocate',1,hook_handle)

  end subroutine deallocate_pdf_sampler_CPU

  elemental subroutine sample_from_pdf_CPU(this, fsd, cdf, val)

    class(pdf_sampler_type), intent(in)  :: this

    ! Fractional standard deviation (0 to 4) and cumulative
    ! distribution function (0 to 1)
    real(jprb),              intent(in)  :: fsd, cdf

    ! Sample from distribution
    real(jprb),              intent(out) :: val

    ! Index to look-up table
    integer    :: ifsd, icdf

    ! Weights in bilinear interpolation
    real(jprb) :: wfsd, wcdf

    ! Bilinear interpolation with bounds
    wcdf = cdf * (this%ncdf-1) + 1.0_jprb
    icdf = max(1, min(int(wcdf), this%ncdf-1))
    wcdf = max(0.0_jprb, min(wcdf - icdf, 1.0_jprb))

    wfsd = (fsd-this%fsd1) * this%inv_fsd_interval + 1.0_jprb
    ifsd = max(1, min(int(wfsd), this%nfsd-1))
    wfsd = max(0.0_jprb, min(wfsd - ifsd, 1.0_jprb))

    val =    (1.0_jprb-wcdf)*(1.0_jprb-wfsd) * this%val(icdf  ,ifsd)   &
         & + (1.0_jprb-wcdf)*          wfsd  * this%val(icdf  ,ifsd+1) &
         & +           wcdf *(1.0_jprb-wfsd) * this%val(icdf+1,ifsd)   &
         & +           wcdf *          wfsd  * this%val(icdf+1,ifsd+1)

  end subroutine sample_from_pdf_CPU

  subroutine sample_from_pdf_masked_CPU(this, nsamp, fsd, cdf, val, mask)

    class(pdf_sampler_type), intent(in)  :: this

    ! Number of samples
    integer,    intent(in) :: nsamp

    ! Fractional standard deviation (0 to 4) and cumulative
    ! distribution function (0 to 1)
    real(jprb), intent(in)  :: fsd(nsamp), cdf(nsamp)

    ! Sample from distribution
    real(jprb), intent(out) :: val(:)

    ! Mask
    logical,    intent(in) :: mask(nsamp)

    ! Loop index
    integer    :: jsamp

    ! Index to look-up table
    integer    :: ifsd, icdf

    ! Weights in bilinear interpolation
    real(jprb) :: wfsd, wcdf

    do jsamp = 1,nsamp
      if (mask(jsamp)) then
        ! Bilinear interpolation with bounds
        wcdf = cdf(jsamp) * (this%ncdf-1) + 1.0_jprb
        icdf = max(1, min(int(wcdf), this%ncdf-1))
        wcdf = max(0.0_jprb, min(wcdf - icdf, 1.0_jprb))

        wfsd = (fsd(jsamp)-this%fsd1) * this%inv_fsd_interval + 1.0_jprb
        ifsd = max(1, min(int(wfsd), this%nfsd-1))
        wfsd = max(0.0_jprb, min(wfsd - ifsd, 1.0_jprb))

        val(jsamp)=(1.0_jprb-wcdf)*(1.0_jprb-wfsd) * this%val(icdf  ,ifsd)   &
             &    +(1.0_jprb-wcdf)*          wfsd  * this%val(icdf  ,ifsd+1) &
             &    +          wcdf *(1.0_jprb-wfsd) * this%val(icdf+1,ifsd)   &
             &    +          wcdf *          wfsd  * this%val(icdf+1,ifsd+1)
      else
        val(jsamp) = 0.0_jprb
      end if
    end do
  end subroutine sample_from_pdf_masked_CPU

  subroutine sample_from_pdf_block_CPU(this, nz, ng, fsd, cdf, val)

    class(pdf_sampler_type), intent(in)  :: this

    ! Number of samples
    integer,    intent(in) :: nz, ng

    ! Fractional standard deviation (0 to 4) and cumulative
    ! distribution function (0 to 1)
    real(jprb), intent(in)  :: fsd(nz), cdf(ng, nz)

    ! Sample from distribution
    real(jprb), intent(out) :: val(:,:)

    ! Loop index
    integer    :: jz, jg

    ! Index to look-up table
    integer    :: ifsd, icdf

    ! Weights in bilinear interpolation
    real(jprb) :: wfsd, wcdf

    do jz = 1,nz
      do jg = 1,ng
        if (cdf(jg, jz) > 0.0_jprb) then
          ! Bilinear interpolation with bounds
          wcdf = cdf(jg,jz) * (this%ncdf-1) + 1.0_jprb
          icdf = max(1, min(int(wcdf), this%ncdf-1))
          wcdf = max(0.0_jprb, min(wcdf - icdf, 1.0_jprb))

          wfsd = (fsd(jz)-this%fsd1) * this%inv_fsd_interval + 1.0_jprb
          ifsd = max(1, min(int(wfsd), this%nfsd-1))
          wfsd = max(0.0_jprb, min(wfsd - ifsd, 1.0_jprb))

          val(jg,jz)=(1.0_jprb-wcdf)*(1.0_jprb-wfsd) * this%val(icdf  ,ifsd)   &
               &    +(1.0_jprb-wcdf)*          wfsd  * this%val(icdf  ,ifsd+1) &
               &    +          wcdf *(1.0_jprb-wfsd) * this%val(icdf+1,ifsd)   &
               &    +          wcdf *          wfsd  * this%val(icdf+1,ifsd+1)
        else
          val(jg,jz) = 0.0_jprb
        end if
      end do
    end do

  end subroutine sample_from_pdf_block_CPU

  subroutine sample_from_pdf_masked_block_CPU(this, nz, ng, fsd, cdf, val, mask)

    class(pdf_sampler_type), intent(in)  :: this

    ! Number of samples
    integer,    intent(in) :: nz, ng

    ! Fractional standard deviation (0 to 4) and cumulative
    ! distribution function (0 to 1)
    real(jprb), intent(in)  :: fsd(nz), cdf(ng, nz)

    ! Sample from distribution
    real(jprb), intent(out) :: val(:,:)

    ! Mask
    logical,    intent(in), optional :: mask(nz)

    ! Loop index
    integer    :: jz, jg

    ! Index to look-up table
    integer    :: ifsd, icdf

    ! Weights in bilinear interpolation
    real(jprb) :: wfsd, wcdf

    do jz = 1,nz

      if (mask(jz)) then

        do jg = 1,ng
          if (cdf(jg, jz) > 0.0_jprb) then
            ! Bilinear interpolation with bounds
            wcdf = cdf(jg,jz) * (this%ncdf-1) + 1.0_jprb
            icdf = max(1, min(int(wcdf), this%ncdf-1))
            wcdf = max(0.0_jprb, min(wcdf - icdf, 1.0_jprb))

            wfsd = (fsd(jz)-this%fsd1) * this%inv_fsd_interval + 1.0_jprb
            ifsd = max(1, min(int(wfsd), this%nfsd-1))
            wfsd = max(0.0_jprb, min(wfsd - ifsd, 1.0_jprb))

            val(jg,jz)=(1.0_jprb-wcdf)*(1.0_jprb-wfsd) * this%val(icdf  ,ifsd)   &
                 &    +(1.0_jprb-wcdf)*          wfsd  * this%val(icdf  ,ifsd+1) &
                 &    +          wcdf *(1.0_jprb-wfsd) * this%val(icdf+1,ifsd)   &
                 &    +          wcdf *          wfsd  * this%val(icdf+1,ifsd+1)
          else
            val(jg,jz) = 0.0_jprb
          end if
        end do

      end if

    end do

  end subroutine sample_from_pdf_masked_block_CPU

end module radiation_pdf_sampler
