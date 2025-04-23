PROGRAM ECRAD_IFS_DRIVER

IMPLICIT NONE

CHARACTER*16 :: CLARCH, CLLACC
LOGICAL :: LACC

CALL GETENV ('ARCH', CLARCH)

LACC = TRIM (CLARCH) == 'GPU'

CALL GETENV ('LACC', CLLACC)
READ (CLLACC, *) LACC

WRITE (0, *) " ARCH = ", CLARCH
WRITE (0, *) " LACC = ", LACC


SELECT CASE (CLARCH)
  CASE ('CPU') 
    CALL ECRAD_IFS_DRIVER_CPU
  CASE ('GPU') 
    CALL ECRAD_IFS_DRIVER_GPU
  CASE DEFAULT
    ERROR STOP "UNEXPECTED ARCH"
END SELECT

CONTAINS


SUBROUTINE ECRAD_IFS_DRIVER_GPU

  ! --------------------------------------------------------
  ! Section 1: Declarations
  ! --------------------------------------------------------
  use parkind1 ! Working/double precision

  use radiation_io
  use radiation_single_level
  use radiation_thermodynamics
  use radiation_gas
  use radiation_cloud
  use radiation_aerosol
  use radiation_flux
  use radiation_save
  use radiation_setup
  use radiation_constants
  use ecrad_driver_config
  use ecrad_driver_read_input
  use easy_netcdf
  use ifs_blocking
  
#ifdef __PGI
  use nvtx
#endif

  implicit none

#include "radiation_scheme.intfb.h"

  ! The NetCDF file containing the input profiles
  type(netcdf_file)         :: file

  ! Configuration for the radiation scheme, IFS style
  type(tradiation)          :: yradiation

  ! Derived types for the inputs to the radiation scheme
  type(single_level_type)   :: single_level
  type(thermodynamics_type) :: thermodynamics
  type(gas_type)            :: gas
  type(cloud_type)          :: cloud
  type(aerosol_type)        :: aerosol

  ! Configuration specific to this driver
  type(driver_config_type)  :: driver_config

  ! Derived type containing outputs from the radiation scheme
  type(flux_type)           :: flux

  ! Additional arrays passed to radiation_scheme
  real(jprb), allocatable, dimension(:) :: ccn_land, ccn_sea, sin_latitude, longitude_rad, land_frac
  real(jprb), allocatable, dimension(:,:) :: pressure_fl, temperature_fl
  real(jprb), allocatable, dimension(:) :: flux_sw_direct_normal, flux_uv, flux_par, flux_par_clear, &
       &  emissivity_out
  real(jprb), allocatable, dimension(:,:) :: flux_diffuse_band, flux_direct_band

  ! Bespoke data types to set-up the blocked memory layout
  type(ifs_config_type)        :: ifs_config
  real(kind=jprb), allocatable :: zrgp(:,:,:) ! monolithic IFS data structure
  integer, allocatable         :: iseed(:,:) ! Seed for random number generator

  integer :: ncol, nlev         ! Number of columns and levels
  integer :: nproma             ! block size

  ! Name of file names specified on command line
  character(len=512) :: file_name
  integer            :: istatus ! Result of command_argument_count

  ! OpenMP functions
  integer, external :: omp_get_thread_num
  real(kind=jprd), external :: omp_get_wtime
  ! Start/stop time in seconds
  real(kind=jprd) :: tstart, tstop, t0

  ! For demonstration of get_sw_weights later on
  ! Ultraviolet weightings
  !integer    :: nweight_uv
  !integer    :: iband_uv(100)
  !real(jprb) :: weight_uv(100)
  ! Photosynthetically active radiation weightings
  !integer    :: nweight_par
  !integer    :: iband_par(100)
  !real(jprb) :: weight_par(100)

  ! Loop index for repeats (for benchmarking)
  integer :: jrepeat

  ! Loop index
  integer :: jrl, ibeg, iend, il, ib

  ! Are any variables out of bounds?
  logical :: is_out_of_bounds

!  integer    :: iband(20), nweights
!  real(jprb) :: weight(20)


  ! --------------------------------------------------------
  ! Section 2: Configure
  ! --------------------------------------------------------

  ! Check program called with correct number of arguments
  if (command_argument_count() < 3) then
    stop 'Usage: ecrad config.nam input_file.nc output_file.nc'
  end if

  ! Use namelist to configure the radiation calculation
  call get_command_argument(1, file_name, status=istatus)
  if (istatus /= 0) then
    stop 'Failed to read name of namelist file as string of length < 512'
  end if

  ! Read "radiation_driver" namelist into radiation driver config type
  call driver_config%read(file_name)
  nproma = driver_config%nblocksize

  if (driver_config%iverbose >= 2) then
    write(nulout,'(a)') '-------------------------- OFFLINE ECRAD RADIATION SCHEME --------------------------'
    write(nulout,'(a)') 'Copyright (C) 2014- ECMWF'
    write(nulout,'(a)') 'Contact: Robin Hogan (r.j.hogan@ecmwf.int)'
    write(nulout,'(a)') 'Floating-point precision: double'
  end if

  ! Albedo/emissivity intervals may be specified like this
  !call config%define_sw_albedo_intervals(6, &
  !     &  [0.25e-6_jprb, 0.44e-6_jprb, 0.69e-6_jprb, &
  !     &     1.19_jprb, 2.38e-6_jprb], [1,2,3,4,5,6], &
  !     &   do_nearest=.false.)
  !call config%define_lw_emiss_intervals(3, &
  !     &  [8.0e-6_jprb, 13.0e-6_jprb], [1,2,1], &
  !     &   do_nearest=.false.)

  ! If monochromatic aerosol properties are required, then the
  ! wavelengths can be specified (in metres) as follows - these can be
  ! whatever you like for the general aerosol optics, but must match
  ! the monochromatic values in the aerosol input file for the older
  ! aerosol optics
  !call config%set_aerosol_wavelength_mono( &
  !     &  [3.4e-07_jprb, 3.55e-07_jprb, 3.8e-07_jprb, 4.0e-07_jprb, 4.4e-07_jprb, &
  !     &   4.69e-07_jprb, 5.0e-07_jprb, 5.32e-07_jprb, 5.5e-07_jprb, 6.45e-07_jprb, &
  !     &   6.7e-07_jprb, 8.0e-07_jprb, 8.58e-07_jprb, 8.65e-07_jprb, 1.02e-06_jprb, &
  !     &   1.064e-06_jprb, 1.24e-06_jprb, 1.64e-06_jprb, 2.13e-06_jprb, 1.0e-05_jprb])

  call yradiation%rad_config%read_GPU(file_name=file_name, lacc=lacc)

  ! Setup aerosols
  if (yradiation%rad_config%use_aerosols) then
    yradiation%yrerad%naermacc = 1 ! MACC-derived aerosol climatology on a NMCLAT x NMCLON grid
  else
    yradiation%yrerad%naermacc = 0
  endif

  ! Setup the radiation scheme: load the coefficients for gas and
  ! cloud optics, currently from RRTMG
  call setup_radiation_scheme_GPU(yradiation, .true., file_name=file_name, lacc=lacc)
  ! Or call without specifying the namelist filename, in which case
  ! the default settings are from yoerad.F90
  !call setup_radiation_scheme(yradiation, .true.)

  ! Demonstration of how to get weights for UV and PAR fluxes
  !if (config%do_sw) then
  !  call config%get_sw_weights(0.2e-6_jprb, 0.4415e-6_jprb,&
  !       &  nweight_uv, iband_uv, weight_uv,&
  !       &  'ultraviolet')
  !  call config%get_sw_weights(0.4e-6_jprb, 0.7e-6_jprb,&
  !       &  nweight_par, iband_par, weight_par,&
  !       &  'photosynthetically active radiation, PAR')
  !end if

  ! --------------------------------------------------------
  ! Section 3: Read input data file
  ! --------------------------------------------------------

  ! Get NetCDF input file name
  call get_command_argument(2, file_name, status=istatus)
  if (istatus /= 0) then
    stop 'Failed to read name of input NetCDF file as string of length < 512'
  end if

  ! Open the file and configure the way it is read
  call file%open(trim(file_name), iverbose=driver_config%iverbose)

  ! Get NetCDF output file name
  call get_command_argument(3, file_name, status=istatus)
  if (istatus /= 0) then
    stop 'Failed to read name of output NetCDF file as string of length < 512'
  end if

  ! 2D arrays are assumed to be stored in the file with height varying
  ! more rapidly than column index. Specifying "true" here transposes
  ! all 2D arrays so that the column index varies fastest within the
  ! program.
  call file%transpose_matrices(.true.)

  ! Read input variables from NetCDF file, noting that cloud overlap
  ! and effective radius are ignored
  call read_input_GPU(file, yradiation%rad_config, driver_config, ncol, nlev, &
       &          single_level, thermodynamics, &
       &          gas, cloud, aerosol, lacc=lacc)

  ! Latitude is used for cloud overlap and ice effective radius
  if (file%exists('lat')) then
    call file%get('lat', sin_latitude)
    sin_latitude = sin(sin_latitude * Pi/180.0_jprb)
  else
    allocate(sin_latitude(ncol))
    sin_latitude = 0.0_jprb
  end if

  if (file%exists('lon')) then
    call file%get('lon', longitude_rad)
    longitude_rad = longitude_rad * Pi/180.0_jprb
  else
    allocate(longitude_rad(ncol))
    longitude_rad = 0.0_jprb
  end if

  ! Close input file
  call file%close()

  ! Convert gas units to mass-mixing ratio
  call gas%set_units_GPU(IMassMixingRatio, lacc=.false.)

  ! Compute seed from skin temperature residual
  !  single_level%iseed = int(1.0e9*(single_level%skin_temperature &
  !       &                            -int(single_level%skin_temperature)))

  ! Set first and last columns to process
  if (driver_config%iendcol < 1 .or. driver_config%iendcol > ncol) then
    driver_config%iendcol = ncol
  end if

  if (driver_config%istartcol > driver_config%iendcol) then
    write(nulout,'(a,i0,a,i0,a,i0,a)') '*** Error: requested column range (', &
         &  driver_config%istartcol, &
         &  ' to ', driver_config%iendcol, ') is out of the range in the data (1 to ', &
         &  ncol, ')'
    stop 1
  end if

  ! --------------------------------------------------------
  ! Section 4: Call radiation scheme
  ! --------------------------------------------------------

  ! Compute saturation with respect to liquid (needed for aerosol
  ! hydration) call
  !  call thermodynamics%calc_saturation_wrt_liquid(driver_config%istartcol,driver_config%iendcol)

  ! Check inputs are within physical bounds, printing message if not
  is_out_of_bounds =     gas%out_of_physical_bounds_GPU(driver_config%istartcol, driver_config%iendcol, &
       &                                            driver_config%do_correct_unphysical_inputs, lacc) &
       & .or.   single_level%out_of_physical_bounds_GPU(driver_config%istartcol, driver_config%iendcol, &
       &                                            driver_config%do_correct_unphysical_inputs, lacc) &
       & .or. thermodynamics%out_of_physical_bounds_GPU(driver_config%istartcol, driver_config%iendcol, &
       &                                            driver_config%do_correct_unphysical_inputs, lacc) &
       & .or.          cloud%out_of_physical_bounds_GPU(driver_config%istartcol, driver_config%iendcol, &
       &                                            driver_config%do_correct_unphysical_inputs, lacc) &
       & .or.        aerosol%out_of_physical_bounds_GPU(driver_config%istartcol, driver_config%iendcol, &
       &                                            driver_config%do_correct_unphysical_inputs, lacc)

  ! Allocate memory for the flux profiles, which may include arrays
  ! of dimension n_bands_sw/n_bands_lw, so must be called after
  ! setup_radiation
  call flux%allocate_GPU(yradiation%rad_config, 1, ncol, nlev, lacc=lacc)

  ! set relevant fluxes to zero
  flux%lw_up(:,:) = 0._jprb
  flux%lw_dn(:,:) = 0._jprb
  flux%sw_up(:,:) = 0._jprb
  flux%sw_dn(:,:) = 0._jprb
  flux%sw_dn_direct(:,:) = 0._jprb
  flux%lw_up_clear(:,:) = 0._jprb
  flux%lw_dn_clear(:,:) = 0._jprb
  flux%sw_up_clear(:,:) = 0._jprb
  flux%sw_dn_clear(:,:) = 0._jprb
  flux%sw_dn_direct_clear(:,:) = 0._jprb

  flux%lw_dn_surf_canopy(:,:) = 0._jprb
  flux%sw_dn_diffuse_surf_canopy(:,:) = 0._jprb
  flux%sw_dn_direct_surf_canopy(:,:) = 0._jprb
  flux%lw_derivatives(:,:) = 0._jprb

  ! Allocate memory for additional arrays
  allocate(ccn_land(ncol))
  allocate(ccn_sea(ncol))
  allocate(land_frac(ncol))
  allocate(pressure_fl(ncol,nlev))
  allocate(temperature_fl(ncol,nlev))
  allocate(flux_sw_direct_normal(ncol))
  allocate(flux_uv(ncol))
  allocate(flux_par(ncol))
  allocate(flux_par_clear(ncol))
  allocate(emissivity_out(ncol))
  allocate(flux_diffuse_band(ncol,yradiation%yrerad%nsw))
  allocate(flux_direct_band(ncol,yradiation%yrerad%nsw))

  ccn_land = yradiation%yrerad%rccnlnd
  ccn_sea = yradiation%yrerad%rccnsea
  pressure_fl = 0.5_jprb * (thermodynamics%pressure_hl(:,1:nlev)+thermodynamics%pressure_hl(:,2:nlev+1))
  temperature_fl = 0.5_jprb * (thermodynamics%temperature_hl(:,1:nlev)+thermodynamics%temperature_hl(:,2:nlev+1))

  ! --------------------------------------------------------
  ! Section 4a: Reshuffle into blocked memory layout
  ! --------------------------------------------------------

  call ifs_setup_indices_GPU(driver_config, ifs_config, yradiation, nlev, lacc=lacc)
  call ifs_copy_inputs_to_blocked_GPU(driver_config, ifs_config, yradiation,&
        & ncol, nlev, single_level, thermodynamics, gas, cloud, aerosol,&
        & sin_latitude, longitude_rad, land_frac, pressure_fl, temperature_fl,&
        & zrgp &
        &, iseed=iseed, lacc=lacc &
        & )

  ! --------------------------------------------------------
  ! Section 4b: Call radiation_scheme with blocked memory data
  ! --------------------------------------------------------

  if (driver_config%iverbose >= 2) then
    write(nulout,'(a)')  'Performing radiative transfer calculations'
  end if

  t0 = omp_get_wtime()

  ! Option of repeating calculation multiple time for more accurate
  ! profiling
  tstart = omp_get_wtime()
  do jrepeat = 1,driver_config%nrepeat
     
#ifdef __PGI
     call nvtxStartRange("ecrad_it")
#endif

!    if (driver_config%do_parallel) then
      ! Run radiation scheme over blocks of columns in parallel

      !$OMP PARALLEL DO SCHEDULE(DYNAMIC,1)&
      !$OMP&PRIVATE(JRL,IBEG,IEND,IL,IB)
      do jrl=1,ncol,nproma
        ibeg=jrl
        iend=min(ibeg+nproma-1,ncol)
        il=iend-ibeg+1
        ib=(jrl-1)/nproma+1

        if (driver_config%iverbose >= 3) then
          write(nulout,'(a,i0,a,i0,a,i0)')  'Thread ', omp_get_thread_num(), &
               &  ' processing columns ', ibeg, '-', iend
        end if

        !$acc data create(zrgp(:,:,ib)) &
        !$acc&     copyin(iseed(:,ib)) IF(lacc) 


        !$acc update device(zrgp(1:il,:,ib)) IF(lacc)
        ! Call the ECRAD radiation scheme
        call radiation_scheme_GPU &
             & (yradiation, &
             &  1, il, nproma, &                       ! startcol, endcol, ncol
             &  nlev, size(aerosol%mixing_ratio,3), &    ! nlev, naerosols
             &  single_level%solar_irradiance, &                               ! solar_irrad
             ! array inputs
             &  zrgp(:,ifs_config%iamu0,ib), zrgp(:,ifs_config%its,ib), &    ! mu0, skintemp
             &  zrgp(:,ifs_config%iald,ib) , zrgp(:,ifs_config%ialp,ib), &    ! albedo_dif, albedo_dir
             &  zrgp(:,ifs_config%iemiss,ib), &                   ! spectral emissivity
             &  zrgp(:,ifs_config%iccnl,ib), zrgp(:,ifs_config%iccno,ib) ,&  ! CCN concentration, land and sea
             &  zrgp(:,ifs_config%igelam,ib),zrgp(:,ifs_config%igemu,ib), &  ! longitude, sine of latitude
             &  zrgp(:,ifs_config%islm,ib), &                     ! land sea mask
             &  zrgp(:,ifs_config%ipr,ib),   zrgp(:,ifs_config%iti,ib),  &   ! full level pressure and temperature
             &  zrgp(:,ifs_config%iaprs,ib), zrgp(:,ifs_config%ihti,ib), &   ! half-level pressure and temperature
             &  zrgp(:,ifs_config%iwv,ib),   zrgp(:,ifs_config%iico2,ib), &
             &  zrgp(:,ifs_config%iich4,ib), zrgp(:,ifs_config%iin2o,ib), &
             &  zrgp(:,ifs_config%ino2,ib),  zrgp(:,ifs_config%ic11,ib), &
             &  zrgp(:,ifs_config%ic12,ib),  zrgp(:,ifs_config%ic22,ib), &
             &  zrgp(:,ifs_config%icl4,ib),  zrgp(:,ifs_config%ioz,ib), &
             &  zrgp(:,ifs_config%iclc,ib),  zrgp(:,ifs_config%ilwa,ib), &
             &  zrgp(:,ifs_config%iiwa,ib),  zrgp(:,ifs_config%irwa,ib), &
             &  zrgp(:,ifs_config%iswa,ib), &
             &  zrgp(:,ifs_config%iaer,ib),  zrgp(:,ifs_config%iaero,ib), &
             ! flux outputs
             &  zrgp(:,ifs_config%ifrso,ib), zrgp(:,ifs_config%ifrth,ib), &
             &  zrgp(:,ifs_config%iswfc,ib), zrgp(:,ifs_config%ilwfc,ib),&
             &  zrgp(:,ifs_config%ifrsod,ib),zrgp(:,ifs_config%ifrted,ib), &
             &  zrgp(:,ifs_config%ifrsodc,ib),zrgp(:,ifs_config%ifrtedc,ib),&
             &  zrgp(:,ifs_config%ifdir,ib), zrgp(:,ifs_config%icdir,ib), &
             &  zrgp(:,ifs_config%isudu,ib), &
             &  zrgp(:,ifs_config%iuvdf,ib), zrgp(:,ifs_config%iparf,ib), &
             &  zrgp(:,ifs_config%iparcf,ib),zrgp(:,ifs_config%itincf,ib), &
             &  zrgp(:,ifs_config%iemit,ib) ,zrgp(:,ifs_config%ilwderivative,ib), &
             &  zrgp(:,ifs_config%iswdiffuseband,ib), zrgp(:,ifs_config%iswdirectband,ib)&
            ! To validate results against standalone ecrad, we overwrite effective
            ! radii, cloud overlap and seed with input values
             &  ,pre_liq=zrgp(:,ifs_config%ire_liq,ib), &
             &  pre_ice=zrgp(:,ifs_config%ire_ice,ib), &
             &  pcloud_overlap=zrgp(:,ifs_config%ioverlap,ib), &
             &  iseed=iseed(:,ib), lacc=lacc &
             & )
          !$acc update host(zrgp(1:il,:,ib)) IF(lacc)
          !$acc end data
      end do
      !$OMP END PARALLEL DO

!    else
      ! Run radiation scheme serially
!      if (driver_config%iverbose >= 3) then
!        write(nulout,'(a,i0,a)')  'Processing ', ncol, ' columns'
!      end if

      ! Call the ECRAD radiation scheme
!      call radiation_scheme(ncol, nlev, driver_config%istartcol, driver_config%iendcol, &
!           &  config, single_level, thermodynamics, gas, cloud, aerosol, flux)

!    end if
     
#ifdef __PGI
     call nvtxEndRange
#endif

  end do

  tstop = omp_get_wtime()
  write(nulout, '(a,g12.5,a)') 'Time elapsed in radiative transfer: ', tstop-tstart, ' seconds'

  tstop = omp_get_wtime()
  write(nulout, '(a,g12.5,a)') 'Time elapsed in radiative transfer including data offload: ', tstop-t0, ' seconds'

  ! --------------------------------------------------------
  ! Section 4c: Copy fluxes from blocked memory data
  ! --------------------------------------------------------

  call ifs_copy_fluxes_from_blocked_GPU(driver_config, ifs_config, yradiation, ncol, nlev,&
          & zrgp, flux, flux_sw_direct_normal, flux_uv, flux_par, flux_par_clear, &
          & emissivity_out, flux_diffuse_band, flux_direct_band, lacc=lacc)

  ! "up" fluxes are actually net fluxes at this point - we modify the
  ! upwelling flux so that net=dn-up, while the TOA and surface
  ! downwelling fluxes are correct.
  flux%sw_up = -flux%sw_up
  flux%sw_up(:,1) = flux%sw_up(:,1)+flux%sw_dn(:,1)
  flux%sw_up(:,nlev+1) = flux%sw_up(:,nlev+1)+flux%sw_dn(:,nlev+1)

  flux%lw_up = -flux%lw_up
  flux%lw_up(:,1) = flux%lw_up(:,1)+flux%lw_dn(:,1)
  flux%lw_up(:,nlev+1) = flux%lw_up(:,nlev+1)+flux%lw_dn(:,nlev+1)

  flux%sw_up_clear = -flux%sw_up_clear
  flux%sw_up_clear(:,1) = flux%sw_up_clear(:,1)+flux%sw_dn_clear(:,1)
  flux%sw_up_clear(:,nlev+1) = flux%sw_up_clear(:,nlev+1)+flux%sw_dn_clear(:,nlev+1)

  flux%lw_up_clear = -flux%lw_up_clear
  flux%lw_up_clear(:,1) = flux%lw_up_clear(:,1)+flux%lw_dn_clear(:,1)
  flux%lw_up_clear(:,nlev+1) = flux%lw_up_clear(:,nlev+1)+flux%lw_dn_clear(:,nlev+1)

  ! --------------------------------------------------------
  ! Section 5: Check and save output
  ! --------------------------------------------------------

  ! This is unreliable because only the net fluxes are valid:
  !is_out_of_bounds = flux%out_of_physical_bounds(driver_config%istartcol, driver_config%iendcol)

  ! Store the fluxes in the output file
  yradiation%rad_config%do_surface_sw_spectral_flux = .false.
  yradiation%rad_config%do_canopy_fluxes_sw = .false.
  yradiation%rad_config%do_canopy_fluxes_lw = .false.

  call save_net_fluxes_GPU(file_name, yradiation%rad_config, thermodynamics, flux, &
       &   iverbose=driver_config%iverbose, is_hdf5_file=driver_config%do_write_hdf5, &
       &   experiment_name=driver_config%experiment_name, &
       &   is_double_precision=driver_config%do_write_double_precision, lacc=lacc)

  if (driver_config%iverbose >= 2) then
    write(nulout,'(a)') '------------------------------------------------------------------------------------'
  end if

END SUBROUTINE

SUBROUTINE ECRAD_IFS_DRIVER_CPU

  ! --------------------------------------------------------
  ! Section 1: Declarations
  ! --------------------------------------------------------
  use parkind1 ! Working/double precision

  use radiation_io
  use radiation_single_level
  use radiation_thermodynamics
  use radiation_gas
  use radiation_cloud
  use radiation_aerosol
  use radiation_flux
  use radiation_save
  use radiation_setup
  use radiation_constants
  use ecrad_driver_config
  use ecrad_driver_read_input
  use easy_netcdf
  use ifs_blocking
  
#ifdef __PGI
  use nvtx
#endif

  implicit none

#include "radiation_scheme.intfb.h"

  ! The NetCDF file containing the input profiles
  type(netcdf_file)         :: file

  ! Configuration for the radiation scheme, IFS style
  type(tradiation)          :: yradiation

  ! Derived types for the inputs to the radiation scheme
  type(single_level_type)   :: single_level
  type(thermodynamics_type) :: thermodynamics
  type(gas_type)            :: gas
  type(cloud_type)          :: cloud
  type(aerosol_type)        :: aerosol

  ! Configuration specific to this driver
  type(driver_config_type)  :: driver_config

  ! Derived type containing outputs from the radiation scheme
  type(flux_type)           :: flux

  ! Additional arrays passed to radiation_scheme
  real(jprb), allocatable, dimension(:) :: ccn_land, ccn_sea, sin_latitude, longitude_rad, land_frac
  real(jprb), allocatable, dimension(:,:) :: pressure_fl, temperature_fl
  real(jprb), allocatable, dimension(:) :: flux_sw_direct_normal, flux_uv, flux_par, flux_par_clear, &
       &  emissivity_out
  real(jprb), allocatable, dimension(:,:) :: flux_diffuse_band, flux_direct_band

  ! Bespoke data types to set-up the blocked memory layout
  type(ifs_config_type)        :: ifs_config
  real(kind=jprb), allocatable :: zrgp(:,:,:) ! monolithic IFS data structure
  integer, allocatable         :: iseed(:,:) ! Seed for random number generator

  integer :: ncol, nlev         ! Number of columns and levels
  integer :: nproma             ! block size

  ! Name of file names specified on command line
  character(len=512) :: file_name
  integer            :: istatus ! Result of command_argument_count

  ! OpenMP functions
  integer, external :: omp_get_thread_num
  real(kind=jprd), external :: omp_get_wtime
  ! Start/stop time in seconds
  real(kind=jprd) :: tstart, tstop, t0

  ! For demonstration of get_sw_weights later on
  ! Ultraviolet weightings
  !integer    :: nweight_uv
  !integer    :: iband_uv(100)
  !real(jprb) :: weight_uv(100)
  ! Photosynthetically active radiation weightings
  !integer    :: nweight_par
  !integer    :: iband_par(100)
  !real(jprb) :: weight_par(100)

  ! Loop index for repeats (for benchmarking)
  integer :: jrepeat

  ! Loop index
  integer :: jrl, ibeg, iend, il, ib

  ! Are any variables out of bounds?
  logical :: is_out_of_bounds

!  integer    :: iband(20), nweights
!  real(jprb) :: weight(20)


  ! --------------------------------------------------------
  ! Section 2: Configure
  ! --------------------------------------------------------

  ! Check program called with correct number of arguments
  if (command_argument_count() < 3) then
    stop 'Usage: ecrad config.nam input_file.nc output_file.nc'
  end if

  ! Use namelist to configure the radiation calculation
  call get_command_argument(1, file_name, status=istatus)
  if (istatus /= 0) then
    stop 'Failed to read name of namelist file as string of length < 512'
  end if

  ! Read "radiation_driver" namelist into radiation driver config type
  call driver_config%read(file_name)
  nproma = driver_config%nblocksize

  if (driver_config%iverbose >= 2) then
    write(nulout,'(a)') '-------------------------- OFFLINE ECRAD RADIATION SCHEME --------------------------'
    write(nulout,'(a)') 'Copyright (C) 2014- ECMWF'
    write(nulout,'(a)') 'Contact: Robin Hogan (r.j.hogan@ecmwf.int)'
    write(nulout,'(a)') 'Floating-point precision: double'
  end if

  ! Albedo/emissivity intervals may be specified like this
  !call config%define_sw_albedo_intervals(6, &
  !     &  [0.25e-6_jprb, 0.44e-6_jprb, 0.69e-6_jprb, &
  !     &     1.19_jprb, 2.38e-6_jprb], [1,2,3,4,5,6], &
  !     &   do_nearest=.false.)
  !call config%define_lw_emiss_intervals(3, &
  !     &  [8.0e-6_jprb, 13.0e-6_jprb], [1,2,1], &
  !     &   do_nearest=.false.)

  ! If monochromatic aerosol properties are required, then the
  ! wavelengths can be specified (in metres) as follows - these can be
  ! whatever you like for the general aerosol optics, but must match
  ! the monochromatic values in the aerosol input file for the older
  ! aerosol optics
  !call config%set_aerosol_wavelength_mono( &
  !     &  [3.4e-07_jprb, 3.55e-07_jprb, 3.8e-07_jprb, 4.0e-07_jprb, 4.4e-07_jprb, &
  !     &   4.69e-07_jprb, 5.0e-07_jprb, 5.32e-07_jprb, 5.5e-07_jprb, 6.45e-07_jprb, &
  !     &   6.7e-07_jprb, 8.0e-07_jprb, 8.58e-07_jprb, 8.65e-07_jprb, 1.02e-06_jprb, &
  !     &   1.064e-06_jprb, 1.24e-06_jprb, 1.64e-06_jprb, 2.13e-06_jprb, 1.0e-05_jprb])

  call yradiation%rad_config%read_CPU(file_name=file_name)

  ! Setup aerosols
  if (yradiation%rad_config%use_aerosols) then
    yradiation%yrerad%naermacc = 1 ! MACC-derived aerosol climatology on a NMCLAT x NMCLON grid
  else
    yradiation%yrerad%naermacc = 0
  endif

  ! Setup the radiation scheme: load the coefficients for gas and
  ! cloud optics, currently from RRTMG
  call setup_radiation_scheme_CPU(yradiation, .true., file_name=file_name)
  ! Or call without specifying the namelist filename, in which case
  ! the default settings are from yoerad.F90
  !call setup_radiation_scheme(yradiation, .true.)

  ! Demonstration of how to get weights for UV and PAR fluxes
  !if (config%do_sw) then
  !  call config%get_sw_weights(0.2e-6_jprb, 0.4415e-6_jprb,&
  !       &  nweight_uv, iband_uv, weight_uv,&
  !       &  'ultraviolet')
  !  call config%get_sw_weights(0.4e-6_jprb, 0.7e-6_jprb,&
  !       &  nweight_par, iband_par, weight_par,&
  !       &  'photosynthetically active radiation, PAR')
  !end if

  ! --------------------------------------------------------
  ! Section 3: Read input data file
  ! --------------------------------------------------------

  ! Get NetCDF input file name
  call get_command_argument(2, file_name, status=istatus)
  if (istatus /= 0) then
    stop 'Failed to read name of input NetCDF file as string of length < 512'
  end if

  ! Open the file and configure the way it is read
  call file%open(trim(file_name), iverbose=driver_config%iverbose)

  ! Get NetCDF output file name
  call get_command_argument(3, file_name, status=istatus)
  if (istatus /= 0) then
    stop 'Failed to read name of output NetCDF file as string of length < 512'
  end if

  ! 2D arrays are assumed to be stored in the file with height varying
  ! more rapidly than column index. Specifying "true" here transposes
  ! all 2D arrays so that the column index varies fastest within the
  ! program.
  call file%transpose_matrices(.true.)

  ! Read input variables from NetCDF file, noting that cloud overlap
  ! and effective radius are ignored
  call read_input_CPU(file, yradiation%rad_config, driver_config, ncol, nlev, &
       &          single_level, thermodynamics, &
       &          gas, cloud, aerosol)

  ! Latitude is used for cloud overlap and ice effective radius
  if (file%exists('lat')) then
    call file%get('lat', sin_latitude)
    sin_latitude = sin(sin_latitude * Pi/180.0_jprb)
  else
    allocate(sin_latitude(ncol))
    sin_latitude = 0.0_jprb
  end if

  if (file%exists('lon')) then
    call file%get('lon', longitude_rad)
    longitude_rad = longitude_rad * Pi/180.0_jprb
  else
    allocate(longitude_rad(ncol))
    longitude_rad = 0.0_jprb
  end if

  ! Close input file
  call file%close()

  ! Convert gas units to mass-mixing ratio
  call gas%set_units_CPU(IMassMixingRatio, lacc=.false.)

  ! Compute seed from skin temperature residual
  !  single_level%iseed = int(1.0e9*(single_level%skin_temperature &
  !       &                            -int(single_level%skin_temperature)))

  ! Set first and last columns to process
  if (driver_config%iendcol < 1 .or. driver_config%iendcol > ncol) then
    driver_config%iendcol = ncol
  end if

  if (driver_config%istartcol > driver_config%iendcol) then
    write(nulout,'(a,i0,a,i0,a,i0,a)') '*** Error: requested column range (', &
         &  driver_config%istartcol, &
         &  ' to ', driver_config%iendcol, ') is out of the range in the data (1 to ', &
         &  ncol, ')'
    stop 1
  end if

  ! --------------------------------------------------------
  ! Section 4: Call radiation scheme
  ! --------------------------------------------------------

  ! Compute saturation with respect to liquid (needed for aerosol
  ! hydration) call
  !  call thermodynamics%calc_saturation_wrt_liquid(driver_config%istartcol,driver_config%iendcol)

  ! Check inputs are within physical bounds, printing message if not
  is_out_of_bounds =     gas%out_of_physical_bounds_CPU(driver_config%istartcol, driver_config%iendcol, &
       &                                            driver_config%do_correct_unphysical_inputs) &
       & .or.   single_level%out_of_physical_bounds_CPU(driver_config%istartcol, driver_config%iendcol, &
       &                                            driver_config%do_correct_unphysical_inputs) &
       & .or. thermodynamics%out_of_physical_bounds_CPU(driver_config%istartcol, driver_config%iendcol, &
       &                                            driver_config%do_correct_unphysical_inputs) &
       & .or.          cloud%out_of_physical_bounds_CPU(driver_config%istartcol, driver_config%iendcol, &
       &                                            driver_config%do_correct_unphysical_inputs) &
       & .or.        aerosol%out_of_physical_bounds_CPU(driver_config%istartcol, driver_config%iendcol, &
       &                                            driver_config%do_correct_unphysical_inputs)

  ! Allocate memory for the flux profiles, which may include arrays
  ! of dimension n_bands_sw/n_bands_lw, so must be called after
  ! setup_radiation
  call flux%allocate_CPU(yradiation%rad_config, 1, ncol, nlev)

  ! set relevant fluxes to zero
  flux%lw_up(:,:) = 0._jprb
  flux%lw_dn(:,:) = 0._jprb
  flux%sw_up(:,:) = 0._jprb
  flux%sw_dn(:,:) = 0._jprb
  flux%sw_dn_direct(:,:) = 0._jprb
  flux%lw_up_clear(:,:) = 0._jprb
  flux%lw_dn_clear(:,:) = 0._jprb
  flux%sw_up_clear(:,:) = 0._jprb
  flux%sw_dn_clear(:,:) = 0._jprb
  flux%sw_dn_direct_clear(:,:) = 0._jprb

  flux%lw_dn_surf_canopy(:,:) = 0._jprb
  flux%sw_dn_diffuse_surf_canopy(:,:) = 0._jprb
  flux%sw_dn_direct_surf_canopy(:,:) = 0._jprb
  flux%lw_derivatives(:,:) = 0._jprb

  ! Allocate memory for additional arrays
  allocate(ccn_land(ncol))
  allocate(ccn_sea(ncol))
  allocate(land_frac(ncol))
  allocate(pressure_fl(ncol,nlev))
  allocate(temperature_fl(ncol,nlev))
  allocate(flux_sw_direct_normal(ncol))
  allocate(flux_uv(ncol))
  allocate(flux_par(ncol))
  allocate(flux_par_clear(ncol))
  allocate(emissivity_out(ncol))
  allocate(flux_diffuse_band(ncol,yradiation%yrerad%nsw))
  allocate(flux_direct_band(ncol,yradiation%yrerad%nsw))

  ccn_land = yradiation%yrerad%rccnlnd
  ccn_sea = yradiation%yrerad%rccnsea
  pressure_fl = 0.5_jprb * (thermodynamics%pressure_hl(:,1:nlev)+thermodynamics%pressure_hl(:,2:nlev+1))
  temperature_fl = 0.5_jprb * (thermodynamics%temperature_hl(:,1:nlev)+thermodynamics%temperature_hl(:,2:nlev+1))

  ! --------------------------------------------------------
  ! Section 4a: Reshuffle into blocked memory layout
  ! --------------------------------------------------------

  call ifs_setup_indices_CPU(driver_config, ifs_config, yradiation, nlev)
  call ifs_copy_inputs_to_blocked_CPU(driver_config, ifs_config, yradiation,&
        & ncol, nlev, single_level, thermodynamics, gas, cloud, aerosol,&
        & sin_latitude, longitude_rad, land_frac, pressure_fl, temperature_fl,&
        & zrgp &
        &, iseed=iseed &
        & )

  ! --------------------------------------------------------
  ! Section 4b: Call radiation_scheme with blocked memory data
  ! --------------------------------------------------------

  if (driver_config%iverbose >= 2) then
    write(nulout,'(a)')  'Performing radiative transfer calculations'
  end if

  t0 = omp_get_wtime()

  ! Option of repeating calculation multiple time for more accurate
  ! profiling
  tstart = omp_get_wtime()
  do jrepeat = 1,driver_config%nrepeat
     
#ifdef __PGI
     call nvtxStartRange("ecrad_it")
#endif

!    if (driver_config%do_parallel) then
      ! Run radiation scheme over blocks of columns in parallel

      !$OMP PARALLEL DO SCHEDULE(DYNAMIC,1)&
      !$OMP&PRIVATE(JRL,IBEG,IEND,IL,IB)
      do jrl=1,ncol,nproma
        ibeg=jrl
        iend=min(ibeg+nproma-1,ncol)
        il=iend-ibeg+1
        ib=(jrl-1)/nproma+1

        if (driver_config%iverbose >= 3) then
          write(nulout,'(a,i0,a,i0,a,i0)')  'Thread ', omp_get_thread_num(), &
               &  ' processing columns ', ibeg, '-', iend
        end if

        !$acc data create(zrgp(:,:,ib)) &
        !$acc&     copyin(iseed(:,ib)) 


        !$acc update device(zrgp(1:il,:,ib))
        ! Call the ECRAD radiation scheme
        call radiation_scheme_CPU &
             & (yradiation, &
             &  1, il, nproma, &                       ! startcol, endcol, ncol
             &  nlev, size(aerosol%mixing_ratio,3), &    ! nlev, naerosols
             &  single_level%solar_irradiance, &                               ! solar_irrad
             ! array inputs
             &  zrgp(:,ifs_config%iamu0,ib), zrgp(:,ifs_config%its,ib), &    ! mu0, skintemp
             &  zrgp(:,ifs_config%iald,ib) , zrgp(:,ifs_config%ialp,ib), &    ! albedo_dif, albedo_dir
             &  zrgp(:,ifs_config%iemiss,ib), &                   ! spectral emissivity
             &  zrgp(:,ifs_config%iccnl,ib), zrgp(:,ifs_config%iccno,ib) ,&  ! CCN concentration, land and sea
             &  zrgp(:,ifs_config%igelam,ib),zrgp(:,ifs_config%igemu,ib), &  ! longitude, sine of latitude
             &  zrgp(:,ifs_config%islm,ib), &                     ! land sea mask
             &  zrgp(:,ifs_config%ipr,ib),   zrgp(:,ifs_config%iti,ib),  &   ! full level pressure and temperature
             &  zrgp(:,ifs_config%iaprs,ib), zrgp(:,ifs_config%ihti,ib), &   ! half-level pressure and temperature
             &  zrgp(:,ifs_config%iwv,ib),   zrgp(:,ifs_config%iico2,ib), &
             &  zrgp(:,ifs_config%iich4,ib), zrgp(:,ifs_config%iin2o,ib), &
             &  zrgp(:,ifs_config%ino2,ib),  zrgp(:,ifs_config%ic11,ib), &
             &  zrgp(:,ifs_config%ic12,ib),  zrgp(:,ifs_config%ic22,ib), &
             &  zrgp(:,ifs_config%icl4,ib),  zrgp(:,ifs_config%ioz,ib), &
             &  zrgp(:,ifs_config%iclc,ib),  zrgp(:,ifs_config%ilwa,ib), &
             &  zrgp(:,ifs_config%iiwa,ib),  zrgp(:,ifs_config%irwa,ib), &
             &  zrgp(:,ifs_config%iswa,ib), &
             &  zrgp(:,ifs_config%iaer,ib),  zrgp(:,ifs_config%iaero,ib), &
             ! flux outputs
             &  zrgp(:,ifs_config%ifrso,ib), zrgp(:,ifs_config%ifrth,ib), &
             &  zrgp(:,ifs_config%iswfc,ib), zrgp(:,ifs_config%ilwfc,ib),&
             &  zrgp(:,ifs_config%ifrsod,ib),zrgp(:,ifs_config%ifrted,ib), &
             &  zrgp(:,ifs_config%ifrsodc,ib),zrgp(:,ifs_config%ifrtedc,ib),&
             &  zrgp(:,ifs_config%ifdir,ib), zrgp(:,ifs_config%icdir,ib), &
             &  zrgp(:,ifs_config%isudu,ib), &
             &  zrgp(:,ifs_config%iuvdf,ib), zrgp(:,ifs_config%iparf,ib), &
             &  zrgp(:,ifs_config%iparcf,ib),zrgp(:,ifs_config%itincf,ib), &
             &  zrgp(:,ifs_config%iemit,ib) ,zrgp(:,ifs_config%ilwderivative,ib), &
             &  zrgp(:,ifs_config%iswdiffuseband,ib), zrgp(:,ifs_config%iswdirectband,ib)&
            ! To validate results against standalone ecrad, we overwrite effective
            ! radii, cloud overlap and seed with input values
             &  ,pre_liq=zrgp(:,ifs_config%ire_liq,ib), &
             &  pre_ice=zrgp(:,ifs_config%ire_ice,ib), &
             &  pcloud_overlap=zrgp(:,ifs_config%ioverlap,ib), &
             &  iseed=iseed(:,ib) &
             & )
          !$acc update host(zrgp(1:il,:,ib))
          !$acc end data
      end do
      !$OMP END PARALLEL DO

!    else
      ! Run radiation scheme serially
!      if (driver_config%iverbose >= 3) then
!        write(nulout,'(a,i0,a)')  'Processing ', ncol, ' columns'
!      end if

      ! Call the ECRAD radiation scheme
!      call radiation_scheme(ncol, nlev, driver_config%istartcol, driver_config%iendcol, &
!           &  config, single_level, thermodynamics, gas, cloud, aerosol, flux)

!    end if
     
#ifdef __PGI
     call nvtxEndRange
#endif

  end do

  tstop = omp_get_wtime()
  write(nulout, '(a,g12.5,a)') 'Time elapsed in radiative transfer: ', tstop-tstart, ' seconds'

  tstop = omp_get_wtime()
  write(nulout, '(a,g12.5,a)') 'Time elapsed in radiative transfer including data offload: ', tstop-t0, ' seconds'

  ! --------------------------------------------------------
  ! Section 4c: Copy fluxes from blocked memory data
  ! --------------------------------------------------------

  call ifs_copy_fluxes_from_blocked_CPU(driver_config, ifs_config, yradiation, ncol, nlev,&
          & zrgp, flux, flux_sw_direct_normal, flux_uv, flux_par, flux_par_clear, &
          & emissivity_out, flux_diffuse_band, flux_direct_band)

  ! "up" fluxes are actually net fluxes at this point - we modify the
  ! upwelling flux so that net=dn-up, while the TOA and surface
  ! downwelling fluxes are correct.
  flux%sw_up = -flux%sw_up
  flux%sw_up(:,1) = flux%sw_up(:,1)+flux%sw_dn(:,1)
  flux%sw_up(:,nlev+1) = flux%sw_up(:,nlev+1)+flux%sw_dn(:,nlev+1)

  flux%lw_up = -flux%lw_up
  flux%lw_up(:,1) = flux%lw_up(:,1)+flux%lw_dn(:,1)
  flux%lw_up(:,nlev+1) = flux%lw_up(:,nlev+1)+flux%lw_dn(:,nlev+1)

  flux%sw_up_clear = -flux%sw_up_clear
  flux%sw_up_clear(:,1) = flux%sw_up_clear(:,1)+flux%sw_dn_clear(:,1)
  flux%sw_up_clear(:,nlev+1) = flux%sw_up_clear(:,nlev+1)+flux%sw_dn_clear(:,nlev+1)

  flux%lw_up_clear = -flux%lw_up_clear
  flux%lw_up_clear(:,1) = flux%lw_up_clear(:,1)+flux%lw_dn_clear(:,1)
  flux%lw_up_clear(:,nlev+1) = flux%lw_up_clear(:,nlev+1)+flux%lw_dn_clear(:,nlev+1)

  ! --------------------------------------------------------
  ! Section 5: Check and save output
  ! --------------------------------------------------------

  ! This is unreliable because only the net fluxes are valid:
  !is_out_of_bounds = flux%out_of_physical_bounds(driver_config%istartcol, driver_config%iendcol)

  ! Store the fluxes in the output file
  yradiation%rad_config%do_surface_sw_spectral_flux = .false.
  yradiation%rad_config%do_canopy_fluxes_sw = .false.
  yradiation%rad_config%do_canopy_fluxes_lw = .false.

  call save_net_fluxes_CPU(file_name, yradiation%rad_config, thermodynamics, flux, &
       &   iverbose=driver_config%iverbose, is_hdf5_file=driver_config%do_write_hdf5, &
       &   experiment_name=driver_config%experiment_name, &
       &   is_double_precision=driver_config%do_write_double_precision)

  if (driver_config%iverbose >= 2) then
    write(nulout,'(a)') '------------------------------------------------------------------------------------'
  end if

END SUBROUTINE

END PROGRAM ECRAD_IFS_DRIVER
