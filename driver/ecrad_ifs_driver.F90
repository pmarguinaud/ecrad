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
  real(jprb), allocatable, dimension(:,:) :: pressure_fl, temperature_fl, zeros
  real(jprb), allocatable, dimension(:,:,:) :: tegen_aerosol
  real(jprb), allocatable, dimension(:) :: flux_sw_direct_normal, flux_uv, flux_par, flux_par_clear, &
       &  flux_incoming, emissivity_out
  real(jprb), allocatable, dimension(:,:) :: flux_diffuse_band, flux_direct_band
  real(jprb), allocatable, dimension(:,:) :: cloud_fraction, cloud_q_liq, cloud_q_ice

  integer :: ncol, nlev         ! Number of columns and levels
  integer :: istartcol, iendcol ! Range of columns to process

  ! Name of file names specified on command line
  character(len=512) :: file_name
  integer            :: istatus ! Result of command_argument_count

  ! For parallel processing of multiple blocks
  integer :: jblock, nblock ! Block loop index and number

  ! OpenMP functions
  integer, external :: omp_get_thread_num
  real(kind=jprd), external :: omp_get_wtime
  ! Start/stop time in seconds
  real(kind=jprd) :: tstart, tstop

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
  allocate(zeros(ncol,nlev))
  allocate(tegen_aerosol(ncol,6,nlev))
  allocate(flux_sw_direct_normal(ncol))
  allocate(flux_uv(ncol))
  allocate(flux_par(ncol))
  allocate(flux_par_clear(ncol))
  allocate(flux_incoming(ncol))
  allocate(emissivity_out(ncol))
  allocate(flux_diffuse_band(ncol,yradiation%yrerad%nsw))
  allocate(flux_direct_band(ncol,yradiation%yrerad%nsw))
  allocate(cloud_fraction(ncol,nlev))
  allocate(cloud_q_liq(ncol,nlev))
  allocate(cloud_q_ice(ncol,nlev))

  ccn_land = yradiation%yrerad%rccnlnd
  ccn_sea = yradiation%yrerad%rccnsea
  tegen_aerosol = 0.0_jprb
  pressure_fl = 0.5_jprb * (thermodynamics%pressure_hl(:,1:nlev)+thermodynamics%pressure_hl(:,2:nlev+1))
  temperature_fl = 0.5_jprb * (thermodynamics%temperature_hl(:,1:nlev)+thermodynamics%temperature_hl(:,2:nlev+1))
  zeros = 0.0_jprb ! Dummy snow/rain water mixing ratios

  if (yradiation%rad_config%do_clouds) then
    cloud_fraction = cloud%fraction
    cloud_q_liq = cloud%q_liq
    cloud_q_ice = cloud%q_ice
  else
    cloud_fraction = 0.0_jprb
    cloud_q_liq = 0.0_jprb
    cloud_q_ice = 0.0_jprb
  endif

  if (driver_config%iverbose >= 2) then
    write(nulout,'(a)')  'Performing radiative transfer calculations'
  end if

  ! Option of repeating calculation multiple time for more accurate
  ! profiling
  tstart = omp_get_wtime()
  do jrepeat = 1,driver_config%nrepeat

!    if (driver_config%do_parallel) then
      ! Run radiation scheme over blocks of columns in parallel

      ! Compute number of blocks to process
      nblock = (driver_config%iendcol - driver_config%istartcol &
           &  + driver_config%nblocksize) / driver_config%nblocksize

      !$OMP PARALLEL DO PRIVATE(istartcol, iendcol) SCHEDULE(RUNTIME)
      do jblock = 1, nblock
        ! Specify the range of columns to process.
        istartcol = (jblock-1) * driver_config%nblocksize &
             &    + driver_config%istartcol
        iendcol = min(istartcol + driver_config%nblocksize - 1, &
             &        driver_config%iendcol)

        if (driver_config%iverbose >= 3) then
          write(nulout,'(a,i0,a,i0,a,i0)')  'Thread ', omp_get_thread_num(), &
               &  ' processing columns ', istartcol, '-', iendcol
        end if

        ! Call the ECRAD radiation scheme; note that we are simply
        ! passing arrays in rather than ecRad structures, which are
        ! used here just for convenience
        call radiation_scheme_GPU(yradiation, istartcol, iendcol, ncol, nlev, size(aerosol%mixing_ratio,3), &
             &  single_level%solar_irradiance, single_level%cos_sza, single_level%skin_temperature, &
             &  single_level%sw_albedo, single_level%sw_albedo_direct, single_level%lw_emissivity, &
             &  ccn_land, ccn_sea, longitude_rad, sin_latitude, land_frac, pressure_fl, temperature_fl, &
             &  thermodynamics%pressure_hl, thermodynamics%temperature_hl, &
             &  gas%mixing_ratio(:,:,IH2O), gas%mixing_ratio(:,:,ICO2), &
             &  gas%mixing_ratio(:,:,ICH4), gas%mixing_ratio(:,:,IN2O), gas%mixing_ratio(:,:,INO2), &
             &  gas%mixing_ratio(:,:,ICFC11), gas%mixing_ratio(:,:,ICFC12), gas%mixing_ratio(:,:,IHCFC22), &
             &  gas%mixing_ratio(:,:,ICCl4), gas%mixing_ratio(:,:,IO3), cloud_fraction, cloud_q_liq, &
             &  cloud_q_ice, zeros, zeros, tegen_aerosol, aerosol%mixing_ratio, flux%sw_up, flux%lw_up, &
             &  flux%sw_up_clear, flux%lw_up_clear, flux%sw_dn(:,nlev+1), flux%lw_dn(:,nlev+1), &
             &  flux%sw_dn_clear(:,nlev+1), flux%lw_dn_clear(:,nlev+1), &
             &  flux%sw_dn_direct(:,nlev+1), flux%sw_dn_direct_clear(:,nlev+1), flux_sw_direct_normal, &
             &  flux_uv, flux_par, &
             &  flux_par_clear, flux%sw_dn(:,1), emissivity_out, flux%lw_derivatives, flux_diffuse_band, &
             &  flux_direct_band &
            ! To validate results against standalone ecrad, we overwrite effective
            ! radii, cloud overlap and seed with input values
             &  ,pre_liq=cloud%re_liq, pre_ice=cloud%re_ice, &
             &  pcloud_overlap=cloud%overlap_param, &
             &  iseed=single_level%iseed, lacc=lacc &
             & )
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

  end do

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

  tstop = omp_get_wtime()
  write(nulout, '(a,g12.5,a)') 'Time elapsed in radiative transfer: ', tstop-tstart, ' seconds'

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
  real(jprb), allocatable, dimension(:,:) :: pressure_fl, temperature_fl, zeros
  real(jprb), allocatable, dimension(:,:,:) :: tegen_aerosol
  real(jprb), allocatable, dimension(:) :: flux_sw_direct_normal, flux_uv, flux_par, flux_par_clear, &
       &  flux_incoming, emissivity_out
  real(jprb), allocatable, dimension(:,:) :: flux_diffuse_band, flux_direct_band
  real(jprb), allocatable, dimension(:,:) :: cloud_fraction, cloud_q_liq, cloud_q_ice

  integer :: ncol, nlev         ! Number of columns and levels
  integer :: istartcol, iendcol ! Range of columns to process

  ! Name of file names specified on command line
  character(len=512) :: file_name
  integer            :: istatus ! Result of command_argument_count

  ! For parallel processing of multiple blocks
  integer :: jblock, nblock ! Block loop index and number

  ! OpenMP functions
  integer, external :: omp_get_thread_num
  real(kind=jprd), external :: omp_get_wtime
  ! Start/stop time in seconds
  real(kind=jprd) :: tstart, tstop

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
  allocate(zeros(ncol,nlev))
  allocate(tegen_aerosol(ncol,6,nlev))
  allocate(flux_sw_direct_normal(ncol))
  allocate(flux_uv(ncol))
  allocate(flux_par(ncol))
  allocate(flux_par_clear(ncol))
  allocate(flux_incoming(ncol))
  allocate(emissivity_out(ncol))
  allocate(flux_diffuse_band(ncol,yradiation%yrerad%nsw))
  allocate(flux_direct_band(ncol,yradiation%yrerad%nsw))
  allocate(cloud_fraction(ncol,nlev))
  allocate(cloud_q_liq(ncol,nlev))
  allocate(cloud_q_ice(ncol,nlev))

  ccn_land = yradiation%yrerad%rccnlnd
  ccn_sea = yradiation%yrerad%rccnsea
  tegen_aerosol = 0.0_jprb
  pressure_fl = 0.5_jprb * (thermodynamics%pressure_hl(:,1:nlev)+thermodynamics%pressure_hl(:,2:nlev+1))
  temperature_fl = 0.5_jprb * (thermodynamics%temperature_hl(:,1:nlev)+thermodynamics%temperature_hl(:,2:nlev+1))
  zeros = 0.0_jprb ! Dummy snow/rain water mixing ratios

  if (yradiation%rad_config%do_clouds) then
    cloud_fraction = cloud%fraction
    cloud_q_liq = cloud%q_liq
    cloud_q_ice = cloud%q_ice
  else
    cloud_fraction = 0.0_jprb
    cloud_q_liq = 0.0_jprb
    cloud_q_ice = 0.0_jprb
  endif

  if (driver_config%iverbose >= 2) then
    write(nulout,'(a)')  'Performing radiative transfer calculations'
  end if

  ! Option of repeating calculation multiple time for more accurate
  ! profiling
  tstart = omp_get_wtime()
  do jrepeat = 1,driver_config%nrepeat

!    if (driver_config%do_parallel) then
      ! Run radiation scheme over blocks of columns in parallel

      ! Compute number of blocks to process
      nblock = (driver_config%iendcol - driver_config%istartcol &
           &  + driver_config%nblocksize) / driver_config%nblocksize

      !$OMP PARALLEL DO PRIVATE(istartcol, iendcol) SCHEDULE(RUNTIME)
      do jblock = 1, nblock
        ! Specify the range of columns to process.
        istartcol = (jblock-1) * driver_config%nblocksize &
             &    + driver_config%istartcol
        iendcol = min(istartcol + driver_config%nblocksize - 1, &
             &        driver_config%iendcol)

        if (driver_config%iverbose >= 3) then
          write(nulout,'(a,i0,a,i0,a,i0)')  'Thread ', omp_get_thread_num(), &
               &  ' processing columns ', istartcol, '-', iendcol
        end if

        ! Call the ECRAD radiation scheme; note that we are simply
        ! passing arrays in rather than ecRad structures, which are
        ! used here just for convenience
        call radiation_scheme_CPU(yradiation, istartcol, iendcol, ncol, nlev, size(aerosol%mixing_ratio,3), &
             &  single_level%solar_irradiance, single_level%cos_sza, single_level%skin_temperature, &
             &  single_level%sw_albedo, single_level%sw_albedo_direct, single_level%lw_emissivity, &
             &  ccn_land, ccn_sea, longitude_rad, sin_latitude, land_frac, pressure_fl, temperature_fl, &
             &  thermodynamics%pressure_hl, thermodynamics%temperature_hl, &
             &  gas%mixing_ratio(:,:,IH2O), gas%mixing_ratio(:,:,ICO2), &
             &  gas%mixing_ratio(:,:,ICH4), gas%mixing_ratio(:,:,IN2O), gas%mixing_ratio(:,:,INO2), &
             &  gas%mixing_ratio(:,:,ICFC11), gas%mixing_ratio(:,:,ICFC12), gas%mixing_ratio(:,:,IHCFC22), &
             &  gas%mixing_ratio(:,:,ICCl4), gas%mixing_ratio(:,:,IO3), cloud_fraction, cloud_q_liq, &
             &  cloud_q_ice, zeros, zeros, tegen_aerosol, aerosol%mixing_ratio, flux%sw_up, flux%lw_up, &
             &  flux%sw_up_clear, flux%lw_up_clear, flux%sw_dn(:,nlev+1), flux%lw_dn(:,nlev+1), &
             &  flux%sw_dn_clear(:,nlev+1), flux%lw_dn_clear(:,nlev+1), &
             &  flux%sw_dn_direct(:,nlev+1), flux%sw_dn_direct_clear(:,nlev+1), flux_sw_direct_normal, &
             &  flux_uv, flux_par, &
             &  flux_par_clear, flux%sw_dn(:,1), emissivity_out, flux%lw_derivatives, flux_diffuse_band, &
             &  flux_direct_band &
            ! To validate results against standalone ecrad, we overwrite effective
            ! radii, cloud overlap and seed with input values
             &  ,pre_liq=cloud%re_liq, pre_ice=cloud%re_ice, &
             &  pcloud_overlap=cloud%overlap_param, &
             &  iseed=single_level%iseed &
             & )
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

  end do

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

  tstop = omp_get_wtime()
  write(nulout, '(a,g12.5,a)') 'Time elapsed in radiative transfer: ', tstop-tstart, ' seconds'

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
