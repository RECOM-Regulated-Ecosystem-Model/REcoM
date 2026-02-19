module recom_atbox_module
    interface
        subroutine recom_atbox(MPI_COMM_FESOM, myDim_nod2D, eDim_nod2D, ulevels_nod2D, areasvol, dt)
            use recom_declarations, only: wp

            integer,       intent(in) :: MPI_COMM_FESOM, myDim_nod2D, eDim_nod2D
            real(kind=WP), intent(in)                 :: dt
            integer,       intent(in), dimension(:)   :: ulevels_nod2D
            real(kind=WP), intent(in), dimension(:,:) :: areasvol
        end subroutine
    end interface
end module recom_atbox_module

    subroutine recom_atbox(MPI_COMM_FESOM, myDim_nod2D, eDim_nod2D, ulevels_nod2D, areasvol, dt)
!     Simple 0-d box model to calculate the temporal evolution of atmospheric CO2.
!     Initially the box model was part of module recom_ciso. Now it can be run also
!     without carbon isotopes (ciso==.false.)
!     mbutzin, 2021-07-08
      use REcoM_GloVar
      use recom_config
      use recom_ciso
      use recom_extra, only: integrate_nod_2d_recom
      use recom_declarations, only: wp

      implicit none

      integer,       intent(in) :: MPI_COMM_FESOM, myDim_nod2D, eDim_nod2D
      real(kind=WP), intent(in)                 :: dt
      integer,       intent(in), dimension(:)   :: ulevels_nod2D
      real(kind=WP), intent(in), dimension(:,:) :: areasvol

      integer                  :: n, elem, elnodes(3),n1
      real(kind=WP), parameter :: mol_allatm = 1.7726e20  ! atmospheric inventory of all compounds (mol)
      real(kind=WP)            :: total_co2flux,    &     ! (mol / s)
                                  total_co2flux_13, &     ! (mol / s) carbon-13
                                  total_co2flux_14        ! (mol / s) radiocarbon

!     Globally integrated air-sea CO2 flux (mol / s)
      total_co2flux = 0.
      call integrate_nod_2d_recom(0.001 * GloCO2flux_seaicemask, total_CO2flux, MPI_COMM_FESOM, &
                                  myDim_nod2D, eDim_nod2D, ulevels_nod2D, areasvol)

!     Atmospheric carbon budget (mol)
!     mass of the dry atmosphere = 5.1352e18 kg (Trenberth & Smith 2005, doi:10.1175/JCLI-3299.1)
!     mean density of air = 0.02897 kg / mol (https://nssdc.gsfc.nasa.gov/planetary/factsheet/earthfact.html)
!     => total molecular inventory of the dry atmosphere: moles_atm = 1.7726e20 mol == constant.
!     mol_co2atm    = mol_co2atm    - total_co2flux    * dt
!     Atmospheric mixing ratios in ppm
!     x_co2atm(1)    = mol_co2atm    / mol_allatm * 1.e6 ! ppm
      x_co2atm(1) = x_co2atm(1) - total_co2flux / mol_allatm * dt * 1.e6
      x_co2atm    = x_co2atm(1)

      if (ciso) then
!       Consider 13CO2 (and maybe also 14CO2)

!       Globally integrated air-sea 13CO2 flux (mol / s)
        total_co2flux_13 = 0.
        call integrate_nod_2d_recom(0.001 * GloCO2flux_seaicemask_13, total_co2flux_13,     &
                                    MPI_COMM_FESOM, myDim_nod2D, eDim_nod2D, ulevels_nod2D, &
                                    areasvol)

!       Atmospheric carbon-13 budget (mol)
!       mol_co2atm_13 = mol_co2atm_13 - total_co2flux_13 * dt
!       Budget in terms of the 13C / 12C volume mixing ratio
!       x_co2atm_13(1) = mol_co2atm_13 / mol_allatm * 1.e6
        x_co2atm_13(1) = x_co2atm_13(1) - total_co2flux_13 / mol_allatm * dt * 1.e6
        x_co2atm_13    = x_co2atm_13(1)

        if (ciso_14) then
          total_co2flux_14 = 0.  ! globally integrated air-sea 14CO2 flux (mol / s)
          call integrate_nod_2d_recom(0.001 * GloCO2flux_seaicemask_14, total_co2flux_14,     &
                                      MPI_COMM_FESOM, myDim_nod2D, eDim_nod2D, ulevels_nod2D, &
                                      areasvol)
!         Atmospheric radiocarbon budget in mol:
!         mol_co2atm_14 = mol_co2atm_14  + dt * (cosmic_14(1) - mol_co2atm_14 * lambda_14 - total_co2flux_14)
!                       = (mol_co2atm_14  + dt * (cosmic_14(1) - total_co2flux_14)) / (1 + lambda_14 * dt)
!         Budget in terms of the 14C / 12C volume mixing ratio
          x_co2atm_14(1) = (x_co2atm_14(1) + dt * (cosmic_14(1) - total_co2flux_14) / mol_allatm * 1.e6) / &
                           (1 + lambda_14 * dt)
          x_co2atm_14    = x_co2atm_14(1)

!         Adjust cosmogenic 14C production (mol / s) in spinup runs,
          r_atm_14        = x_co2atm_14(1) / x_co2atm(1)
!         r_atm_spinup_14 is calculated once-only in subroutine recom_init
          if (atbox_spinup .and. abs(r_atm_14 - r_atm_spinup_14) > 0.001) then
            cosmic_14(1) = cosmic_14(1) * (r_atm_spinup_14 / r_atm_14)
!           cosmic_14(1) = cosmic_14(1) * (1 + 0.01 * (r_atm_14_spinup / r_atm_14))
          end if
          cosmic_14 = cosmic_14(1)
        endif
      end if

      return
    end subroutine recom_atbox
