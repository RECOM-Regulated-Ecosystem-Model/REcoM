module recom_extra

    use recom_declarations, only: wp

    implicit none
    private

    public :: integrate_nod_2D_recom
    public :: krill_resp
    public :: Cobeta
    public :: Depth_calculations
    public :: calculate_solar_declination
    public :: solar_incidence_cosine

contains

    !===============================================================================
    ! Depth_calculations
    !   Computes layer thicknesses, reciprocal thicknesses, flux-point depths, and
    !   background sinking-velocity profiles for a single water column at node n.
    !
    !   Arguments:
    !     n          - node index
    !     nzmin      - top active level (1 for open ocean, >1 for ice-shelf cavity)
    !     nn         - bottom active level (= nlevels_nod2D(n) - 1)
    !     wf         - sinking velocity at flux interfaces [m/d], dimension (nl, 6)
    !     zf         - depth at flux interfaces [m, negative], dimension (nl)
    !     thick      - layer thickness [m], dimension (nl-1)
    !     recipthick - 1/thick [1/m], dimension (nl-1)
    !     nl         - number of vertical levels (interfaces)
    !     hnode      - layer thickness at nodes [m], dimension (nl-1, node)
    !     zbar_3d_n  - interface depths at nodes [m, negative], dimension (nl, node)
    !
    !   Cavity convention (nzmin > 1):
    !     - wf, thick, recipthick, and zf are only filled for the active range
    !       [nzmin, nn]; values outside are left at zero.
    !     - Surface boundary wf(nzmin, :) = 0 replaces wf(1, :) = 0.
    !     - Bottom boundary wf(nn+1, :) = 0.
    !===============================================================================
    subroutine Depth_calculations(n, nzmin, nn, wf, zf, thick, recipthick, nl, hnode, zbar_3d_n)

        use recom_config, only: ivcoc, ivdia, ivdet, ivdetsc, ivpha, ivphy, &
               VCalc, VDet, VDet_zoo2, VCocco, VPhaeo, VDia, VPhy

        implicit none

        ! Input parameters
        integer, intent(in) :: n       ! Node index
        integer, intent(in) :: nzmin   ! Top active level (1 = open ocean, >1 = cavity)
        integer, intent(in) :: nn      ! Bottom active level (= nlevels_nod2D(n) - 1)
        integer, intent(in) :: nl      ! Number of vertical levels
        real(kind=wp), intent(in), dimension(:, :) :: hnode, zbar_3d_n

        ! Output arrays
        ! Sinking velocities at flux interfaces [m/d]; index 1 = ivphy..ivpha
        real(kind=wp), dimension(nl, 6), intent(out) :: wf

        ! Depths at flux interfaces [m, negative]
        real(kind=wp), dimension(nl), intent(out) :: zf

        ! Layer thicknesses [m] and their reciprocals [1/m]
        real(kind=wp), dimension(nl - 1), intent(out) :: thick
        real(kind=wp), dimension(nl - 1), intent(out) :: recipthick

        ! Local variables
        integer :: k

        !---------------------------------------------------------------------------
        ! Initialise all output arrays to zero so unused levels stay clean.
        ! This also sets the surface (nzmin) and bottom (nn+1) boundary
        ! interfaces of wf to zero: no flux enters from above the top active
        ! level or exits below the bottom.
        !---------------------------------------------------------------------------
        wf         = 0.0_WP
        zf         = 0.0_WP
        thick      = 0.0_WP
        recipthick = 0.0_WP

        !---------------------------------------------------------------------------
        ! Background sinking velocities at interior flux interfaces
        !   Range nzmin+1:nn covers interfaces between active layers.
        !---------------------------------------------------------------------------
        wf(nzmin+1:nn, ivphy)   = VPhy      ! Small phytoplankton
        wf(nzmin+1:nn, ivdia)   = VDia      ! Diatoms
        wf(nzmin+1:nn, ivdet)   = VDet      ! Detritus class 1
        wf(nzmin+1:nn, ivdetsc) = VDet_zoo2 ! Detritus class 2
        wf(nzmin+1:nn, ivcoc)   = VCocco    ! Coccolithophores
        wf(nzmin+1:nn, ivpha)   = VPhaeo    ! Phaeocystis

        !if (allow_var_sinking) then
        !!YY: use Vcalc instead of Vdet_a, only needed for calculating calc_diss
        !!    Cavity-aware: interior interfaces only; bottom face (nn+1) stays 0.
        !!    Requires VCalc from recom_config.
        !    wf(nzmin+1:nn, ivdet) = VCalc * abs(zbar_3d_n(nzmin+1:nn, n)) + VDet
        !end if

        !---------------------------------------------------------------------------
        ! Layer thickness and its reciprocal over the active column
        !---------------------------------------------------------------------------
        do k = nzmin, nn
            thick(k) = hnode(k, n)
            if (hnode(k, n) > 0.0_WP) then
                recipthick(k) = 1.0_WP / hnode(k, n)
            else
                recipthick(k) = 0.0_WP
            end if
        end do

        !---------------------------------------------------------------------------
        ! Flux-interface depths over the active column (including bottom face)
        !---------------------------------------------------------------------------
        do k = nzmin, nn + 1
            zf(k) = zbar_3d_n(k, n)
        end do

    end subroutine Depth_calculations

    !===============================================================================
    ! calculate_solar_declination
    !   Solar declination [rad] after Paltridge & Platt (1976).
    !
    !   Note: P&P define the day angle with day 0 = 1 January. Here daynew = 1
    !   on 1 January, i.e. the angle is shifted by one day (<= ~0.4 deg in
    !   declination). Kept for consistency with legacy REcoM results.
    !===============================================================================
    subroutine calculate_solar_declination(daynew, ndpyr, declination)

        use recom_declarations, only: pi

        implicit none

        integer, intent(in) :: daynew, ndpyr
        real(kind=wp), intent(out) :: declination

        real(kind=wp) :: yearfrac ! Fraction of year [0 1]
        real(kind=wp) :: yDay     ! Year fraction in radians [0 2*pi]

        yearfrac = mod(real(daynew, WP), real(ndpyr, WP)) / real(ndpyr, WP)
        yDay = 2.0_WP * pi * yearfrac

        declination = 0.006918_WP &
                - 0.399912_WP * cos(yDay) &
                + 0.070257_WP * sin(yDay) &
                - 0.006758_WP * cos(2.0_WP * yDay) &
                + 0.000907_WP * sin(2.0_WP * yDay) &
                - 0.002697_WP * cos(3.0_WP * yDay) &
                + 0.001480_WP * sin(3.0_WP * yDay)
    end subroutine calculate_solar_declination

    !===============================================================================
    ! solar_incidence_cosine
    !   Cosine of the noon angle of incidence below the sea surface, after
    !   refraction (Snell's law, nWater = 1.33).
    !
    !   latitude and declination in radians.
    !   sin(lat)*sin(dec) + cos(lat)*cos(dec) = cos(lat - dec) = cos(noon zenith).
    !
    !   Note: in polar night cosAngleNoon < 0 (sun below horizon at noon). The
    !   squared form still returns a value >= sqrt(1 - 1/nWater**2) ~ 0.66; this
    !   is harmless as long as incoming PAR is zero there.
    !===============================================================================
    function solar_incidence_cosine(latitude, declination) result(cos_refr)

        implicit none

        real(kind=wp), intent(in) :: latitude, declination
        real(kind=wp) :: cos_refr

        ! Local variables
        real(kind=wp) :: cosAngleNoon

        ! Constants
        real(kind=wp), parameter :: nWater = 1.33_WP ! Refractive index of water

        cosAngleNoon = sin(latitude) * sin(declination) &
                + cos(latitude) * cos(declination)
        cos_refr = sqrt(1.0_WP - (1.0_WP - cosAngleNoon ** 2) / nWater ** 2)
    end function solar_incidence_cosine

    !===============================================================================
    ! Cobeta
    !   Computes cosAI(n): cosine of the angle of incidence of sunlight at the
    !   sea surface after refraction at the air-water interface, for each node.
    !
    !   Uses the solar declination formula of Paltridge & Platt (1976) and
    !   Snell's law with refractive index nWater = 1.33.
    !
    !   This routine is geometry-only and is cavity-safe - it does not access
    !   vertical levels.
    !
    !   Only owned nodes (1:myDim_nod2D) are filled; halo nodes of cosAI are
    !   not set here.
    !
    !   Reference:
    !     Paltridge & Platt (1976), Radiative Processes in Meteorology and
    !     Climatology, Developments in Atmospheric Sciences, vol. 5,
    !     Elsevier, ISBN 0-444-41444-4.
    !===============================================================================
    subroutine Cobeta(daynew, ndpyr, myDim_nod2D, geo_coord_nod2D)

        use REcoM_GloVar, only: cosAI

        implicit none

        integer, intent(in) :: daynew, ndpyr, myDim_nod2D
        real(kind=WP), intent(in), dimension(:, :) :: geo_coord_nod2D

        ! Local variables
        real(kind=wp) :: declination ! Declination of the sun at present time [rad]
        integer :: n

        call calculate_solar_declination(daynew, ndpyr, declination)

        do n = 1, myDim_nod2D
            cosAI(n) = solar_incidence_cosine(geo_coord_nod2D(2, n), declination)
        end do
    end subroutine Cobeta

    !================================================================================
    ! krill_resp
    !   Seasonal respiration modifier for the second zooplankton group (krill).
    !   Piecewise-linear in day of year, -0.5 in local winter, 0 in local summer,
    !   with 45-day linear transitions. Assumes a 365/366-day calendar.
    !
    !   Writes the module variable res_zoo2_a (REcoM_LocVar). Not thread-safe:
    !   call per node, immediately before use, outside OpenMP-parallel loops.
    !================================================================================
    subroutine krill_resp(daynew, node_latitude)
        use REcoM_LocVar, only: res_zoo2_a

        implicit none

        ! Input parameters
        integer,       intent(in)  :: daynew
        real(kind=WP), intent(in)  :: node_latitude

        real(kind=wp), parameter :: slope = 1.0_wp / 90.0_wp

        if (node_latitude < 0.0_WP) then ! Southern Hemisphere
            if (daynew <= 105) then
                res_zoo2_a = 0.0_WP
            else if (daynew <= 150) then                        ! 105 -> 150: 0 -> -0.5
                res_zoo2_a = -slope * daynew + 7.0_WP / 6.0_WP
            else if (daynew < 250) then                         ! austral winter
                res_zoo2_a = -0.5_WP
            else if (daynew <= 295) then                        ! 250 -> 295: -0.5 -> 0
                res_zoo2_a = slope * daynew - 59.0_WP / 18.0_WP
            else
                res_zoo2_a = 0.0_WP
            end if
        else                             ! Northern Hemisphere (incl. equator)
            if (daynew <= 65) then                              ! boreal winter
                res_zoo2_a = -0.5_WP
            else if (daynew <= 110) then                        ! 65 -> 110: -0.5 -> 0
                res_zoo2_a = slope * daynew - 22.0_WP / 18.0_WP
            else if (daynew < 285) then
                res_zoo2_a = 0.0_WP
            else if (daynew <= 330) then                        ! 285 -> 330: 0 -> -0.5
                res_zoo2_a = -slope * daynew + 57.0_WP / 18.0_WP
            else                                                ! boreal winter
                res_zoo2_a = -0.5_WP
            end if
        end if
    end subroutine krill_resp

    !================================================================================
    ! integrate_nod_2D_recom
    !   Global area integral of a 2D nodal field over owned nodes, using the
    !   area of the top active level (cavity-aware).
    !
    !   Note: MPI_DOUBLE_PRECISION assumes wp is real64.
    !================================================================================
    subroutine integrate_nod_2D_recom(data, int2D, MPI_COMM_FESOM, myDim_nod2D, ulevels_nod2D, &
            areasvol)

        use mpi

        implicit none

        integer, intent(in) :: MPI_COMM_FESOM, myDim_nod2D
        real(kind=WP), intent(out) :: int2D

        integer, intent(in), dimension(:) :: ulevels_nod2D
        real(kind=WP), intent(in), dimension(:) :: data
        real(kind=WP), intent(in), dimension(:, :) :: areasvol

        integer :: row, MPIerr
        real(kind=WP) :: lval

        lval = 0.0_WP
        do row = 1, myDim_nod2D
            lval = lval + data(row) * areasvol(ulevels_nod2D(row), row)
        end do

        int2D = 0.0_WP
        call MPI_Allreduce(lval, int2D, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
                MPI_COMM_FESOM, MPIerr)
    end subroutine integrate_nod_2D_recom

end module recom_extra
