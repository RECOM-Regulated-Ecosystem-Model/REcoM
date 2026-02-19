module recom_extra
    interface
        subroutine Depth_calculations(n, nn, wf, zf, thick, recipthick, myDim_nod2D, eDim_nod2D, nl, hnode, zbar_3d_n)
            use recom_declarations, only: wp

            ! Input parameters
            integer, intent(in)                        :: n           ! Current node
            integer, intent(in)                        :: nn	    ! Total number of vertical nodes
            integer, intent(in)                        :: myDim_nod2D, eDim_nod2D, nl
            real(kind=WP), intent(in), dimension(:,:)  :: hnode, zbar_3d_n

            ! Output arrays
            real(kind=8), dimension(nl,6), intent(out) :: wf          ! [m/day] Flux velocities at the border of the control volumes
            real(kind=8), dimension(nl),   intent(out) :: zf          ! [m] Depth of vertical fluxes
            real(kind=8), dimension(nl-1), intent(out) :: thick       ! [m] Distance between two nodes = layer thickness
            real(kind=8), dimension(nl-1), intent(out) :: recipthick  ! [1/m] Reciprocal thickness
        end subroutine Depth_calculations

        subroutine Cobeta(daynew, ndpyr, myDim_nod2D, eDim_nod2D, geo_coord_nod2D)
            use recom_declarations, only: wp
            integer, intent(in)                       :: daynew, ndpyr, myDim_nod2D, eDim_nod2D
            real(kind=WP), intent(in), dimension(:,:) :: geo_coord_nod2D
        end subroutine Cobeta

        subroutine krill_resp(n, daynew, myDim_nod2D, eDim_nod2D, geo_coord_nod2D)
            use recom_declarations, only: wp
            integer, intent(in)                       :: n, daynew
            integer, intent(in)                       :: myDim_nod2D, eDim_nod2D
            real(kind=WP), intent(in), dimension(:,:) :: geo_coord_nod2D
        end subroutine krill_resp

        subroutine integrate_nod_2D_recom(data, int2D, MPI_COMM_FESOM, myDim_nod2D, eDim_nod2D, ulevels_nod2D, areasvol)
            use recom_declarations, only: wp

            implicit none

            integer,       intent(in)    :: MPI_COMM_FESOM, myDim_nod2D, eDim_nod2D
            real(kind=WP), intent(inout) :: int2D

            integer,       intent(in), dimension(:)   :: ulevels_nod2D
            real(kind=WP), intent(in), dimension(:)   :: data
            real(kind=WP), intent(in), dimension(:,:) :: areasvol
        end subroutine
    end interface
end module recom_extra

!===============================================================================
! Subroutine for calculating flux-depth and thickness of control volumes
!===============================================================================
subroutine Depth_calculations(n, nn, wf, zf, thick, recipthick, myDim_nod2D, eDim_nod2D, nl, hnode, zbar_3d_n)
    use recom_config
    use recom_declarations, only: wp

    implicit none

    ! Input parameters
    integer, intent(in)                       :: n           ! Current node
    integer, intent(in)                       :: nn	    ! Total number of vertical nodes
    integer, intent(in)                       :: myDim_nod2D, eDim_nod2D, nl
    real(kind=WP), intent(in), dimension(:,:) :: hnode, zbar_3d_n

    ! Output arrays
    real(kind=8), dimension(nl,6), intent(out)  :: wf          ! [m/day] Flux velocities at the border of the control volumes
    real(kind=8), dimension(nl),   intent(out)  :: zf          ! [m] Depth of vertical fluxes
    real(kind=8), dimension(nl-1), intent(out)  :: thick       ! [m] Distance between two nodes = layer thickness
    real(kind=8), dimension(nl-1), intent(out)  :: recipthick  ! [1/m] Reciprocal thickness

    ! Local variables
    integer                                     :: k           ! Layer index

! ======================================================================================
!! zbar(nl) allocate the array for storing the standard depths (depth of layers)
!! zbar is negative
!! Z(nl-1)  mid-depths of layers

    ! zbar_n: depth of layers due to ale thinkness variations at every node n
!    allocate(zbar_n(nl))
!    allocate(zbar_3d_n(nl,myDim_nod2D+eDim_nod2D))

    ! Z_n: mid depth of layers due to ale thinkness variations at every node n
!    allocate(Z_n(nl-1))
!    allocate(Z_3d_n(nl-1,myDim_nod2D+eDim_nod2D))
! ============================================================================== modular

    !! Background sinking speed
    wF(2:Nn, ivphy)   = VPhy      ! Phytoplankton sinking velocity
    wF(2:Nn, ivdia)   = VDia      ! Diatoms sinking velocity
    wF(2:Nn, ivdet)   = VDet      ! Detritus sinking velocity
    wF(2:Nn, ivdetsc) = VDet_zoo2 ! Second detritus sinking velocity
    wF(2:Nn, ivcoc)   = VCocco    ! Coccolithophores sinking velocity
    wF(2:Nn, ivpha)   = VPhaeo    ! Phaeocystis sinking velocity

    !! Boundary conditions (surface and bottom)
    wF(1,:)          = 0.d0
    wF(Nn+1,:)       = 0.d0

!if (allow_var_sinking) then
!!    wF(2:Nn+1,ivdet) = Vdet_a * abs(zbar_n(2:Nn+1)) + VDet
!!YY: use Vcalc instead of Vdet_a, only needed for calculating calc_diss
!    wF(2:Nn+1,ivdet) = Vcalc * abs(zbar_3d_n(2:Nn+1, n)) + VDet
!end if

    !! Calculate layer thickness and reciprocal of it
    thick   = 0.0_WP
    recipthick = 0.0_WP
    zf = 0.0_WP

    do k=1, nn
        thick(k) = hnode(k,n)
        if (hnode(k,n) > 0.0_WP) then
!        if thick(k) > 0.0_WP) then        ! alternative
         recipthick(k) = 1.0/hnode(k,n)
        else
         recipthick(k) = 0.0_WP
        end if
    end do

    !! set layer depth (negative)
    do k=1, nn+1
        zf(k)=zbar_3d_n(k,n)
    end do

end subroutine Depth_calculations

!===============================================================================
! Subroutine for calculating cos(AngleOfIncidence)
!===============================================================================
subroutine Cobeta(daynew, ndpyr, myDim_nod2D, eDim_nod2D, geo_coord_nod2D)
    use REcoM_GloVar
    use recom_declarations, only: wp, pi

    implicit none

    ! Local variables
    real(kind=8)                                     :: yearfrac              ! Fraction of year [0 1]
    real(kind=8)                                     :: yDay                  ! Year fraction in radians [0 2*pi]
    real(kind=8)                                     :: declination   = 0.d0  ! Declination of the sun at present lat and time
    real(kind=8)                                     :: CosAngleNoon  = 0.d0  ! Cosine of Angle of Incidence at noon
    integer                                          :: n

    ! Constants
    real(kind=8), parameter                          :: nWater        = 1.33  ! Refractive indices of water

    integer, intent(in)                       :: daynew, ndpyr, myDim_nod2D, eDim_nod2D
    real(kind=WP), intent(in), dimension(:,:) :: geo_coord_nod2D

!! find day (****NOTE for year starting in winter*****)
!! Paltridge, G. W. and C. M. R. Platt, Radiative Processes
!! in Meteorology and Climatology, Developments in
!! Atmospheric Sciences, vol. 5, Elsevier Scientific
!! Publishing Company, Amsterdam, Oxford,
!! New York, 1976, ISBN 0-444-41444-4.

    !! Calculate solar declination using Paltridge & Platt (1976) formula
    yearfrac    = mod(real(daynew), real(ndpyr)) / real(ndpyr)
    yDay        = 2.0d0 * pi * yearfrac

    declination = 0.006918                   &
         	- 0.399912 * cos(    yDay) &
     	        + 0.070257 * sin(    yDay) &
          	- 0.006758 * cos(2 * yDay) &
        	+ 0.000907 * sin(2 * yDay) &
        	- 0.002697 * cos(3 * yDay) &
        	+ 0.001480 * sin(3 * yDay)

    !! Calculate cosine of angle of incidence for each node
    do n = 1, myDim_nod2D
        cosAngleNoon = sin(geo_coord_nod2D(2, n)) * sin(declination) &
                     + cos(geo_coord_nod2D(2, n)) * cos(declination)

        cosAI(n)     = sqrt(1.0d0 - (1.0d0 - cosAngleNoon**2) / nWater**2)
    end do
end subroutine Cobeta
!================================================================================
! Calculating second zooplankton respiration rates
!================================================================================
subroutine krill_resp(n, daynew, myDim_nod2D, eDim_nod2D, geo_coord_nod2D)
    use REcoM_declarations
    use REcoM_LocVar
    use REcoM_GloVar

    implicit none

    ! Input parameters
    integer, intent(in)                       :: n, daynew
    integer, intent(in)                       :: myDim_nod2D, eDim_nod2D
    real(kind=WP), intent(in), dimension(:,:) :: geo_coord_nod2D

   ! Values from FESOM

   if (geo_coord_nod2D(2,n)<0.0_WP) then  !SH
      if ((daynew .le. 105)) then
       res_zoo2_a = 0.d0
      else if((105 .le. daynew).and.(daynew .le. 150)) then
       res_zoo2_a = (-1./90.)*daynew +7./6.
      else if((150 .lt. daynew).and.(daynew .lt. 250)) then
       res_zoo2_a = -0.5
      else if((250 .le. daynew).and.(daynew .le. 295)) then
       res_zoo2_a = (1/90.)*daynew - 59./18.
      else if((daynew .gt. 295)) then
       res_zoo2_a = 0.d0
      end if
   else
      if ((daynew .le. 65)) then
       res_zoo2_a = -0.5
      else if((285 .le. daynew).and.(daynew .le. 330)) then
       res_zoo2_a = (-1./90.)*daynew +57./18.
      else if((330 .lt. daynew)) then
       res_zoo2_a = -0.5
      else if((65 .le. daynew).and.(daynew .le. 110)) then
       res_zoo2_a = (1/90.)*daynew - 22./18.
      else if((110 .lt. daynew).and.(daynew .lt. 285)) then
       res_zoo2_a = 0.d0
      end if
   endif
 end subroutine krill_resp


subroutine integrate_nod_2D_recom(data, int2D, MPI_COMM_FESOM, myDim_nod2D, eDim_nod2D, ulevels_nod2D, areasvol)
    use recom_declarations, only: wp
    use mpi

    implicit none

    integer,       intent(in)    :: MPI_COMM_FESOM, myDim_nod2D, eDim_nod2D
    real(kind=WP), intent(inout) :: int2D

    integer,       intent(in), dimension(:)   :: ulevels_nod2D
    real(kind=WP), intent(in), dimension(:)   :: data
    real(kind=WP), intent(in), dimension(:,:) :: areasvol

    integer       :: row, MPIERR
    real(kind=WP) :: lval

    lval = 0.0_WP
    do row = 1, myDim_nod2D
       lval = lval + data(row) * areasvol(ulevels_nod2D(row), row)
    end do

    int2D = 0.0_WP
    call MPI_Allreduce(lval, int2D, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
                       MPI_COMM_FESOM, MPIerr)
end subroutine integrate_nod_2D_recom
