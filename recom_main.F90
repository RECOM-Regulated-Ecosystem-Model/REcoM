! 24.03.2023
! OG
!===============================================================================
! Main REcoM
module recom_interface
    interface
            subroutine recom(ice_data_values, nl, ulevels_nod2D, nlevels_nod2D, hnode,           &
                             z_3d_n, zbar_3d_n, geo_coord_nod2D, ocean_area, areasvol, myDim_nod2d,      &
                             eDim_nod2D, mype, MPI_COMM_FESOM, tracers_info, num_tracers, tra_recom_sms, &
                             npes, sn, rn, s_mpitype_nod2D, r_mpitype_nod2D, s_mpitype_nod3D, r_mpitype_nod3D, &
                             sPE, rPE, requests, nreq, dt, daynew, month, mstep, ndpyr, yearold, timenew, rad, kappa, &
                             press_air, u_wind, v_wind, shortwave)
            use recom_glovar
            use recom_declarations, only: wp

            integer, intent(in)                              :: nl, myDim_nod2d, eDim_nod2D
            integer, intent(in)                              :: mype, MPI_COMM_FESOM, num_tracers
            integer, intent(in)                              :: mstep, daynew, yearold, month, ndpyr
            integer, intent(in), dimension(:)                :: ulevels_nod2D, nlevels_nod2D
            real(kind=WP), intent(inout)                     :: dt, kappa, timenew
            real(kind=WP), intent(in)                        :: ocean_area, rad
            real(kind=WP), intent(in), dimension(:)          :: ice_data_values, press_air
            real(kind=WP), intent(in), dimension(:)          :: u_wind, v_wind, shortwave
            real(kind=WP), intent(in), dimension(:, :)       :: hnode, z_3d_n, zbar_3d_n
            real(kind=WP), intent(in), dimension(:, :)       :: geo_coord_nod2D, areasvol
            real(kind=WP), intent(inout), dimension(:, :, :) :: tra_recom_sms

            integer, intent(in)                                 :: sn, rn, npes
            integer, intent(inout)                              :: nreq
            integer, intent(in),    dimension(:)                :: sPE, rPE
            integer, intent(inout), dimension(:)                :: requests
            integer, intent(in),    dimension(:),       pointer :: s_mpitype_nod2D, r_mpitype_nod2D
            integer, intent(in),    dimension(:, :, :), pointer :: s_mpitype_nod3D, r_mpitype_nod3D

            type(tracers_info_type), intent(in) :: tracers_info
        end subroutine
    end interface
end module

module bio_fluxes_interface
    interface
        subroutine bio_fluxes(alkalinity, MPI_COMM_FESOM, myDim_nod2D, eDim_nod2D, ocean_area, ulevels_nod2D, areasvol)
            use recom_declarations, only: wp

            integer, intent(in)               :: MPI_COMM_FESOM, myDim_nod2D, eDim_nod2D
            integer, intent(in), dimension(:) :: ulevels_nod2D

            real(kind=WP), intent(in)                 :: ocean_area
            real(kind=WP), intent(in), dimension(:,:) :: areasvol
            real(kind=WP), intent(in), dimension(:,:) :: alkalinity
      end subroutine
    end interface
end module

subroutine recom(ice_data_values, nl, ulevels_nod2D, nlevels_nod2D, hnode,           &
                 z_3d_n, zbar_3d_n, geo_coord_nod2D, ocean_area, areasvol, myDim_nod2d,      &
                 eDim_nod2D, mype, MPI_COMM_FESOM, tracers_info, num_tracers, tra_recom_sms, &
                 npes, sn, rn, s_mpitype_nod2D, r_mpitype_nod2D, s_mpitype_nod3D, r_mpitype_nod3D, &
                 sPE, rPE, requests, nreq, dt, daynew, month, mstep, ndpyr, yearold, timenew, rad, kappa, &
                 press_air, u_wind, v_wind, shortwave)

    use recom_g_comm_auto, only: recom_exchange_nod

    use recom_declarations
    use bio_fluxes_interface
    use recom_locvar
    use recom_glovar
    use recom_config
    use recom_ciso
    use recom_diags_management
    use recom_forcing_module
    use recom_atbox_module

    implicit none

    integer, intent(in)                              :: nl, myDim_nod2d, eDim_nod2D
    integer, intent(in)                              :: mype, MPI_COMM_FESOM, num_tracers
    integer, intent(in)                              :: mstep, daynew, yearold, month, ndpyr
    integer, intent(in), dimension(:)                :: ulevels_nod2D, nlevels_nod2D
    real(kind=WP), intent(inout)                     :: dt, kappa, timenew
    real(kind=WP), intent(in)                        :: ocean_area, rad
    real(kind=WP), intent(in), dimension(:)          :: ice_data_values, press_air
    real(kind=WP), intent(in), dimension(:)          :: u_wind, v_wind, shortwave
    real(kind=WP), intent(in), dimension(:, :)       :: hnode, z_3d_n, zbar_3d_n
    real(kind=WP), intent(in), dimension(:, :)       :: geo_coord_nod2D, areasvol
    real(kind=WP), intent(inout), dimension(:, :, :) :: tra_recom_sms

    ! These should all go into a dedicated REcoM type
    integer, intent(in)                                 :: sn, rn, npes
    integer, intent(inout)                              :: nreq
    integer, intent(in),    dimension(:)                :: sPE, rPE
    integer, intent(inout), dimension(:)                :: requests
    integer, intent(in),    dimension(:),       pointer :: s_mpitype_nod2D, r_mpitype_nod2D
    integer, intent(in),    dimension(:, :, :), pointer :: s_mpitype_nod3D, r_mpitype_nod3D

    type(tracers_info_type), intent(in) :: tracers_info

    !___________________________________________________________________________

! ======================================================================================
!! Depth information

!! zbar(nl) allocate the array for storing the standard depths, it is negativ
!! Z(nl-1)  mid-depths of cells

!! max. number of levels at node n
!! nzmax = nlevels_nod2D(n)
!! u_ice and v_ice are at nodes
!! u_w, v_w are at nodes (interpolated from elements)
!! u_wind and v_wind are always at nodes
! ======================================================================================

    real(kind=8)               :: SW, Loc_slp
    integer                    :: tr_num
    integer                    :: nz, n, nzmin, nzmax
    integer                    :: idiags

    real(kind=8)               :: Sali
    logical                    :: do_update = .false.

    real(kind=8),  allocatable :: Temp(:), Sali_depth(:), zr(:), PAR(:)
    real(kind=8),  allocatable :: C(:,:)

    !! * Mocsy *
    real(kind=8),  allocatable :: CO2_watercolumn(:)
    real(kind=8),  allocatable :: pH_watercolumn(:)
    real(kind=8),  allocatable :: pCO2_watercolumn(:)
    real(kind=8),  allocatable :: HCO3_watercolumn(:)

    !! * Diss *
    real(kind=8),  allocatable :: CO3_watercolumn(:)
    real(kind=8),  allocatable :: OmegaC_watercolumn(:)
    real(kind=8),  allocatable :: kspc_watercolumn(:)
    real(kind=8),  allocatable :: rhoSW_watercolumn(:)
    real(kind=WP)              :: ttf_rhs_bak (nl-1, num_tracers) ! local variable

    allocate(Temp(nl-1), Sali_depth(nl-1), zr(nl-1) , PAR(nl-1))
    allocate(C(nl-1, bgc_num))
    allocate(CO2_watercolumn(nl-1), pH_watercolumn(nl-1), pCO2_watercolumn(nl-1) , HCO3_watercolumn(nl-1))
    allocate(CO3_watercolumn(nl-1), OmegaC_watercolumn(nl-1), kspc_watercolumn(nl-1) , rhoSW_watercolumn(nl-1))

    !< ice concentration [0 to 1]

    !< alkalinity restoring to climatology
    !< virtual flux is possible

    if (restore_alkalinity) call bio_fluxes(tracers_info%data_pointers(2+ialk)%tracer_data(:,:), MPI_COMM_FESOM, &
                                            myDim_nod2D, eDim_nod2D, ocean_area,  &
                                            ulevels_nod2D, areasvol)
    if (recom_debug .and. mype==0) print *, achar(27)//'[36m'//'     --> bio_fluxes'//achar(27)//'[0m'

  if (use_atbox) then    ! MERGE
! Prognostic atmospheric isoCO2
    call recom_atbox(MPI_COMM_FESOM, myDim_nod2D, &
                     eDim_nod2D, ulevels_nod2D,  &
                     areasvol, dt)
!   optional I/O of isoCO2 and inferred cosmogenic 14C production; this may cost some CPU time
    if (ciso .and. ciso_14) then
      if ((daynew == ndpyr) .and. (timenew==86400.)) then
         do_update=.true.
      else
         do_update=.false.
      endif

      if (do_update .and. mype==0) write (*, fmt = '(a50,2x,i6,4(2x,f6.2))') &
                                         'Year, xCO2 (ppm), cosmic 14C flux (at / cm² / s):', &
                                          yearold, x_co2atm(1), x_co2atm_13(1), x_co2atm_14(1), cosmic_14(1) * production_rate_to_flux_14
    end if
  end if
! ======================================================================================
!********************************* LOOP STARTS *****************************************

    do n=1, myDim_nod2D  ! needs exchange_nod in the end
!     if (ulevels_nod2D(n)>1) cycle
!       nzmin = ulevels_nod2D(n)

        !!---- Number of vertical layers
        nzmax = nlevels_nod2D(n)-1

        !!---- This is needed for piston velocity
        Loc_ice_conc = ice_data_values(n)

        !!---- Mean sea level pressure
#if defined (__oasis)
!!      MB: This is an ad-hoc patch for AWIESM-2.1 and needs to be improved:
!!      We should consider air pressure provided by ECHAM.
        Loc_slp = pa2atm
#else
        Loc_slp = press_air(n)
#endif

        !!---- Benthic layers
        LocBenthos(1:benthos_num) = Benthos(n,1:benthos_num)

        !!---- Local conc of [H+]-ions from last time step. Decleared and saved in LocVar.
        !!---- used as first guess for H+ conc. in subroutine CO2flux (provided by recom_init)
        Hplus = GloHplus(n)

        !!---- Interpolated wind from atmospheric forcing
#if defined (__oasis)
        !! Derive 10m-wind speed from wind stress fields, see module recom_ciso.
        !! This is an ad-hoc solution as long as 10m-winds are not handled from OASIS.
        Uloc = wind_10(stress_atmoce_x(n), stress_atmoce_y(n))
#else
        ULoc = sqrt(u_wind(n)**2+v_wind(n)**2)
#endif

        !!---- Atmospheric CO2 in LocVar
        LocAtmCO2 = AtmCO2(month)

! Update of prognostic atmospheric CO2 values
     if (use_atbox) then
       LocAtmCO2                   = x_co2atm(1)
       if (ciso) then
         LocAtmCO2_13              = x_co2atm_13(1)
         if (ciso_14) LocAtmCO2_14 = x_co2atm_14(1)
       end if
     else
! Consider prescribed atmospheric CO2 values
       if (ciso) then
         LocAtmCO2_13              = AtmCO2_13(month)
         if (ciso_14) then
!          Latitude of nodal point n
           lat_val = geo_coord_nod2D(2,n) / rad
!          Zonally binned NH / SH / TZ 14CO2 input values
           LocAtmCO2_14 = AtmCO2_14(lat_zone(lat_val), month)
         end if
       end if
     end if  ! use_atbox

     if (ciso) then
       r_atm_13                    = LocAtmCO2_13(1) / LocAtmCO2(1)
       if (ciso_14) r_atm_14       = LocAtmCO2_14(1) / LocAtmCO2(1)
     end if

        !!---- Shortwave penetration
        SW = parFrac * shortwave(n)
        SW = SW * (1.d0 - ice_data_values(n))

        !!---- Temperature in water column
        Temp(1:nzmax) = tracers_info%data_pointers(1)%tracer_data(1:nzmax, n)

        !!---- Surface salinity
        Sali                = tracers_info%data_pointers(2)%tracer_data(1, n)
        Sali_depth(1:nzmax) = tracers_info%data_pointers(2)%tracer_data(1:nzmax, n)


        !!---- CO2 in the watercolumn

        !! * Mocsy *
        CO2_watercolumn(1:nzmax)    = CO23D(1:nzmax, n)
        pH_watercolumn(1:nzmax)     = pH3D(1:nzmax, n)
        pCO2_watercolumn(1:nzmax)   = pCO23D(1:nzmax, n)
        HCO3_watercolumn(1:nzmax)   = HCO33D(1:nzmax, n)
        !! * Diss *
        CO3_watercolumn(1:nzmax)    = CO33D(1:nzmax, n)
        OmegaC_watercolumn(1:nzmax) = OmegaC3D(1:nzmax, n)
        kspc_watercolumn(1:nzmax)   = kspc3D(1:nzmax, n)
        rhoSW_watercolumn(1:nzmax)  = rhoSW3D(1:nzmax, n)

        !!---- Biogeochemical tracers
        do tr_num = num_tracers-bgc_num+1, num_tracers
            C(1:nzmax, tr_num-2) = tracers_info%data_pointers(tr_num)%tracer_data(1:nzmax, n)
        end do

        ttf_rhs_bak = 0.0


        do tr_num=1, num_tracers
            if (tracers_info%ltra_diag(tr_num)) then
                ttf_rhs_bak(1:nzmax,tr_num) = tracers_info%data_pointers(tr_num)%tracer_data(1:nzmax, n)
            end if
        end do

        !!---- Depth of the nodes in the water column
        zr(1:nzmax) = Z_3d_n(1:nzmax, n)

        !!---- The PAR in the local water column is initialized
        PAR(1:nzmax) = 0.d0

        !!---- ice_data_values(row): Ice concentration in the local node
        FeDust = GloFeDust(n) * (1.d0 - ice_data_values(n)) * dust_sol
        NDust  = GloNDust(n)  * (1.d0 - ice_data_values(n))

        if (Diags) then
            ! Allocate and initialize all diagnostic arrays for a water column
            call allocate_and_init_diags(nl)
        end if
!
!        !! * Allocate 3D diagnostics *
!            allocate(vertrespmeso(nl-1))
!            vertrespmeso  = 0.d0
!
!if (enable_3zoo2det) then
!            allocate(vertrespmacro(nl-1), vertrespmicro(nl-1))
!            vertrespmacro = 0.d0
!            vertrespmicro = 0.d0
!endif
!            allocate(vertcalcdiss(nl-1), vertcalcif(nl-1))
!            vertcalcdiss = 0.d0
!            vertcalcif   = 0.d0
!
!            allocate(vertaggn(nl-1), vertaggd(nl-1))
!            vertaggn = 0.d0
!            vertaggd = 0.d0
!
!            allocate(vertdocexn(nl-1), vertdocexd(nl-1))
!            vertdocexn = 0.d0
!            vertdocexd = 0.d0
!
!            allocate(vertrespn(nl-1), vertrespd(nl-1))
!            vertrespn = 0.d0
!            vertrespd = 0.d0
!
!            allocate(VTPhyCO2(nl-1), VTDiaCO2(nl-1))
!            VTPhyCO2 = 0.d0
!            VTDiaCO2 = 0.d0
!
!            allocate(VTCphotLigLim_phyto(nl-1), VTCphotLigLim_diatoms(nl-1))
!            VTCphotLigLim_phyto = 0.d0
!            VTCphotLigLim_diatoms = 0.d0
!
!            allocate(VTCphot_phyto(nl-1), VTCphot_diatoms(nl-1))
!            VTCphot_phyto = 0.d0
!            VTCphot_diatoms = 0.d0
!
!if (enable_coccos) then
!            allocate(VTTemp_diatoms(nl-1), VTTemp_phyto(nl-1))
!            VTTemp_diatoms = 0.d0
!            VTTemp_phyto = 0.d0
!
!            allocate(VTqlimitFac_phyto(nl-1), VTqlimitFac_diatoms(nl-1))
!            VTqlimitFac_phyto = 0.d0
!            VTqlimitFac_diatoms = 0.d0
!
!            allocate(VTSi_assimDia(nl-1))
!            VTSi_assimDia = 0.d0
!
!            allocate(vertaggc(nl-1), vertdocexc(nl-1), vertrespc(nl-1))
!            vertaggc = 0.d0
!            vertdocexc = 0.d0
!            vertrespc = 0.d0
!
!            allocate(vertaggp(nl-1), vertdocexp(nl-1), vertrespp(nl-1)) ! Phaeocystis
!            vertaggp = 0.d0
!            vertdocexp = 0.d0
!            vertrespp = 0.d0
!
!            allocate(VTTemp_cocco(nl-1), VTTemp_phaeo(nl-1))
!            VTTemp_cocco = 0.d0
!            VTTemp_phaeo = 0.d0
!
!            allocate(VTCoccoCO2(nl-1), VTPhaeoCO2(nl-1))
!            VTCoccoCO2 = 0.d0
!            VTPhaeoCO2 = 0.d0
!
!            allocate(VTqlimitFac_cocco(nl-1), VTqlimitFac_phaeo(nl-1))
!            VTqlimitFac_cocco  = 0.d0
!            VTqlimitFac_phaeo  = 0.d0
!
!            allocate(VTCphotLigLim_cocco(nl-1), VTCphotLigLim_phaeo(nl-1))
!            VTCphotLigLim_cocco  = 0.d0
!            VTCphotLigLim_phaeo  = 0.d0
!
!            allocate(VTCphot_cocco(nl-1), VTCphot_phaeo(nl-1))
!            VTCphot_cocco  = 0.d0
!            VTCphot_phaeo  = 0.d0
!
!
!endif
!
!            !! * Allocate 2D diagnostics *
!            allocate(vertNPPn(nl-1), vertGPPn(nl-1), vertNNAn(nl-1), vertChldegn(nl-1))
!            vertNPPn = 0.d0
!            vertGPPn = 0.d0
!            vertNNAn = 0.d0
!            vertChldegn  = 0.d0
!
!            allocate(vertNPPd(nl-1), vertGPPd(nl-1), vertNNAd(nl-1), vertChldegd(nl-1))
!            vertNPPd = 0.d0
!            vertGPPd = 0.d0
!            vertNNAd = 0.d0
!            vertChldegd  = 0.d0
!
!if (enable_coccos) then
!            allocate(vertNPPc(nl-1), vertGPPc(nl-1), vertNNAc(nl-1), vertChldegc(nl-1))
!            vertNPPc = 0.d0
!            vertGPPc = 0.d0
!            vertNNAc = 0.d0
!            vertChldegc = 0.d0
!
!            allocate(vertNPPp(nl-1), vertGPPp(nl-1), vertNNAp(nl-1), vertChldegp(nl-1))
!            vertNPPp = 0.d0
!            vertGPPp = 0.d0
!            vertNNAp = 0.d0
!            vertChldegp = 0.d0
!endif

        if (recom_debug .and. mype==0) print *, achar(27)//'[36m'//'     --> REcoM_Forcing'//achar(27)//'[0m'

! ======================================================================================
!******************************** RECOM FORCING ****************************************
        call REcoM_Forcing(n, nzmax, C, SW, Loc_slp, Temp, Sali, Sali_depth, &
                           CO2_watercolumn,         & ! NEW MOCSY CO2 for the whole watercolumn
                           pH_watercolumn,          & ! NEW MOCSY pH for the whole watercolumn
                           pCO2_watercolumn,        & ! NEW MOCSY pCO2 for the whole watercolumn
                           HCO3_watercolumn,        & ! NEW MOCSY HCO3 for the whole watercolumn
                           CO3_watercolumn,         & ! NEW DISS CO3 for the whole watercolumn
                           OmegaC_watercolumn,      & ! NEW DISS OmegaC for the whole watercolumn
                           kspc_watercolumn,        & ! NEW DISS stoichiometric solubility product for calcite [mol^2/kg^2]
                           rhoSW_watercolumn,       & ! NEW DISS in-situ density of seawater [mol/m^3]
                           PAR, MPI_COMM_FESOM, mype, myDim_nod2D, &
                           eDim_nod2D, nl, hnode, zbar_3d_n,    &
                           geo_coord_nod2D, daynew, ndpyr, dt, kappa, mstep, rad)

        do tr_num = num_tracers-bgc_num+1, num_tracers !bgc_num+2
            tracers_info%data_pointers(tr_num)%tracer_data(1:nzmax, n) = C(1:nzmax, tr_num-2)
        end do

        ! recom_sms

           do tr_num=1, num_tracers
               if (tracers_info%ltra_diag(tr_num)) then
                   tra_recom_sms(1:nzmax,n,tr_num) = tracers_info%data_pointers(tr_num)%tracer_data(1:nzmax, n) - ttf_rhs_bak(1:nzmax,tr_num)
             !if (mype==0)  print *,  tra_recom_sms(:,:,tr_num)
               end if

           end do

        !!---- Local variables that have been changed during the time-step are stored so they can be saved
        Benthos(n,1:benthos_num) = LocBenthos(1:benthos_num)
        GlodecayBenthos(n, 1:benthos_num) = decayBenthos(1:benthos_num)/SecondsPerDay ! convert from [mmol/m2/d] to [mmol/m2/s]

        if (recom_debug .and. mype==0) print *, achar(27)//'[36m'//'     --> ciso after REcoM_Forcing'//achar(27)//'[0m'

        if (Diags) then
           call update_2d_diags(n)
           call update_3d_diags(n, nzmax)
           ! Deallocate vertical tracer array
           call deallocate_diags()
        endif
!            !! * Update 2D diagnostics *
!            NPPn(n) = locNPPn
!            NPPd(n) = locNPPd
!            GPPn(n) = locGPPn
!            GPPd(n) = locGPPd
!            NNAn(n) = locNNAn
!            NNAd(n) = locNNAd
!            Chldegn(n) = locChldegn
!            Chldegd(n) = locChldegd
!
!if (enable_coccos) then
!            NPPc(n) = locNPPc
!            GPPc(n) = locGPPc
!            NNAc(n) = locNNAc
!            Chldegc(n) = locChldegc
!            NPPp(n) = locNPPp
!            GPPp(n) = locGPPp
!            NNAp(n) = locNNAp
!            Chldegp(n) = locChldegp
!endif
!
!            !! * Update 3D diagnostics *
!            respmeso     (1:nzmax,n) = vertrespmeso     (1:nzmax)
!if (enable_3zoo2det) then
!            respmacro    (1:nzmax,n) = vertrespmacro    (1:nzmax)
!            respmicro    (1:nzmax,n) = vertrespmicro    (1:nzmax)
!endif
!            calcdiss     (1:nzmax,n) = vertcalcdiss     (1:nzmax)
!            calcif       (1:nzmax,n) = vertcalcif       (1:nzmax)
!
!            aggn         (1:nzmax,n) = vertaggn         (1:nzmax)
!            docexn       (1:nzmax,n) = vertdocexn       (1:nzmax)
!            respn        (1:nzmax,n) = vertrespn        (1:nzmax)
!            NPPn3D       (1:nzmax,n) = vertNPPn         (1:nzmax)
!
!            aggd         (1:nzmax,n) = vertaggd         (1:nzmax)
!            docexd       (1:nzmax,n) = vertdocexd       (1:nzmax)
!            respd        (1:nzmax,n) = vertrespd        (1:nzmax)
!            NPPd3D       (1:nzmax,n) = vertNPPd         (1:nzmax)
!
!            TPhyCO2             (1:nzmax,n) = VTPhyCO2                  (1:nzmax)
!            TDiaCO2             (1:nzmax,n) = VTDiaCO2                  (1:nzmax)
!            TCphotLigLim_phyto  (1:nzmax,n) = VTCphotLigLim_phyto       (1:nzmax)
!            TCphotLigLim_diatoms(1:nzmax,n) = VTCphotLigLim_diatoms     (1:nzmax)
!            TCphot_phyto        (1:nzmax,n) = VTCphot_phyto             (1:nzmax)
!            TCphot_diatoms      (1:nzmax,n) = VTCphot_diatoms           (1:nzmax)
!
!if (enable_coccos) then
!
!            TTemp_phyto         (1:nzmax,n) = VTTemp_phyto              (1:nzmax)
!            TqlimitFac_phyto    (1:nzmax,n) = VTqlimitFac_phyto         (1:nzmax)
!            TTemp_diatoms       (1:nzmax,n) = VTTemp_diatoms            (1:nzmax) !! NEW from here tracking vars
!            TqlimitFac_diatoms  (1:nzmax,n) = VTqlimitFac_diatoms       (1:nzmax)
!            TSi_assimDia        (1:nzmax,n) = VTSi_assimDia             (1:nzmax)
!
!            aggc         (1:nzmax,n) = vertaggc         (1:nzmax)
!            docexc       (1:nzmax,n) = vertdocexc       (1:nzmax)
!            respc        (1:nzmax,n) = vertrespc        (1:nzmax)
!            NPPc3D       (1:nzmax,n) = vertNPPc         (1:nzmax)
!
!            aggp         (1:nzmax,n) = vertaggp         (1:nzmax)
!            docexp       (1:nzmax,n) = vertdocexp       (1:nzmax)
!            respp        (1:nzmax,n) = vertrespp        (1:nzmax)
!            NPPp3D       (1:nzmax,n) = vertNPPp         (1:nzmax)
!
!            TTemp_cocco         (1:nzmax,n) = VTTemp_cocco              (1:nzmax)
!            TCoccoCO2           (1:nzmax,n) = VTCoccoCO2                (1:nzmax)
!            TqlimitFac_cocco    (1:nzmax,n) = VTqlimitFac_cocco         (1:nzmax)
!            TCphotLigLim_cocco  (1:nzmax,n) = VTCphotLigLim_cocco       (1:nzmax)
!            TCphot_cocco        (1:nzmax,n) = VTCphot_cocco             (1:nzmax)
!
!            TTemp_phaeo         (1:nzmax,n) = VTTemp_phaeo              (1:nzmax)
!            TPhaeoCO2           (1:nzmax,n) = VTPhaeoCO2                (1:nzmax)
!            TqlimitFac_phaeo    (1:nzmax,n) = VTqlimitFac_phaeo         (1:nzmax)
!            TCphotLigLim_phaeo  (1:nzmax,n) = VTCphotLigLim_phaeo       (1:nzmax)
!            TCphot_phaeo        (1:nzmax,n) = VTCphot_phaeo             (1:nzmax)
!
!    endif
!
!if (recom_debug .and. mype==0) print *, achar(27)//'[36m'//'     --> ciso after REcoM_Forcing'//achar(27)//'[0m'
!
!            !! * Deallocating 2D diagnostics *
!            deallocate(vertNPPn, vertGPPn, vertNNAn, vertChldegn)
!            deallocate(vertNPPd, vertGPPd, vertNNAd, vertChldegd)
!if (enable_coccos) then
!            deallocate(vertNPPc, vertGPPc, vertNNAc, vertChldegc)
!            deallocate(vertNPPp, vertGPPp, vertNNAp, vertChldegp)
!endif
!
!            !! * Deallocating 3D Diagnostics *
!            deallocate(vertrespmeso)
!if (enable_3zoo2det) then
!            deallocate(vertrespmacro, vertrespmicro)
!endif
!            deallocate(vertcalcdiss, vertcalcif)
!            deallocate(vertaggn, vertdocexn, vertrespn)
!            deallocate(vertaggd, vertdocexd, vertrespd)
!            deallocate(VTPhyCO2, VTCphotLigLim_phyto, VTCphot_phyto)
!            deallocate(VTDiaCO2, VTCphotLigLim_diatoms, VTCphot_diatoms)
!
!if (enable_coccos) then
!           deallocate(vertgrazmeso_c)
!            deallocate(VTTemp_phyto, VTqlimitFac_phyto)
!            deallocate(VTTemp_diatoms, VTqlimitFac_diatoms)
!
!            deallocate(VTSi_assimDia)
!            deallocate(vertaggc, vertdocexc, vertrespc)
!            deallocate(vertaggp, vertdocexp, vertrespp)
!
!            deallocate(VTTemp_cocco, VTCoccoCO2, VTqlimitFac_cocco, VTCphotLigLim_cocco, VTCphot_cocco)
!            deallocate(VTTemp_phaeo, VTPhaeoCO2, VTqlimitFac_phaeo, VTCphotLigLim_phaeo, VTCphot_phaeo)
!
!endif

        AtmFeInput(n)            = FeDust
        AtmNInput(n)             = NDust
        GloHplus(n)              = ph(1)
        PistonVelocity(n)        = kw660(1)
        alphaCO2(n)              = K0(1)

        GloPCO2surf(n)           = pco2surf(1)
        GlodPCO2surf(n)          = dpco2surf(1)
        GloCO2flux(n)            = dflux(1)                   !  [mmol/m2/d]
        GloCO2flux_seaicemask(n) = co2flux_seaicemask(1)      !  [mmol/m2/s]
        GloO2flux_seaicemask(n)  = o2flux_seaicemask(1)       !  [mmol/m2/s]
    if (ciso) then
        GloCO2flux_seaicemask_13(n)     = co2flux_seaicemask_13(1)        !  [mmol/m2/s]
        if (ciso_14) then
            GloCO2flux_seaicemask_14(n) = co2flux_seaicemask_14(1)        !  [mmol/m2/s]
        end if
     end if
        GloO2flux(n)             = oflux(1)                   !  [mmol/m2/d]

        PAR3D(1:nzmax,n)         = PAR(1:nzmax)

        !! * Mocsy *
        CO23D(1:nzmax,n)         = CO2_watercolumn(1:nzmax)
        pH3D(1:nzmax,n)          = pH_watercolumn(1:nzmax)
        pCO23D(1:nzmax,n)        = pCO2_watercolumn(1:nzmax)
        HCO33D(1:nzmax,n)        = HCO3_watercolumn(1:nzmax)

        !! * Diss *
        CO33D(1:nzmax,n)         = CO3_watercolumn(1:nzmax)
        OmegaC3D(1:nzmax,n)      = OmegaC_watercolumn(1:nzmax)
        kspc3D(1:nzmax,n)        = kspc_watercolumn(1:nzmax)
        rhoSW3D(1:nzmax,n)       = rhoSW_watercolumn(1:nzmax)
    end do

! ======================================================================================
!************************** EXCHANGE NODAL INFORMATION *********************************

    do tr_num=num_tracers-bgc_num+1, num_tracers
        call recom_exchange_nod(tracers_info%data_pointers(tr_num)%tracer_data(:,:), &
                                npes, sn, rn, MPI_COMM_FESOM, mype, s_mpitype_nod3D, &
                                r_mpitype_nod3D, sPE, rPE, requests, nreq)
    end do

    call recom_exchange_nod(GloPCO2surf, &
                            npes, sn, rn, MPI_COMM_FESOM, mype, s_mpitype_nod2D, &
                            r_mpitype_nod2D, sPE, rPE, requests, nreq)

    call recom_exchange_nod(GlodPCO2surf, &
                            npes, sn, rn, MPI_COMM_FESOM, mype, s_mpitype_nod2D, &
                            r_mpitype_nod2D, sPE, rPE, requests, nreq)

    call recom_exchange_nod(GloCO2flux, &
                            npes, sn, rn, MPI_COMM_FESOM, mype, s_mpitype_nod2D, &
                            r_mpitype_nod2D, sPE, rPE, requests, nreq)

    call recom_exchange_nod(GloCO2flux_seaicemask, &
                            npes, sn, rn, MPI_COMM_FESOM, mype, s_mpitype_nod2D, &
                            r_mpitype_nod2D, sPE, rPE, requests, nreq)

    call recom_exchange_nod(GloO2flux_seaicemask, &
                            npes, sn, rn, MPI_COMM_FESOM, mype, s_mpitype_nod2D, &
                            r_mpitype_nod2D, sPE, rPE, requests, nreq)
    if (ciso) then
        call recom_exchange_nod(GloPCO2surf_13, &
                                npes, sn, rn, MPI_COMM_FESOM, mype, s_mpitype_nod2D, &
                                r_mpitype_nod2D, sPE, rPE, requests, nreq)

        call recom_exchange_nod(GloCO2flux_13, &
                                npes, sn, rn, MPI_COMM_FESOM, mype, s_mpitype_nod2D, &
                                r_mpitype_nod2D, sPE, rPE, requests, nreq)

        call recom_exchange_nod(GloCO2flux_seaicemask_13, &
                                npes, sn, rn, MPI_COMM_FESOM, mype, s_mpitype_nod2D, &
                                r_mpitype_nod2D, sPE, rPE, requests, nreq)
        if (ciso_14) then
            call recom_exchange_nod(GloPCO2surf_14, &
                                    npes, sn, rn, MPI_COMM_FESOM, mype, s_mpitype_nod2D, &
                                    r_mpitype_nod2D, sPE, rPE, requests, nreq)

            call recom_exchange_nod(GloCO2flux_14, &
                                    npes, sn, rn, MPI_COMM_FESOM, mype, s_mpitype_nod2D, &
                                    r_mpitype_nod2D, sPE, rPE, requests, nreq)

            call recom_exchange_nod(GloCO2flux_seaicemask_14, &
                                    npes, sn, rn, MPI_COMM_FESOM, mype, s_mpitype_nod2D, &
                                    r_mpitype_nod2D, sPE, rPE, requests, nreq)
        end if
    end if
    do n=1, benthos_num
        call recom_exchange_nod(Benthos(:,n), &
                                npes, sn, rn, MPI_COMM_FESOM, mype, s_mpitype_nod2D, &
                                r_mpitype_nod2D, sPE, rPE, requests, nreq)
    end do

    if (Diags) then
        call recom_exchange_nod(NPPn, &
                                npes, sn, rn, MPI_COMM_FESOM, mype, s_mpitype_nod2D, &
                                r_mpitype_nod2D, sPE, rPE, requests, nreq)
        call recom_exchange_nod(NPPd, &
                                npes, sn, rn, MPI_COMM_FESOM, mype, s_mpitype_nod2D, &
                                r_mpitype_nod2D, sPE, rPE, requests, nreq)
        call recom_exchange_nod(GPPn, &
                                npes, sn, rn, MPI_COMM_FESOM, mype, s_mpitype_nod2D, &
                                r_mpitype_nod2D, sPE, rPE, requests, nreq)
        call recom_exchange_nod(GPPd, &
                                npes, sn, rn, MPI_COMM_FESOM, mype, s_mpitype_nod2D, &
                                r_mpitype_nod2D, sPE, rPE, requests, nreq)
        call recom_exchange_nod(NNAn, &
                                npes, sn, rn, MPI_COMM_FESOM, mype, s_mpitype_nod2D, &
                                r_mpitype_nod2D, sPE, rPE, requests, nreq)
        call recom_exchange_nod(NNAd, &
                                npes, sn, rn, MPI_COMM_FESOM, mype, s_mpitype_nod2D, &
                                r_mpitype_nod2D, sPE, rPE, requests, nreq)
        call recom_exchange_nod(Chldegn, &
                                npes, sn, rn, MPI_COMM_FESOM, mype, s_mpitype_nod2D, &
                                r_mpitype_nod2D, sPE, rPE, requests, nreq)
        call recom_exchange_nod(Chldegd, &
                                npes, sn, rn, MPI_COMM_FESOM, mype, s_mpitype_nod2D, &
                                r_mpitype_nod2D, sPE, rPE, requests, nreq)
        if (enable_coccos) then
            call recom_exchange_nod(NPPc, &
                                    npes, sn, rn, MPI_COMM_FESOM, mype, s_mpitype_nod2D, &
                                    r_mpitype_nod2D, sPE, rPE, requests, nreq)
            call recom_exchange_nod(GPPc, &
                                    npes, sn, rn, MPI_COMM_FESOM, mype, s_mpitype_nod2D, &
                                    r_mpitype_nod2D, sPE, rPE, requests, nreq)
            call recom_exchange_nod(NNAc, &
                                    npes, sn, rn, MPI_COMM_FESOM, mype, s_mpitype_nod2D, &
                                    r_mpitype_nod2D, sPE, rPE, requests, nreq)
            call recom_exchange_nod(Chldegc, &
                                    npes, sn, rn, MPI_COMM_FESOM, mype, s_mpitype_nod2D, &
                                    r_mpitype_nod2D, sPE, rPE, requests, nreq)
            call recom_exchange_nod(NPPp, &
                                    npes, sn, rn, MPI_COMM_FESOM, mype, s_mpitype_nod2D, &
                                    r_mpitype_nod2D, sPE, rPE, requests, nreq)
            call recom_exchange_nod(GPPp, &
                                    npes, sn, rn, MPI_COMM_FESOM, mype, s_mpitype_nod2D, &
                                    r_mpitype_nod2D, sPE, rPE, requests, nreq)
            call recom_exchange_nod(NNAp, &
                                    npes, sn, rn, MPI_COMM_FESOM, mype, s_mpitype_nod2D, &
                                    r_mpitype_nod2D, sPE, rPE, requests, nreq)
            call recom_exchange_nod(Chldegp, &
                                    npes, sn, rn, MPI_COMM_FESOM, mype, s_mpitype_nod2D, &
                                    r_mpitype_nod2D, sPE, rPE, requests, nreq)
        endif
        call recom_exchange_nod(grazmeso_tot, &
                                npes, sn, rn, MPI_COMM_FESOM, mype, s_mpitype_nod2D, &
                                r_mpitype_nod2D, sPE, rPE, requests, nreq)
        call recom_exchange_nod(grazmeso_n, &
                                npes, sn, rn, MPI_COMM_FESOM, mype, s_mpitype_nod2D, &
                                r_mpitype_nod2D, sPE, rPE, requests, nreq)
        call recom_exchange_nod(grazmeso_d, &
                                npes, sn, rn, MPI_COMM_FESOM, mype, s_mpitype_nod2D, &
                                r_mpitype_nod2D, sPE, rPE, requests, nreq)
        if (enable_coccos) then
            call recom_exchange_nod(grazmeso_c, &
                                    npes, sn, rn, MPI_COMM_FESOM, mype, s_mpitype_nod2D, &
                                    r_mpitype_nod2D, sPE, rPE, requests, nreq)
            call recom_exchange_nod(grazmeso_p, &
                                    npes, sn, rn, MPI_COMM_FESOM, mype, s_mpitype_nod2D, &
                                    r_mpitype_nod2D, sPE, rPE, requests, nreq)
        endif
        call recom_exchange_nod(grazmeso_det, &
                                npes, sn, rn, MPI_COMM_FESOM, mype, s_mpitype_nod2D, &
                                r_mpitype_nod2D, sPE, rPE, requests, nreq)
        if (enable_3zoo2det) then
            call recom_exchange_nod(grazmeso_mic, &
                                    npes, sn, rn, MPI_COMM_FESOM, mype, s_mpitype_nod2D, &
                                    r_mpitype_nod2D, sPE, rPE, requests, nreq)
            call recom_exchange_nod(grazmeso_det2, &
                                    npes, sn, rn, MPI_COMM_FESOM, mype, s_mpitype_nod2D, &
                                    r_mpitype_nod2D, sPE, rPE, requests, nreq)
            call recom_exchange_nod(grazmacro_tot, &
                                    npes, sn, rn, MPI_COMM_FESOM, mype, s_mpitype_nod2D, &
                                    r_mpitype_nod2D, sPE, rPE, requests, nreq)
            call recom_exchange_nod(grazmacro_n, &
                                    npes, sn, rn, MPI_COMM_FESOM, mype, s_mpitype_nod2D, &
                                    r_mpitype_nod2D, sPE, rPE, requests, nreq)
            call recom_exchange_nod(grazmacro_d, &
                                    npes, sn, rn, MPI_COMM_FESOM, mype, s_mpitype_nod2D, &
                                    r_mpitype_nod2D, sPE, rPE, requests, nreq)
            if (enable_coccos) then
                call recom_exchange_nod(grazmacro_c, &
                                        npes, sn, rn, MPI_COMM_FESOM, mype, s_mpitype_nod2D, &
                                        r_mpitype_nod2D, sPE, rPE, requests, nreq)
                call recom_exchange_nod(grazmacro_p, &
                                        npes, sn, rn, MPI_COMM_FESOM, mype, s_mpitype_nod2D, &
                                        r_mpitype_nod2D, sPE, rPE, requests, nreq)
            endif
            call recom_exchange_nod(grazmacro_mes, &
                                    npes, sn, rn, MPI_COMM_FESOM, mype, s_mpitype_nod2D, &
                                    r_mpitype_nod2D, sPE, rPE, requests, nreq)
            call recom_exchange_nod(grazmacro_det, &
                                    npes, sn, rn, MPI_COMM_FESOM, mype, s_mpitype_nod2D, &
                                    r_mpitype_nod2D, sPE, rPE, requests, nreq)
            call recom_exchange_nod(grazmacro_mic, &
                                    npes, sn, rn, MPI_COMM_FESOM, mype, s_mpitype_nod2D, &
                                    r_mpitype_nod2D, sPE, rPE, requests, nreq)
            call recom_exchange_nod(grazmacro_det2, &
                                    npes, sn, rn, MPI_COMM_FESOM, mype, s_mpitype_nod2D, &
                                    r_mpitype_nod2D, sPE, rPE, requests, nreq)
            call recom_exchange_nod(grazmicro_tot, &
                                    npes, sn, rn, MPI_COMM_FESOM, mype, s_mpitype_nod2D, &
                                    r_mpitype_nod2D, sPE, rPE, requests, nreq)
            call recom_exchange_nod(grazmicro_n, &
                                    npes, sn, rn, MPI_COMM_FESOM, mype, s_mpitype_nod2D, &
                                    r_mpitype_nod2D, sPE, rPE, requests, nreq)
            call recom_exchange_nod(grazmicro_d, &
                                    npes, sn, rn, MPI_COMM_FESOM, mype, s_mpitype_nod2D, &
                                    r_mpitype_nod2D, sPE, rPE, requests, nreq)
            if (enable_coccos) then
                call recom_exchange_nod(grazmicro_c, &
                                        npes, sn, rn, MPI_COMM_FESOM, mype, s_mpitype_nod2D, &
                                        r_mpitype_nod2D, sPE, rPE, requests, nreq)
                call recom_exchange_nod(grazmicro_p, &
                                        npes, sn, rn, MPI_COMM_FESOM, mype, s_mpitype_nod2D, &
                                        r_mpitype_nod2D, sPE, rPE, requests, nreq)
            endif
        endif
    endif

    do n=1, benthos_num
        call recom_exchange_nod(GlodecayBenthos(:,n), &
                                npes, sn, rn, MPI_COMM_FESOM, mype, s_mpitype_nod2D, &
                                r_mpitype_nod2D, sPE, rPE, requests, nreq)
    end do

    call recom_exchange_nod(GloHplus, &
                            npes, sn, rn, MPI_COMM_FESOM, mype, s_mpitype_nod2D, &
                            r_mpitype_nod2D, sPE, rPE, requests, nreq)
    call recom_exchange_nod(AtmFeInput, &
                            npes, sn, rn, MPI_COMM_FESOM, mype, s_mpitype_nod2D, &
                            r_mpitype_nod2D, sPE, rPE, requests, nreq)
    call recom_exchange_nod(AtmNInput, &
                            npes, sn, rn, MPI_COMM_FESOM, mype, s_mpitype_nod2D, &
                            r_mpitype_nod2D, sPE, rPE, requests, nreq)

    call recom_exchange_nod(PAR3D, &
                            npes, sn, rn, MPI_COMM_FESOM, mype, s_mpitype_nod3D, &
                            r_mpitype_nod3D, sPE, rPE, requests, nreq)

    call recom_exchange_nod(CO23D, &
                            npes, sn, rn, MPI_COMM_FESOM, mype, s_mpitype_nod3D, &
                            r_mpitype_nod3D, sPE, rPE, requests, nreq)
    call recom_exchange_nod(pH3D, &
                            npes, sn, rn, MPI_COMM_FESOM, mype, s_mpitype_nod3D, &
                            r_mpitype_nod3D, sPE, rPE, requests, nreq)
    call recom_exchange_nod(pCO23D, &
                            npes, sn, rn, MPI_COMM_FESOM, mype, s_mpitype_nod3D, &
                            r_mpitype_nod3D, sPE, rPE, requests, nreq)
    call recom_exchange_nod(HCO33D, &
                            npes, sn, rn, MPI_COMM_FESOM, mype, s_mpitype_nod3D, &
                            r_mpitype_nod3D, sPE, rPE, requests, nreq)
    call recom_exchange_nod(CO33D, &
                            npes, sn, rn, MPI_COMM_FESOM, mype, s_mpitype_nod3D, &
                            r_mpitype_nod3D, sPE, rPE, requests, nreq)
    call recom_exchange_nod(OmegaC3D, &
                            npes, sn, rn, MPI_COMM_FESOM, mype, s_mpitype_nod3D, &
                            r_mpitype_nod3D, sPE, rPE, requests, nreq)
    call recom_exchange_nod(kspc3D, &
                            npes, sn, rn, MPI_COMM_FESOM, mype, s_mpitype_nod3D, &
                            r_mpitype_nod3D, sPE, rPE, requests, nreq)
    call recom_exchange_nod(rhoSW3D, &
                            npes, sn, rn, MPI_COMM_FESOM, mype, s_mpitype_nod3D, &
                            r_mpitype_nod3D, sPE, rPE, requests, nreq)

end subroutine recom

! ======================================================================================
! Alkalinity restoring to climatology
! ======================================================================================
subroutine bio_fluxes(alkalinity, MPI_COMM_FESOM, myDim_nod2D, eDim_nod2D, ocean_area, ulevels_nod2D, areasvol)
    use recom_declarations
    use recom_locvar
    use recom_glovar
    use recom_config
    use recom_extra, only: integrate_nod_2d_recom

    implicit none

    integer, intent(in)               :: MPI_COMM_FESOM, myDim_nod2D, eDim_nod2D
    integer, intent(in), dimension(:) :: ulevels_nod2D

    real(kind=WP), intent(in)                 :: ocean_area
    real(kind=WP), intent(in), dimension(:,:) :: areasvol
    real(kind=WP), intent(in), dimension(:,:) :: alkalinity

    integer       :: n, elem, elnodes(3), n1
    real(kind=WP) :: ralk, net

!___________________________________________________________________
! on freshwater inflow/outflow or virtual alkalinity:
  ! 1. In zlevel & zstar the freshwater flux is applied in the update of the
  ! ssh matrix when solving the continuity equation of vertically
  ! integrated flow. The alkalinity concentration in the first layer will
  ! be then adjusted according to the change in volume.

  ! In this case ralk is forced to be zero by setting ref_alk=0. and ref_alk_local=.false.

  ! 2. In cases where the volume of the upper layer is fixed (i.e. linfs)  the freshwater flux
  ! 'ralk*water_flux(n)' is applied as a virtual alkalinity boundary condition via the vertical
  ! diffusion operator.

  ! --> ralk*water_flux(n) : virtual alkalinity flux
  ! virtual alkalinity flux

!  if (use_virt_alk) then ! OG in case of virtual alkalinity flux
!     ralk=ref_alk
!     do n=1, myDim_nod2D+eDim_nod2D
!        if (ref_alk_local) ralk = tracers%data(2+ialk)%values(1, n)
!        virtual_alk(n)=ralk*water_flux(n)
!     end do
!  end if

    !___________________________________________________________________________
    ! Balance alkalinity restoring to climatology
    do n=1, myDim_nod2D + eDim_nod2D
!        relax_alk(n)=surf_relax_Alk * (Alk_surf(n) - tracers%data(2+ialk)%values(1, n))
!        relax_alk(n)=surf_relax_Alk * (Alk_surf(n) - alkalinity(ulevels_nod2d(n),n)
        relax_alk(n)=surf_relax_Alk * (Alk_surf(n) - alkalinity(1, n))
    end do

  ! 2. virtual alkalinity flux
!  if (use_virt_alk) then ! is already zero otherwise
!     call integrate_nod(virtual_alk, net, partit, mesh)
!     virtual_alk=virtual_alk-net/ocean_area
!  end if


  ! 3. restoring to Alkalinity climatology
    call integrate_nod_2D_recom(relax_alk, net, MPI_COMM_FESOM, myDim_nod2D, eDim_nod2D, ulevels_nod2D, areasvol)

    relax_alk = relax_alk - net / ocean_area  ! at ocean surface layer

end subroutine bio_fluxes
