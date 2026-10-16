!===============================================================================
! MODULE: recom4_dump
!   Reference dumps for RECOM4 (https://github.com/RECOM-Regulated-Ecosystem-Model/RECOM4,
!   milestone M1, issue #8). Records, for a few chosen water columns and step
!   windows, everything REcoM receives and returns, so that RECOM4 can be
!   compared with FESOM2-REcoM column by column, bit for bit.
!
!   Switched on only by an optional file `namelist.recom4dump` in the run
!   directory. Without that file nothing happens. The dump only reads model
!   variables; it never changes them.
!
!   namelist.recom4dump, example:
!     &dump_settings
!       dump_dir  = './recom4_dump/'          ! must exist, ends with '/'
!       col_name  = 'tropical', 'polar', 'shelf'
!       col_lon   = -140.0, -40.0, 4.0        ! degrees east
!       col_lat   =    0.0, -65.0, 55.0       ! degrees north
!       win_first = 1, 1081                   ! step windows (mstep), inclusive
!       win_last  = 10, 1090
!     /
!
!   File format
!   -----------
!   One binary stream file per column and recording point:
!     <dump_dir><column>_step.bin       around REcoM_Forcing (recom4_step)
!     <dump_dir><column>_transport.bin  around the transport calls (later step)
!     <dump_dir><column>_surface.bin    surface boundary condition (later step)
!   Each file is a sequence of self-describing field records, appended step
!   by step, in native byte order (little-endian on albedo):
!     name    32 characters, blank padded (e.g. 'in_temp')
!     mstep   int32, the step of this record
!     n1, n2  int32, the shape of the field (scalars: 1, 1)
!     data    n1*n2 float64, column-major (Fortran) order
!   Integers (levels, day of year) are stored as float64; they are exact.
!   Python: name = f.read(32); mstep, n1, n2 = np.fromfile(f, '<i4', 3);
!           data = np.fromfile(f, '<f8', n1*n2).reshape((n1, n2), order='F')
!
!   Step 1: read the namelist, choose the nodes, report them.
!   Step 2: field records around REcoM_Forcing.
!   Step 3: field records around the REcoM part of FESOM's tracer step
!           (benthic return, sinking, sinking into the benthos), per tracer.
!           Per-tracer names carry the tracer ID, e.g. 'in_tr_1001'.
!   Step 4: surface boundary condition per REcoM tracer (bc_surface),
!           split into the terms recom4_step receives, 'surf_<ID>'.
!===============================================================================
module recom4_dump
    use recom_declarations, only: wp
    use recom_config, only: benthos_num
    use recom_glovar, only: Benthos, GlodecayBenthos
    use iso_fortran_env, only: int32, real64
    use mpi
    implicit none
    private

    public :: recom4_dump_init, recom4_dump_step, recom4_dump_column
    public :: recom4_dump_open, recom4_dump_close, recom4_dump_write
    public :: recom4_dump_transport_in, recom4_dump_transport_bf, recom4_dump_transport_sink
    public :: recom4_dump_surface

    ! Recording points
    integer, parameter, public :: DUMP_STEP = 1, DUMP_TRANSPORT = 2, DUMP_SURFACE = 3
    integer, parameter :: nphase = 3
    character(len=*), parameter :: phase_name(nphase) = [character(len=9) :: 'step', 'transport', 'surface']

    interface recom4_dump_write
        module procedure write_r0, write_r1, write_r2, write_i0
    end interface recom4_dump_write

    integer, parameter :: maxcol = 8
    integer, parameter :: maxwin = 4

    logical, save, public :: dump_enabled = .false.

    character(len=512), save :: dump_dir = ''
    integer, save            :: ncol = 0
    integer, save            :: nwin = 0
    character(len=32), save  :: col_name(maxcol) = ''
    real(kind=wp), save      :: col_lon(maxcol) = 0.0_wp
    real(kind=wp), save      :: col_lat(maxcol) = 0.0_wp
    integer, save            :: win_first(maxwin) = 0
    integer, save            :: win_last(maxwin) = 0

    ! Per column: owning MPI rank and local node index on that rank
    integer, save :: col_owner(maxcol) = -1
    integer, save :: col_node(maxcol) = 0
    integer, save :: my_rank = -1

    ! File currently open for writing, and whether a file was started yet
    integer, save :: cur_unit = -1
    integer, save :: cur_step = 0
    logical, save :: started(maxcol, nphase) = .false.
    ! Last step for which the per-step records of a column were written
    integer, save :: last_step(maxcol, nphase) = -1

contains

    !===========================================================================
    ! recom4_dump_init
    !   Called once at the end of recom_init. Reads namelist.recom4dump if it
    !   exists and chooses, per column, the open-ocean node (no cavity) nearest
    !   to the given coordinates over all ranks.
    !===========================================================================
    subroutine recom4_dump_init(ulevels_nod2D, nlevels_nod2D, geo_coord_nod2D, &
            myDim_nod2D, mype, comm)
        integer, intent(in)       :: ulevels_nod2D(:), nlevels_nod2D(:)
        real(kind=wp), intent(in) :: geo_coord_nod2D(:, :)   ! (lon, lat) in radians
        integer, intent(in)       :: myDim_nod2D, mype, comm

        namelist /dump_settings/ dump_dir, col_name, col_lon, col_lat, win_first, win_last

        character(len=*), parameter :: nml_file = 'namelist.recom4dump'
        real(kind=wp), parameter :: deg = 180.0_wp / (4.0_wp * atan(1.0_wp))
        real(kind=wp), parameter :: earth_radius_km = 6371.0_wp

        logical :: exists
        integer :: unit, ios, ierr, i, n
        real(kind=wp) :: lon0, lat0, dist, dmin
        real(kind=wp) :: pair_in(2), pair_out(2)
        integer :: nbest

        my_rank = mype
        dump_enabled = .false.

        inquire(file=nml_file, exist=exists)
        if (.not. exists) return

        open(newunit=unit, file=nml_file, status='old', action='read', iostat=ios)
        if (ios == 0) read(unit, nml=dump_settings, iostat=ios)
        if (ios /= 0) then
            if (mype == 0) write(*, '(a,i0)') 'recom4_dump: ERROR reading ' // nml_file // ', iostat = ', ios
            call MPI_Abort(comm, 1, ierr)
        end if
        close(unit)

        ncol = count(len_trim(col_name) > 0)
        nwin = count(win_first > 0)
        if (ncol == 0 .or. nwin == 0) then
            if (mype == 0) write(*, '(a)') 'recom4_dump: no columns or no step windows, dump off'
            return
        end if

        do i = 1, ncol
            lon0 = col_lon(i) / deg
            lat0 = col_lat(i) / deg

            dmin = huge(1.0_wp)
            nbest = 0
            do n = 1, myDim_nod2D
                if (ulevels_nod2D(n) /= 1) cycle              ! cavity node
                if (nlevels_nod2D(n) - 1 < 1) cycle           ! no wet level
                dist = great_circle(lon0, lat0, geo_coord_nod2D(1, n), geo_coord_nod2D(2, n))
                if (dist < dmin) then
                    dmin = dist
                    nbest = n
                end if
            end do

            ! Overall nearest node; on equal distance the lowest rank wins
            pair_in(1) = dmin
            pair_in(2) = real(mype, wp)
            call MPI_Allreduce(pair_in, pair_out, 1, MPI_2DOUBLE_PRECISION, MPI_MINLOC, comm, ierr)
            col_owner(i) = nint(pair_out(2))

            if (mype == col_owner(i)) then
                col_node(i) = nbest
                write(*, '(a,a,a,f10.4,a,f9.4,a,i0,a,f10.3,a,i0,a,i0)') &
                    'recom4_dump: column ', trim(col_name(i)), &
                    ' -> node at lon ', geo_coord_nod2D(1, nbest) * deg, &
                    ' lat ', geo_coord_nod2D(2, nbest) * deg, &
                    ', levels ', nlevels_nod2D(nbest) - 1, &
                    ', distance ', dmin * earth_radius_km, ' km, rank ', mype, &
                    ', local node ', nbest
            end if
        end do

        dump_enabled = .true.
        if (mype == 0) then
            write(*, '(a,a)') 'recom4_dump: output directory ', trim(dump_dir)
            do i = 1, nwin
                write(*, '(a,i0,a,i0)') 'recom4_dump: steps ', win_first(i), ' to ', win_last(i)
            end do
        end if
    end subroutine recom4_dump_init

    !===========================================================================
    ! recom4_dump_step: .true. if step mstep lies in one of the windows
    !===========================================================================
    logical function recom4_dump_step(mstep)
        integer, intent(in) :: mstep
        integer :: w

        recom4_dump_step = .false.
        if (.not. dump_enabled) return
        do w = 1, nwin
            if (mstep >= win_first(w) .and. mstep <= win_last(w)) recom4_dump_step = .true.
        end do
    end function recom4_dump_step

    !===========================================================================
    ! recom4_dump_column: index of the column whose node is local node n on
    ! this rank, 0 if none
    !===========================================================================
    integer function recom4_dump_column(n)
        integer, intent(in) :: n
        integer :: i

        recom4_dump_column = 0
        if (.not. dump_enabled) return
        do i = 1, ncol
            if (col_owner(i) == my_rank .and. col_node(i) == n) then
                recom4_dump_column = i
                return
            end if
        end do
    end function recom4_dump_column

    !===========================================================================
    ! recom4_dump_open: start writing records of step mstep for column icol at
    ! recording point iphase. The first call per column and point creates the
    ! file (replacing an old one); later calls append.
    !===========================================================================
    subroutine recom4_dump_open(icol, iphase, mstep)
        integer, intent(in) :: icol, iphase, mstep

        character(len=1024) :: path
        integer :: ios

        path = trim(dump_dir) // trim(col_name(icol)) // '_' // trim(phase_name(iphase)) // '.bin'
        if (.not. started(icol, iphase)) then
            open(newunit=cur_unit, file=path, access='stream', form='unformatted', &
                 status='replace', action='write', iostat=ios)
            started(icol, iphase) = .true.
            write(*, '(a,a,a,i0)') 'recom4_dump: writing ', trim(path), ' from step ', mstep
        else
            open(newunit=cur_unit, file=path, access='stream', form='unformatted', &
                 status='old', position='append', action='write', iostat=ios)
        end if
        if (ios /= 0) then
            write(*, '(a,a,a,i0)') 'recom4_dump: ERROR opening ', trim(path), ', iostat = ', ios
            write(*, '(a)') 'recom4_dump: does dump_dir exist and end with / ?'
            call MPI_Abort(MPI_COMM_WORLD, 1, ios)
        end if
        cur_step = mstep
    end subroutine recom4_dump_open

    subroutine recom4_dump_close()
        if (cur_unit >= 0) close(cur_unit)
        cur_unit = -1
    end subroutine recom4_dump_close

    !===========================================================================
    ! Field records: name, step, shape, float64 data (see the file format above)
    !===========================================================================
    subroutine write_header(name, n1, n2)
        character(len=*), intent(in) :: name
        integer, intent(in) :: n1, n2
        character(len=32) :: tag

        tag = name
        write(cur_unit) tag, int(cur_step, int32), int(n1, int32), int(n2, int32)
    end subroutine write_header

    subroutine write_r0(name, x)
        character(len=*), intent(in) :: name
        real(kind=wp), intent(in) :: x
        call write_header(name, 1, 1)
        write(cur_unit) real(x, real64)
    end subroutine write_r0

    subroutine write_i0(name, i)
        character(len=*), intent(in) :: name
        integer, intent(in) :: i
        call write_header(name, 1, 1)
        write(cur_unit) real(i, real64)
    end subroutine write_i0

    subroutine write_r1(name, a)
        character(len=*), intent(in) :: name
        real(kind=wp), intent(in) :: a(:)
        call write_header(name, size(a), 1)
        write(cur_unit) real(a, real64)
    end subroutine write_r1

    subroutine write_r2(name, a)
        character(len=*), intent(in) :: name
        real(kind=wp), intent(in) :: a(:, :)
        call write_header(name, size(a, 1), size(a, 2))
        write(cur_unit) real(a, real64)
    end subroutine write_r2

    !===========================================================================
    ! Step 3: the REcoM part of FESOM's tracer step (diff_tracers_ale).
    !   Called for REcoM tracers only, with FESOM's whole arrays (level, node);
    !   each routine picks the dumped nodes owned by this rank. No MPI calls,
    !   so ranks without a dumped node simply return.
    !
    !   recom4_dump_transport_in   before the REcoM block. On the first REcoM
    !       tracer of a step it also writes the column geometry, temperature,
    !       salinity, the benthos and its decay (what recom4_transport gets
    !       besides the tracer). Then in_tr_<ID>.
    !   recom4_dump_transport_bf   after the benthic return was added:
    !       out_dtr_bf_<ID>, out_tr_bf_<ID>
    !   recom4_dump_transport_sink after sinking was added:
    !       out_vert_sink_<ID>, out_str_bf_<ID>, out_tr_sink_<ID>,
    !       out_benthos_<ID> (sinking into the benthos changes it per tracer)
    !===========================================================================
    subroutine recom4_dump_transport_in(tr_id, mstep, dt, ulevels_nod2D, nlevels_nod2D, &
            nod_in_elem2D_num, nod_in_elem2D, nlevels, elem_area, area, areasvol, &
            hnode, hnode_new, zbar_3d_n, z_3d_n, geo_coord_nod2D, temp, salt, tr)
        integer, intent(in)       :: tr_id, mstep
        real(kind=wp), intent(in) :: dt
        integer, intent(in)       :: ulevels_nod2D(:), nlevels_nod2D(:)
        integer, intent(in)       :: nod_in_elem2D_num(:), nod_in_elem2D(:, :), nlevels(:)
        real(kind=wp), intent(in) :: elem_area(:), area(:, :), areasvol(:, :)
        real(kind=wp), intent(in) :: hnode(:, :), hnode_new(:, :), zbar_3d_n(:, :), z_3d_n(:, :)
        real(kind=wp), intent(in) :: geo_coord_nod2D(:, :), temp(:, :), salt(:, :), tr(:, :)

        integer :: icol, n, nzmin, nzmax, j, e
        real(kind=wp) :: seafloor(size(area, 1))

        if (.not. recom4_dump_step(mstep)) return
        do icol = 1, ncol
            n = local_node(icol)
            if (n == 0) cycle
            nzmin = ulevels_nod2D(n)
            nzmax = nlevels_nod2D(n) - 1
            call recom4_dump_open(icol, DUMP_TRANSPORT, mstep)
            if (last_step(icol, DUMP_TRANSPORT) /= mstep) then
                last_step(icol, DUMP_TRANSPORT) = mstep
                ! seafloor area under each layer, as in diff_ver_recom_expl
                seafloor = 0.0_wp
                do j = 1, nod_in_elem2D_num(n)
                    e = nod_in_elem2D(j, n)
                    seafloor(nlevels(e) - 1) = seafloor(nlevels(e) - 1) + elem_area(e) / 3.0_wp
                end do
                call recom4_dump_write('in_nzmin', nzmin)
                call recom4_dump_write('in_nzmax', nzmax)
                call recom4_dump_write('in_mstep', mstep)
                call recom4_dump_write('in_dt', dt)
                call recom4_dump_write('in_lat', geo_coord_nod2D(2, n))
                call recom4_dump_write('in_hnode', hnode(nzmin:nzmax, n))
                call recom4_dump_write('in_hnode_new', hnode_new(nzmin:nzmax, n))
                call recom4_dump_write('in_zbar', zbar_3d_n(nzmin:nzmax + 1, n))
                call recom4_dump_write('in_zmid', z_3d_n(nzmin:nzmax, n))
                call recom4_dump_write('in_area', area(nzmin:nzmax + 1, n))
                call recom4_dump_write('in_areasvol', areasvol(nzmin:nzmax, n))
                call recom4_dump_write('in_seafloor', seafloor(nzmin:nzmax))
                call recom4_dump_write('in_temp', temp(nzmin:nzmax, n))
                call recom4_dump_write('in_salt', salt(nzmin:nzmax, n))
                call recom4_dump_write('in_benthos', Benthos(n, 1:benthos_num))
                call recom4_dump_write('in_decay_benthos', GlodecayBenthos(n, 1:benthos_num))
            end if
            call recom4_dump_write(tagged('in_tr', tr_id), tr(nzmin:nzmax, n))
            call recom4_dump_close()
        end do
    end subroutine recom4_dump_transport_in

    subroutine recom4_dump_transport_bf(tr_id, mstep, ulevels_nod2D, nlevels_nod2D, dtr_bf, tr)
        integer, intent(in)       :: tr_id, mstep
        integer, intent(in)       :: ulevels_nod2D(:), nlevels_nod2D(:)
        real(kind=wp), intent(in) :: dtr_bf(:, :), tr(:, :)

        integer :: icol, n, nzmin, nzmax

        if (.not. recom4_dump_step(mstep)) return
        do icol = 1, ncol
            n = local_node(icol)
            if (n == 0) cycle
            nzmin = ulevels_nod2D(n)
            nzmax = nlevels_nod2D(n) - 1
            call recom4_dump_open(icol, DUMP_TRANSPORT, mstep)
            call recom4_dump_write(tagged('out_dtr_bf', tr_id), dtr_bf(nzmin:nzmax, n))
            call recom4_dump_write(tagged('out_tr_bf', tr_id), tr(nzmin:nzmax, n))
            call recom4_dump_close()
        end do
    end subroutine recom4_dump_transport_bf

    subroutine recom4_dump_transport_sink(tr_id, mstep, ulevels_nod2D, nlevels_nod2D, &
            vert_sink, str_bf, tr)
        integer, intent(in)       :: tr_id, mstep
        integer, intent(in)       :: ulevels_nod2D(:), nlevels_nod2D(:)
        real(kind=wp), intent(in) :: vert_sink(:, :), str_bf(:, :), tr(:, :)

        integer :: icol, n, nzmin, nzmax

        if (.not. recom4_dump_step(mstep)) return
        do icol = 1, ncol
            n = local_node(icol)
            if (n == 0) cycle
            nzmin = ulevels_nod2D(n)
            nzmax = nlevels_nod2D(n) - 1
            call recom4_dump_open(icol, DUMP_TRANSPORT, mstep)
            call recom4_dump_write(tagged('out_vert_sink', tr_id), vert_sink(nzmin:nzmax, n))
            call recom4_dump_write(tagged('out_str_bf', tr_id), str_bf(nzmin:nzmax, n))
            call recom4_dump_write(tagged('out_tr_sink', tr_id), tr(nzmin:nzmax, n))
            call recom4_dump_write(tagged('out_benthos', tr_id), Benthos(n, 1:benthos_num))
            call recom4_dump_close()
        end do
    end subroutine recom4_dump_transport_sink

    !===========================================================================
    ! Step 4: surface boundary condition of one REcoM tracer at node n
    !   terms = own, virtual, restoring, river, erosion (each per second, the
    !   factors is_riverinput / is_erosioninput / useRivFe already applied),
    !   then FESOM's bc_surface (= dt * their sum). Record 'surf_<ID>', 6 values.
    !===========================================================================
    subroutine recom4_dump_surface(n, mstep, tr_id, terms)
        integer, intent(in)       :: n, mstep, tr_id
        real(kind=wp), intent(in) :: terms(6)

        integer :: icol

        if (.not. dump_enabled) return
        icol = recom4_dump_column(n)
        if (icol == 0) return
        if (.not. recom4_dump_step(mstep)) return
        call recom4_dump_open(icol, DUMP_SURFACE, mstep)
        call recom4_dump_write(tagged('surf', tr_id), terms)
        call recom4_dump_close()
    end subroutine recom4_dump_surface

    !===========================================================================
    ! local_node: local node index of column icol if this rank owns it, else 0
    !===========================================================================
    integer function local_node(icol)
        integer, intent(in) :: icol
        local_node = 0
        if (col_owner(icol) == my_rank) local_node = col_node(icol)
    end function local_node

    !===========================================================================
    ! tagged: record name with the tracer ID appended, e.g. 'in_tr_1001'
    !===========================================================================
    function tagged(base, id) result(name)
        character(len=*), intent(in) :: base
        integer, intent(in) :: id
        character(len=32) :: name
        character(len=12) :: num
        write(num, '(i0)') id
        name = trim(base) // '_' // trim(num)
    end function tagged

    !===========================================================================
    ! Great-circle angle between two points (radians in, radians out)
    !===========================================================================
    pure real(kind=wp) function great_circle(lon1, lat1, lon2, lat2)
        real(kind=wp), intent(in) :: lon1, lat1, lon2, lat2
        real(kind=wp) :: a

        a = sin(0.5_wp * (lat2 - lat1))**2 &
            + cos(lat1) * cos(lat2) * sin(0.5_wp * (lon2 - lon1))**2
        great_circle = 2.0_wp * asin(sqrt(min(1.0_wp, a)))
    end function great_circle

end module recom4_dump
