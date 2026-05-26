Module modio

    Use modconst
    Use hdf5

    Integer(kind=4)   :: mynpart       ! local particle number
    Integer(kind=4)   :: ndim          ! number of dimensions (always 3)
    Integer(kind=4)   :: lmin          ! log2(particle_grid_nx) - kept for compatibility
    Integer(kind=4)   :: nres          ! particle grid resolution in each dimension (= particle_grid_nx)
    Integer(kind=PRI) :: nptot         ! total DM particle number (post star filtering)
    Integer(kind=PRI) :: ngrid         ! total number of grid points: ngrid = nres ** 3
    Real   (kind=SP)  :: xmin, xmax, ymin, ymax, zmin, zmax  ! min and max (x,y,z) for each process (normalized [0,1])
    Integer(kind=4)   :: nstar         ! number of star particles on this rank after splitting
    Integer(kind=PRI) :: nstot         ! total number of star particles across ranks

    ! Snapshot physical metadata (read from .ini)
    Real(kind=DP)     :: Lbox          ! physical box size (assumed cubic — uses x range)
    Real(kind=DP)     :: phys_xmin, phys_xmax  ! mesh.xmin / xmax in physical units
    Integer(kind=4)   :: particle_nx   ! particle_grid.nx from .ini

    Contains


    ! Build an output path under outdir. Returns trim(outdir)//'/'//trim(name),
    ! unless outdir is empty or '.', in which case returns trim(name).
    ! A trailing slash on outdir is tolerated.
    Function out_path(name) Result(p)
        Use modparam, Only: outdir
        Character(len=*), intent(in) :: name
        Character(len=1200) :: p
        Integer :: n
        n = len_trim(outdir)
        If (n == 0 .or. trim(outdir) == '.') Then
            p = trim(name)
        Else If (outdir(n:n) == '/') Then
            p = trim(outdir)//trim(name)
        Else
            p = trim(outdir)//'/'//trim(name)
        End If
    End Function out_path


    ! Write input parameters in file .inp
    Subroutine outputparameters()
        Use modparam
        Use modmpicom
        Implicit none

        Character(len=1200) :: fileopa

        fileopa = out_path(trim(root)//'.inp')
        Open(Unit=Uopa,file=fileopa)
        Write(Uopa,*) 'Nb of processes:'
        Write(Uopa,*) procNB
        Write(Uopa,*) ''
        Write(Uopa,*) 'Input parameters:'
        Write(Uopa,*) root
        Write(Uopa,*) outdir
        Write(Uopa,*) code_index
        Write(Uopa,*) pathinput
        Write(Uopa,*) namepart
        Write(Uopa,*) nameinfo
        Write(Uopa,*) grpsize
        Write(Uopa,*) perco
        Write(Uopa,*) Mmin
        Write(Uopa,*) Mmax
        Write(Uopa,*) star
        Write(Uopa,*) metal
        Write(Uopa,*) outcube
        Write(Uopa,*) dofof
        Write(Uopa,*) readfromcube
        Write(Uopa,*) dotimings
        Write(Uopa,*) invent_ids

        Close(Uopa)

    End Subroutine outputparameters


    ! Minimal parser for a dyablo .ini config file.
    ! Extracts the [mesh] bounds and [particle_grid].nx that we need to normalize
    ! coordinates and set the FoF search grid resolution.
    Subroutine read_dyablo_ini(filename, mxmin, mxmax, mymin, mymax, mzmin, mzmax, nx_out, ok)

        Implicit none

        Character(len=*),  Intent(in)  :: filename
        Real(kind=DP),     Intent(out) :: mxmin, mxmax, mymin, mymax, mzmin, mzmax
        Integer(kind=4),   Intent(out) :: nx_out
        Logical,           Intent(out) :: ok

        Character(len=512) :: line, raw, section, key, val
        Integer :: ios, eq_pos, brk_pos, semi_pos, hash_pos, u, n
        Integer(kind=4) :: amr_bx, amr_level_min
        Logical :: have_xmin, have_xmax, have_ymin, have_ymax, have_zmin, have_zmax, have_nx
        Logical :: have_bx, have_level_min

        ok = .false.
        have_xmin = .false.; have_xmax = .false.
        have_ymin = .false.; have_ymax = .false.
        have_zmin = .false.; have_zmax = .false.
        have_nx   = .false.
        have_bx   = .false.; have_level_min = .false.
        nx_out = -1
        amr_bx = -1; amr_level_min = -1
        mxmin = 0.d0; mxmax = 1.d0
        mymin = 0.d0; mymax = 1.d0
        mzmin = 0.d0; mzmax = 1.d0

        Open(unit=77, file=filename, status='old', action='read', iostat=ios)
        If (ios /= 0) Return
        u = 77

        section = ''

        Do
            Read(u, '(A)', iostat=ios) raw
            If (ios /= 0) Exit

            line = raw
            ! strip everything after ; or # (inline comments)
            semi_pos = index(line, ';')
            If (semi_pos > 0) line = line(1:semi_pos-1)
            hash_pos = index(line, '#')
            If (hash_pos > 0) line = line(1:hash_pos-1)

            line = adjustl(line)
            n = len_trim(line)
            If (n == 0) Cycle

            ! section header [name]
            If (line(1:1) == '[') Then
                brk_pos = index(line, ']')
                If (brk_pos > 2) section = line(2:brk_pos-1)
                Cycle
            End If

            ! key = value
            eq_pos = index(line, '=')
            If (eq_pos == 0) Cycle

            key = adjustl(line(1:eq_pos-1))
            val = adjustl(line(eq_pos+1:))
            ! trim trailing whitespace inside key/val
            key = trim(key)
            val = trim(val)

            Select Case (trim(section))
            Case ('mesh')
                Select Case (trim(key))
                Case ('xmin'); Read(val, *, iostat=ios) mxmin; If (ios == 0) have_xmin = .true.
                Case ('xmax'); Read(val, *, iostat=ios) mxmax; If (ios == 0) have_xmax = .true.
                Case ('ymin'); Read(val, *, iostat=ios) mymin; If (ios == 0) have_ymin = .true.
                Case ('ymax'); Read(val, *, iostat=ios) mymax; If (ios == 0) have_ymax = .true.
                Case ('zmin'); Read(val, *, iostat=ios) mzmin; If (ios == 0) have_zmin = .true.
                Case ('zmax'); Read(val, *, iostat=ios) mzmax; If (ios == 0) have_zmax = .true.
                End Select
            Case ('amr')
                ! Fallback: derive particle_grid nx from bx * 2^level_min
                Select Case (trim(key))
                Case ('bx')
                    Read(val, *, iostat=ios) amr_bx
                    If (ios == 0) have_bx = .true.
                Case ('level_min')
                    Read(val, *, iostat=ios) amr_level_min
                    If (ios == 0) have_level_min = .true.
                End Select
            Case ('particle_grid')
                If (trim(key) == 'nx') Then
                    Read(val, *, iostat=ios) nx_out
                    If (ios == 0) have_nx = .true.
                End If
            End Select
        End Do

        Close(u)

        ! If [particle_grid].nx was absent, derive it from [amr].bx * 2^level_min
        If (.not. have_nx .and. have_bx .and. have_level_min) Then
            nx_out = amr_bx * (2 ** amr_level_min)
            have_nx = .true.
            Print *, '   [particle_grid].nx not found; derived from [amr]: bx=', amr_bx, &
                     ' level_min=', amr_level_min, ' -> nx=', nx_out
        End If

        ok = have_xmin .and. have_xmax .and. have_ymin .and. have_ymax .and. &
             have_zmin .and. have_zmax .and. have_nx

    End Subroutine read_dyablo_ini


    ! Read a 1D float64 dataset slice [offset0, offset0+nread) from an open file
    Subroutine read_h5_slice_dp1(file_id, dsetname, offset0, nread, buf, ok)

        Implicit none

        Integer(HID_T),    Intent(in)  :: file_id
        Character(len=*),  Intent(in)  :: dsetname
        Integer(HSIZE_T),  Intent(in)  :: offset0, nread
        Real(kind=DP),     Intent(out) :: buf(nread)
        Logical,           Intent(out) :: ok

        Integer(HID_T)   :: dset_id, file_space, mem_space
        Integer(HSIZE_T) :: offset(1), counts(1), mdims(1)
        Integer          :: hdferr

        ok = .false.
        Call h5dopen_f(file_id, dsetname, dset_id, hdferr)
        If (hdferr /= 0) Return
        Call h5dget_space_f(dset_id, file_space, hdferr)
        offset(1) = offset0
        counts(1) = nread
        Call h5sselect_hyperslab_f(file_space, H5S_SELECT_SET_F, offset, counts, hdferr)
        mdims(1) = nread
        Call h5screate_simple_f(1, mdims, mem_space, hdferr)
        Call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, buf, mdims, hdferr, mem_space, file_space)
        Call h5sclose_f(mem_space, hdferr)
        Call h5sclose_f(file_space, hdferr)
        Call h5dclose_f(dset_id, hdferr)
        ok = (hdferr == 0)

    End Subroutine read_h5_slice_dp1


    ! Read a hyperslab from a 2D float64 dataset of shape (N, 3) on disk.
    ! Fortran sees it as (3, N) due to column-major storage, so we read
    ! a (3, nread) hyperslab.
    Subroutine read_h5_slice_dp2_coord(file_id, dsetname, offset0, nread, buf, ok)

        Implicit none

        Integer(HID_T),    Intent(in)  :: file_id
        Character(len=*),  Intent(in)  :: dsetname
        Integer(HSIZE_T),  Intent(in)  :: offset0, nread
        Real(kind=DP),     Intent(out) :: buf(3, nread)
        Logical,           Intent(out) :: ok

        Integer(HID_T)   :: dset_id, file_space, mem_space
        Integer(HSIZE_T) :: offset(2), counts(2), mdims(2)
        Integer          :: hdferr

        ok = .false.
        Call h5dopen_f(file_id, dsetname, dset_id, hdferr)
        If (hdferr /= 0) Return
        Call h5dget_space_f(dset_id, file_space, hdferr)
        ! Fortran indices: dim 1 = fast = (x,y,z), dim 2 = slow = particle index
        offset(1) = 0
        offset(2) = offset0
        counts(1) = 3
        counts(2) = nread
        Call h5sselect_hyperslab_f(file_space, H5S_SELECT_SET_F, offset, counts, hdferr)
        mdims(1) = 3
        mdims(2) = nread
        Call h5screate_simple_f(2, mdims, mem_space, hdferr)
        Call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, buf, mdims, hdferr, mem_space, file_space)
        Call h5sclose_f(mem_space, hdferr)
        Call h5sclose_f(file_space, hdferr)
        Call h5dclose_f(dset_id, hdferr)
        ok = (hdferr == 0)

    End Subroutine read_h5_slice_dp2_coord


    ! Quickselect-based median for a Real(SP) array (modifies arr).
    Real(kind=SP) Function median_inplace_sp(arr, n) Result(med)
        Implicit none
        Integer(kind=4), Intent(in) :: n
        Real(kind=SP),   Intent(inout) :: arr(n)
        Integer :: lo, hi, mid, i, j
        Real(kind=SP) :: pivot, tmp

        lo = 1
        hi = n
        mid = (n + 1) / 2

        Do While (lo < hi)
            pivot = arr((lo + hi) / 2)
            i = lo
            j = hi
            Do
                Do While (arr(i) < pivot)
                    i = i + 1
                End Do
                Do While (arr(j) > pivot)
                    j = j - 1
                End Do
                If (i <= j) Then
                    tmp = arr(i); arr(i) = arr(j); arr(j) = tmp
                    i = i + 1
                    j = j - 1
                End If
                If (i > j) Exit
            End Do
            If (j < mid) lo = i
            If (i > mid) hi = j
        End Do

        med = arr(mid)
    End Function median_inplace_sp


    ! Read particles from a dyablo HDF5 snapshot, filter DM from stars,
    ! invent globally-unique sequential IDs, and distribute geographically
    ! across the MPI cartesian process grid.
    !
    ! Replaces the original ramses_lecture; downstream modules (modfofpara,
    ! modio.outputstruct, etc.) see exactly the same arrays as before:
    !   x(3,mynpart), v(3,mynpart), id(mynpart), id4(mynpart), mp(mynpart),
    !   tp(mynpart) [= 0 for DM], plus the star-only xstar/vstar/idstar/...
    Subroutine dyablo_lecture()
        Use modvariable
        Use modparam
        Use modmpicom
        Use modtiming
        Implicit none

        Character(len=1200)             :: nomfich, ininom
        Character(len=9)               :: tmpstr1, tmpstr2

        Integer(kind=4)                :: i, j, idim
        Integer(kind=4)                :: destCoord(3)
        Integer(kind=4)                :: nrecv, recvpoint
        Integer(kind=4)                :: nsd, n_i, n_j, n_k, ind
        Integer(kind=4)                :: prov, dest
        Integer(kind=4)                :: mpistat(Mpi_Status_Size)
        Integer(kind=4)                :: errcode
        Integer(kind=4)                :: ipart, ndm, istar, idm
        Integer(kind=4)                :: nstar_tot
        Integer(kind=PRI)              :: tmplongint, idmin, idpart
        Integer(kind=4), allocatable   :: npartv(:), npartvloc(:)
        Integer(kind=PRI), allocatable :: tabtmp(:), tabnpart(:)
        Real(kind=SP)                  :: deltasd

        ! HDF5 handles / state
        Integer(HID_T)                 :: file_id
        Integer(HID_T)                 :: dset_id, file_space
        Integer(HSIZE_T)               :: file_dims(1), maxdims(1), file_offset, nread
        Integer                        :: hdferr
        Logical                        :: ok, has_id_dataset
        Integer(kind=PRI)              :: N_total_in_file

        ! Buffers for the local chunk read from disk
        Real(kind=DP), allocatable     :: dcoords(:,:)
        Real(kind=DP), allocatable     :: dvx(:), dvy(:), dvz(:), dmass(:), dbirth(:)
        Integer(kind=PRI4), allocatable :: id4_from_file(:)

        ! Temporary distribution buffers (analog of the RAMSES temporaries)
        Real(kind=SP), allocatable     :: tmpx(:,:), tmpv(:,:)
        Real(kind=SP), allocatable     :: tmpmp(:), tmptp(:)
        Integer(kind=PRI), allocatable :: tmpi(:)
        Integer(kind=PRI4), allocatable :: tmpi4(:)
        Real(kind=SP), allocatable     :: tmpsendx(:,:), tmpsendv(:,:)
        Real(kind=SP), allocatable     :: tmpsendmp(:), tmpsendtp(:)
        Integer(kind=PRI), allocatable :: tmpsendi(:)
        Integer(kind=PRI4), allocatable :: tmpsendi4(:)

        ! Mass-median scratch
        Real(kind=SP), allocatable     :: mass_sample(:)
        Real(kind=SP)                  :: local_median, global_median
        Real(kind=SP), allocatable     :: per_rank_med(:)


        ! ------------------------------------------------------------------
        ! 1. Process 0 reads the .ini config to get Lbox and particle_grid_nx.
        ! ------------------------------------------------------------------
        time0 = Mpi_Wtime()

        ! Broadcast a 4-int header: ndim, particle_nx, nres, (placeholder)
        h_length = 4 * bit_size(particle_nx) / 8
        Allocate(header(0:h_length-1))
        h_pos = 0

        If (procID == 0) Then
            Block
                Real(kind=DP) :: ax, bx, ay, by, az, bz
                If (len_trim(pathinput) > 0) Then
                    ininom = trim(pathinput)//'/'//trim(nameinfo)
                Else
                    ininom = trim(nameinfo)
                End If
                Print *, 'Reading dyablo .ini config: ', trim(ininom)
                Call read_dyablo_ini(ininom, ax, bx, ay, by, az, bz, particle_nx, ok)
                If (.not. ok) Then
                    Print *, '** Failed to parse dyablo .ini file: ', trim(ininom)
                    Print *, '   Need [mesh].x/y/zmin,max and [particle_grid].nx **'
                    errcode = 3
                    Call Mpi_Abort(Mpi_Comm_World, errcode, mpierr)
                End If
                phys_xmin = ax
                phys_xmax = bx
                Lbox = bx - ax
                If (abs((by - ay) - Lbox) > 1.d-6 * Lbox .or. &
                    abs((bz - az) - Lbox) > 1.d-6 * Lbox) Then
                    Print *, '** WARNING: non-cubic mesh; using x-range as Lbox =', Lbox
                End If
                nres = particle_nx
                lmin = nint(log(real(particle_nx, kind=DP)) / log(2.d0))
                ndim = 3
                Write(*, '(A, ES20.10)') ' Lbox (physical units):     ', Lbox
                Write(*, '(A, I12)')     ' particle_grid.nx (nres):   ', particle_nx
                Write(*, '(A, I12)')     ' lmin (= log2 particle_nx): ', lmin
                Write(Ulog,*) 'Lbox = ', Lbox, ' particle_nx = ', particle_nx, ' lmin = ', lmin
            End Block

            Call Mpi_Pack(ndim,        1, Mpi_Integer, header, h_length, h_pos, Mpi_Comm_World, mpierr)
            Call Mpi_Pack(particle_nx, 1, Mpi_Integer, header, h_length, h_pos, Mpi_Comm_World, mpierr)
            Call Mpi_Pack(nres,        1, Mpi_Integer, header, h_length, h_pos, Mpi_Comm_World, mpierr)
            Call Mpi_Pack(lmin,        1, Mpi_Integer, header, h_length, h_pos, Mpi_Comm_World, mpierr)
        End If

        Call Mpi_Bcast(header, h_length, Mpi_Packed, 0, Mpi_Comm_World, mpierr)
        If (procID /= 0) Then
            Call Mpi_Unpack(header, h_length, h_pos, ndim,        1, Mpi_Integer, Mpi_Comm_World, mpierr)
            Call Mpi_Unpack(header, h_length, h_pos, particle_nx, 1, Mpi_Integer, Mpi_Comm_World, mpierr)
            Call Mpi_Unpack(header, h_length, h_pos, nres,        1, Mpi_Integer, Mpi_Comm_World, mpierr)
            Call Mpi_Unpack(header, h_length, h_pos, lmin,        1, Mpi_Integer, Mpi_Comm_World, mpierr)
        End If
        Call Mpi_Bcast(Lbox,      1, Mpi_Double_Precision, 0, Mpi_Comm_World, mpierr)
        Call Mpi_Bcast(phys_xmin, 1, Mpi_Double_Precision, 0, Mpi_Comm_World, mpierr)
        Deallocate(header)

        ngrid = int(nres, kind=8) ** 3

        ! ------------------------------------------------------------------
        ! 2. Open the HDF5 snapshot, determine N_total, divide work.
        ! ------------------------------------------------------------------
        If (len_trim(pathinput) > 0) Then
            nomfich = trim(pathinput)//'/'//trim(namepart)
        Else
            nomfich = trim(namepart)
        End If
        If (procID == 0) Print *, 'Opening dyablo HDF5 particle file: ', trim(nomfich)

        Call h5open_f(hdferr)
        Call h5fopen_f(nomfich, H5F_ACC_RDONLY_F, file_id, hdferr)
        If (hdferr /= 0) Call EmergencyStop('Cannot open HDF5 file '//trim(nomfich), 4)

        Call h5dopen_f(file_id, 'mass', dset_id, hdferr)
        Call h5dget_space_f(dset_id, file_space, hdferr)
        Call h5sget_simple_extent_dims_f(file_space, file_dims, maxdims, hdferr)
        Call h5sclose_f(file_space, hdferr)
        Call h5dclose_f(dset_id, hdferr)
        N_total_in_file = int(file_dims(1), kind=PRI)

        ! Check for /id dataset (dyablo currently does not write it).
        Call h5lexists_f(file_id, 'id', has_id_dataset, hdferr)
        If (procID == 0) Then
            Write(*, '(A, I0)') ' Particles in file: ', N_total_in_file
            If (has_id_dataset) Then
                Print *, ' /id dataset present in snapshot.'
            Else
                Print *, ' /id dataset NOT present in snapshot - will invent IDs.'
            End If
            If (invent_ids .and. has_id_dataset) Then
                Print *, ' invent_ids=.true. requested - ignoring snapshot /id dataset.'
            End If
        End If

        ! Split N_total across procNB ranks: each rank reads [file_offset, file_offset+nread)
        Block
            Integer(kind=PRI) :: base, rem
            base = N_total_in_file / procNB
            rem  = N_total_in_file - base * procNB
            If (int(procID, kind=PRI) < rem) Then
                file_offset = base * int(procID, kind=PRI) + int(procID, kind=PRI)
                nread       = base + 1
            Else
                file_offset = base * int(procID, kind=PRI) + rem
                nread       = base
            End If
        End Block

        ! ------------------------------------------------------------------
        ! 3. Read this rank's slice: coords (DP), vx/vy/vz (DP), mass, birth_time.
        ! ------------------------------------------------------------------
        Allocate(dcoords(3, nread))
        Allocate(dvx(nread), dvy(nread), dvz(nread))
        Allocate(dmass(nread), dbirth(nread))

        Call read_h5_slice_dp2_coord(file_id, 'coordinates', file_offset, nread, dcoords, ok)
        If (.not. ok) Call EmergencyStop('Cannot read /coordinates from HDF5 file', 5)

        Block
            Logical :: has_vel
            Call h5lexists_f(file_id, 'vx', has_vel, hdferr)
            If (has_vel) Then
                Call read_h5_slice_dp1(file_id, 'vx', file_offset, nread, dvx, ok)
                Call read_h5_slice_dp1(file_id, 'vy', file_offset, nread, dvy, ok)
                Call read_h5_slice_dp1(file_id, 'vz', file_offset, nread, dvz, ok)
            Else
                If (procID == 0) Print *, ' vx/vy/vz absent in snapshot - velocities set to zero in strct output.'
                dvx = 0.d0; dvy = 0.d0; dvz = 0.d0
            End If
        End Block

        Call read_h5_slice_dp1(file_id, 'mass',       file_offset, nread, dmass,  ok)
        Call read_h5_slice_dp1(file_id, 'birth_time', file_offset, nread, dbirth, ok)

        If (has_id_dataset .and. .not. invent_ids) Then
            Allocate(id4_from_file(nread))
            Call h5dopen_f(file_id, 'id', dset_id, hdferr)
            Call h5dget_space_f(dset_id, file_space, hdferr)
            Block
                Integer(HSIZE_T) :: offset(1), counts(1), mdims(1)
                Integer(HID_T)   :: mem_space
                offset(1) = file_offset
                counts(1) = nread
                Call h5sselect_hyperslab_f(file_space, H5S_SELECT_SET_F, offset, counts, hdferr)
                mdims(1) = nread
                Call h5screate_simple_f(1, mdims, mem_space, hdferr)
                Call h5dread_f(dset_id, H5T_NATIVE_INTEGER, id4_from_file, mdims, hdferr, mem_space, file_space)
                Call h5sclose_f(mem_space, hdferr)
            End Block
            Call h5sclose_f(file_space, hdferr)
            Call h5dclose_f(dset_id, hdferr)
        End If

        Call h5fclose_f(file_id, hdferr)
        Call h5close_f(hdferr)

        mynpart = int(nread, kind=4)

        ! ------------------------------------------------------------------
        ! 4. Auto-detect the DM particle mass (median across all ranks).
        !    DM-dominated assumption: median of all masses = DM mass.
        ! ------------------------------------------------------------------
        Allocate(mass_sample(mynpart))
        Do i = 1, mynpart
            mass_sample(i) = real(dmass(i), kind=SP)
        End Do
        If (mynpart > 0) Then
            local_median = median_inplace_sp(mass_sample, mynpart)
        Else
            local_median = 0.0_SP
        End If
        Deallocate(mass_sample)

        ! Collect each rank's local median on rank 0, take median of those.
        Allocate(per_rank_med(procNB))
        Call Mpi_Gather(local_median, 1, Mpi_Real, per_rank_med, 1, Mpi_Real, 0, Mpi_Comm_World, mpierr)
        If (procID == 0) Then
            global_median = median_inplace_sp(per_rank_med, procNB)
            Write(*, '(A, ES20.10)') ' Auto-detected DM particle mass (median): ', global_median
            Write(Ulog,*) 'Auto-detected DM particle mass:', global_median
        End If
        Deallocate(per_rank_med)
        Call Mpi_Bcast(global_median, 1, Mpi_Real, 0, Mpi_Comm_World, mpierr)
        massdmpart = global_median

        ! ------------------------------------------------------------------
        ! 5. Build per-particle SP arrays in normalized coords [0,1].
        !    Set tp = 0 for DM (matches the (tp==0 .and. mp==massdmpart) DM
        !    filter the downstream code uses), tp = |birth_time|+1 otherwise.
        ! ------------------------------------------------------------------
        Allocate(tmpx(3, mynpart), tmpv(3, mynpart))
        Allocate(tmpi(mynpart), tmpi4(mynpart))
        Allocate(tmpmp(mynpart), tmptp(mynpart))

        Do i = 1, mynpart
            tmpx(1, i) = real((dcoords(1, i) - phys_xmin) / Lbox, kind=SP)
            tmpx(2, i) = real((dcoords(2, i) - phys_xmin) / Lbox, kind=SP)
            tmpx(3, i) = real((dcoords(3, i) - phys_xmin) / Lbox, kind=SP)
            tmpv(1, i) = real(dvx(i), kind=SP)
            tmpv(2, i) = real(dvy(i), kind=SP)
            tmpv(3, i) = real(dvz(i), kind=SP)
            tmpmp(i)   = real(dmass(i), kind=SP)
            If (tmpmp(i) == massdmpart) Then
                tmptp(i) = 0.0_SP
            Else
                tmptp(i) = real(abs(dbirth(i)), kind=SP) + 1.0_SP
            End If
        End Do
        Deallocate(dcoords, dvx, dvy, dvz, dmass, dbirth)

        ! Clamp positions to [0,1) — periodic wrap of x==1.0.
        Do i = 1, mynpart
            Do idim = 1, 3
                If (tmpx(idim, i) >= 1.0_SP) tmpx(idim, i) = tmpx(idim, i) - 1.0_SP
                If (tmpx(idim, i) <  0.0_SP) tmpx(idim, i) = tmpx(idim, i) + 1.0_SP
            End Do
        End Do

        ! ------------------------------------------------------------------
        ! 6. Invent globally-unique sequential IDs. id4 = original /id if available,
        !    otherwise mirror the invented id (downstream code assumes id4 exists).
        ! ------------------------------------------------------------------
        Allocate(tabtmp(0:procNB-1), tabnpart(0:procNB-1))
        tabtmp = 0
        tabnpart = 0
        tabtmp(procID) = int(mynpart, kind=PRI)
        Call Mpi_AllReduce(tabtmp, tabnpart, procNB, MPI_PRI, Mpi_Sum, Mpi_Comm_World, mpierr)

        idmin = 0
        Do i = 0, procID - 1
            idmin = idmin + tabnpart(i)
        End Do

        Do i = 1, mynpart
            tmpi(i) = idmin + int(i, kind=PRI)
        End Do

        If (has_id_dataset .and. .not. invent_ids) Then
            tmpi4(1:mynpart) = id4_from_file(1:mynpart)
            Deallocate(id4_from_file)
        Else
            ! Mirror invented IDs (truncate to PRI4 — fine for <2^31 particles).
            Do i = 1, mynpart
                tmpi4(i) = int(tmpi(i), kind=PRI4)
            End Do
        End If

        If (procID == 0) Then
            Write(*, '(A, I0)') ' Sequential IDs created. Local range 1..N: ', mynpart
            Write(Ulog,*) 'Created invented IDs across ranks; total = ', sum(tabnpart)
        End If
        Deallocate(tabtmp, tabnpart)

        ! Timing
        If (dotimings) Then
            Call Mpi_Barrier(MPICube, mpierr)
            timeInt = Mpi_Wtime()
            tInitRead = timeInt - time0
            tReadfile = 0.d0
        End If

        ! ------------------------------------------------------------------
        ! 7. Geographic distribution across the MPI cartesian grid (same logic
        !    as in the RAMSES path: each rank ends up owning particles in its
        !    own [xmin..xmax) x [ymin..ymax) x [zmin..zmax) sub-cube).
        ! ------------------------------------------------------------------
        nsd = int(procNB ** (1.0/3.0))
        deltasd = 1.0_SP / real(nsd, kind=SP)

        If (procID == 0) Then
            Write(*, *) 'Number of subdomains per dim: ', nsd
            Write(*, *) 'Size of each subdomain      : ', deltasd
        End If

        Allocate(npartv(procNB), npartvloc(procNB))
        npartv = 0
        npartvloc = 0

        Do i = 1, mynpart
            Do idim = 1, 3
                If (tmpx(idim, i) == 1.0_SP) tmpx(idim, i) = 0.0_SP
            End Do
            n_i = int(tmpx(1, i) / deltasd)
            n_j = int(tmpx(2, i) / deltasd)
            n_k = int(tmpx(3, i) / deltasd)
            If (n_i >= nsd) n_i = nsd - 1
            If (n_j >= nsd) n_j = nsd - 1
            If (n_k >= nsd) n_k = nsd - 1
            ind = nsd*nsd * n_i + nsd * n_j + n_k + 1
            npartvloc(ind) = npartvloc(ind) + 1
        End Do

        Call Mpi_AllReduce(npartvloc, npartv, procNB, Mpi_Integer, Mpi_Sum, Mpi_Comm_World, mpierr)

        If (dotimings) Then
            tReadFile = Mpi_Wtime() - timeInt
            timeInt = Mpi_Wtime()
        End If

        Allocate(x(3, npartv(procID+1)))
        Allocate(v(3, npartv(procID+1)))
        Allocate(id(npartv(procID+1)))
        Allocate(id4(npartv(procID+1)))
        Allocate(mp(npartv(procID+1)))
        Allocate(tp(npartv(procID+1)))
        ! metal is always .false. for dyablo, but allocate zp if requested for compat
        If (metal) Allocate(zp(npartv(procID+1)))

        recvpoint = 1

        processus: Do i = 1, procNB - 1
            dest = mod(procID + i, procNB)
            prov = mod(procID + procNB - i, procNB)
            Call Mpi_Cart_coords(MPICube, dest, 3, destCoord, mpierr)

            Call Mpi_Isend(npartvloc(dest+1), 1, Mpi_Integer, dest, procID, MPICube, mpireqs1, mpierr)
            Call Mpi_Irecv(nrecv,             1, Mpi_Integer, prov, prov,   MPICube, mpireqr1, mpierr)

            xmin =  destCoord(1)      * deltasd
            xmax = (destCoord(1) + 1) * deltasd
            ymin =  destCoord(2)      * deltasd
            ymax = (destCoord(2) + 1) * deltasd
            zmin =  destCoord(3)      * deltasd
            zmax = (destCoord(3) + 1) * deltasd

            If (npartvloc(dest+1) > 0) Then
                Allocate(tmpsendx(3, npartvloc(dest+1)))
                Allocate(tmpsendv(3, npartvloc(dest+1)))
                Allocate(tmpsendi(   npartvloc(dest+1)))
                Allocate(tmpsendi4(  npartvloc(dest+1)))
                Allocate(tmpsendmp(  npartvloc(dest+1)))
                Allocate(tmpsendtp(  npartvloc(dest+1)))

                ind = 1
                Do j = 1, mynpart
                    If (tmpx(1,j) >= xmin .and. tmpx(1,j) < xmax .and. &
                        tmpx(2,j) >= ymin .and. tmpx(2,j) < ymax .and. &
                        tmpx(3,j) >= zmin .and. tmpx(3,j) < zmax) Then
                        tmpsendx(:,ind) = tmpx(:,j)
                        tmpsendv(:,ind) = tmpv(:,j)
                        tmpsendi(ind)   = tmpi(j)
                        tmpsendi4(ind)  = tmpi4(j)
                        tmpsendmp(ind)  = tmpmp(j)
                        tmpsendtp(ind)  = tmptp(j)
                        ind = ind + 1
                    End If
                End Do
                If (ind /= npartvloc(dest+1) + 1) Then
                    Call EmergencyStop('Particle binning mismatch in dyablo_lecture (send count)', 2)
                End If
            End If

            Call Mpi_Wait(mpireqs1, mpistat, mpierr)
            Call Mpi_Wait(mpireqr1, mpistat, mpierr)

            If (npartvloc(dest+1) /= 0) Then
                Call Mpi_Isend(tmpsendx,  3*npartvloc(dest+1), Mpi_Real,    dest, procID,           MPICube, mpireqs1, mpierr)
                Call Mpi_Isend(tmpsendv,  3*npartvloc(dest+1), Mpi_Real,    dest, dest,             MPICube, mpireqs2, mpierr)
                Call Mpi_Isend(tmpsendi,    npartvloc(dest+1), MPI_PRI,     dest, procID,           MPICube, mpireqs3, mpierr)
                Call Mpi_Isend(tmpsendi4,   npartvloc(dest+1), MPI_PRI4,    dest, procID,           MPICube, mpireqs4, mpierr)
                Call Mpi_Isend(tmpsendmp,   npartvloc(dest+1), Mpi_Real,    dest, procID+1*procNB,  MPICube, mpireqs5, mpierr)
                Call Mpi_Isend(tmpsendtp,   npartvloc(dest+1), Mpi_Real,    dest, procID+2*procNB,  MPICube, mpireqs6, mpierr)
            End If
            If (nrecv /= 0) Then
                Call Mpi_Irecv(x(1,recvpoint), 3*nrecv, Mpi_Real, prov, prov,            MPICube, mpireqr1, mpierr)
                Call Mpi_Irecv(v(1,recvpoint), 3*nrecv, Mpi_Real, prov, procID,          MPICube, mpireqr2, mpierr)
                Call Mpi_Irecv(id(recvpoint),    nrecv, MPI_PRI,  prov, prov,            MPICube, mpireqr3, mpierr)
                Call Mpi_Irecv(id4(recvpoint),   nrecv, MPI_PRI4, prov, prov,            MPICube, mpireqr4, mpierr)
                Call Mpi_Irecv(mp(recvpoint),    nrecv, Mpi_Real, prov, prov+1*procNB,   MPICube, mpireqr5, mpierr)
                Call Mpi_Irecv(tp(recvpoint),    nrecv, Mpi_Real, prov, prov+2*procNB,   MPICube, mpireqr6, mpierr)
                recvpoint = recvpoint + nrecv
            End If

            If (npartvloc(dest+1) /= 0) Then
                Call Mpi_Wait(mpireqs1, mpistat, mpierr); Deallocate(tmpsendx)
                Call Mpi_Wait(mpireqs2, mpistat, mpierr); Deallocate(tmpsendv)
                Call Mpi_Wait(mpireqs3, mpistat, mpierr); Deallocate(tmpsendi)
                Call Mpi_Wait(mpireqs4, mpistat, mpierr); Deallocate(tmpsendi4)
                Call Mpi_Wait(mpireqs5, mpistat, mpierr); Deallocate(tmpsendmp)
                Call Mpi_Wait(mpireqs6, mpistat, mpierr); Deallocate(tmpsendtp)
            End If
            If (nrecv /= 0) Then
                Call Mpi_Wait(mpireqr1, mpistat, mpierr)
                Call Mpi_Wait(mpireqr2, mpistat, mpierr)
                Call Mpi_Wait(mpireqr3, mpistat, mpierr)
                Call Mpi_Wait(mpireqr4, mpistat, mpierr)
                Call Mpi_Wait(mpireqr5, mpistat, mpierr)
                Call Mpi_Wait(mpireqr6, mpistat, mpierr)
            End If
        End Do processus

        ! Copy this rank's own resident particles to the final arrays
        xmin =  CubeCoord(1)      * deltasd
        xmax = (CubeCoord(1) + 1) * deltasd
        ymin =  CubeCoord(2)      * deltasd
        ymax = (CubeCoord(2) + 1) * deltasd
        zmin =  CubeCoord(3)      * deltasd
        zmax = (CubeCoord(3) + 1) * deltasd

        ind = 0
        Do j = 1, mynpart
            If (tmpx(1,j) >= xmin .and. tmpx(1,j) < xmax .and. &
                tmpx(2,j) >= ymin .and. tmpx(2,j) < ymax .and. &
                tmpx(3,j) >= zmin .and. tmpx(3,j) < zmax) Then
                x(:, recvpoint + ind) = tmpx(:, j)
                v(:, recvpoint + ind) = tmpv(:, j)
                id(recvpoint + ind)   = tmpi(j)
                id4(recvpoint + ind)  = tmpi4(j)
                mp(recvpoint + ind)   = tmpmp(j)
                tp(recvpoint + ind)   = tmptp(j)
                ind = ind + 1
            End If
        End Do

        If (recvpoint + ind /= npartv(procID+1) + 1) Then
            Write(tmpstr1, '(I9.9)') recvpoint + ind
            Write(tmpstr2, '(I9.9)') npartv(procID+1) + 1
            Call EmergencyStop('Wrong particles count after send/recv: '//tmpstr1//' vs '//tmpstr2, 2)
        End If

        mynpart = npartv(procID+1)

        Deallocate(tmpx, tmpv, tmpi, tmpi4, tmpmp, tmptp)
        Deallocate(npartv, npartvloc)

        If (dotimings) Then
            tTailPart = Mpi_Wtime() - timeInt
            tReadRA   = Mpi_Wtime() - time0
        End If

        Call Mpi_Barrier(MPICube, mpierr)

        ! ------------------------------------------------------------------
        ! 8. Split DM and stars (same logic as in ramses_lecture).
        !    After this, x/v/id/id4/mp/tp contain only DM particles, while
        !    xstar/vstar/... hold the stars on this rank.
        ! ------------------------------------------------------------------
        multimass_star = .true.
        multimass_dm   = .true.

        If (star) Then
            nstar = 0
            ndm   = 0
            Do ipart = 1, mynpart
                If (tp(ipart) == 0.0_SP .and. mp(ipart) == massdmpart) Then
                    ndm = ndm + 1
                Else
                    nstar = nstar + 1
                End If
            End Do

            Write(*, *) procID, ' counted dm/star: ', ndm, nstar

            If (nstar > 0) Then
                ! split x
                Allocate(tmpx(3, mynpart)); tmpx = x; Deallocate(x)
                Allocate(x(3, ndm)); Allocate(xstar(3, nstar))
                istar = 0; idm = 0
                Do ipart = 1, mynpart
                    If (tp(ipart) == 0.0_SP .and. mp(ipart) == massdmpart) Then
                        idm = idm + 1
                        x(:, idm) = tmpx(:, ipart)
                    Else
                        istar = istar + 1
                        xstar(:, istar) = tmpx(:, ipart)
                    End If
                End Do
                Deallocate(tmpx)

                ! split v
                Allocate(tmpv(3, mynpart)); tmpv = v; Deallocate(v)
                Allocate(v(3, ndm)); Allocate(vstar(3, nstar))
                istar = 0; idm = 0
                Do ipart = 1, mynpart
                    If (tp(ipart) == 0.0_SP .and. mp(ipart) == massdmpart) Then
                        idm = idm + 1
                        v(:, idm) = tmpv(:, ipart)
                    Else
                        istar = istar + 1
                        vstar(:, istar) = tmpv(:, ipart)
                    End If
                End Do
                Deallocate(tmpv)

                ! split id, id4
                Allocate(tmpi(mynpart));  tmpi  = id;  Deallocate(id)
                Allocate(tmpi4(mynpart)); tmpi4 = id4; Deallocate(id4)
                Allocate(id(ndm)); Allocate(id4(ndm)); Allocate(idstar(nstar))
                istar = 0; idm = 0
                Do ipart = 1, mynpart
                    If (tp(ipart) == 0.0_SP .and. mp(ipart) == massdmpart) Then
                        idm = idm + 1
                        id(idm)  = tmpi(ipart)
                        id4(idm) = tmpi4(ipart)
                    Else
                        istar = istar + 1
                        idstar(istar) = tmpi(ipart)
                    End If
                End Do
                Deallocate(tmpi, tmpi4)

                ! split mp
                Allocate(tmpmp(mynpart)); tmpmp = mp; Deallocate(mp)
                Allocate(mp(ndm)); Allocate(mpstar(nstar))
                istar = 0; idm = 0
                Do ipart = 1, mynpart
                    If (tp(ipart) == 0.0_SP .and. tmpmp(ipart) == massdmpart) Then
                        idm = idm + 1
                        mp(idm) = tmpmp(ipart)
                    Else
                        istar = istar + 1
                        mpstar(istar) = tmpmp(ipart)
                    End If
                End Do
                Deallocate(tmpmp)

                ! split tp -> tpstar
                Allocate(tmptp(mynpart)); tmptp = tp; Deallocate(tp)
                Allocate(tp(ndm)); Allocate(tpstar(nstar))
                istar = 0; idm = 0
                Do ipart = 1, mynpart
                    If (tmptp(ipart) == 0.0_SP) Then
                        idm = idm + 1
                        tp(idm) = tmptp(ipart)
                    Else
                        istar = istar + 1
                        tpstar(istar) = tmptp(ipart)
                    End If
                End Do
                Deallocate(tmptp)

                mynpart = ndm
            End If
        Else
            nstar = 0
        End If

        Call Mpi_Barrier(MPICube, mpierr)
        tmplongint = mynpart
        Call Mpi_AllReduce(tmplongint, nptot, 1, MPI_PRI, Mpi_Sum, Mpi_Comm_World, mpierr)
        tmplongint = nstar
        Call Mpi_AllReduce(tmplongint, nstot, 1, MPI_PRI, Mpi_Sum, Mpi_Comm_World, mpierr)
        If (procID == 0) Then
            Write(*, *) 'There are ', nptot, ' DM particles'
            Write(*, *) 'There are ', nstot, ' star particles'
            Write(Ulog,*) 'There are ', nptot, ' DM particles'
            Write(Ulog,*) 'There are ', nstot, ' star particles'
        End If

        ! Re-sequence DM particle IDs to [1..nptot] after star filtering.
        ! IDs were invented for all particles (DM + stars) before filtering,
        ! so some DM particles may carry IDs > nptot, breaking procMasse's
        ! smin/smax scan ranges in modfofpara.f90.
        Allocate(tabtmp(0:procNB-1), tabnpart(0:procNB-1))
        tabtmp = 0
        tabnpart = 0
        tabtmp(procID) = int(mynpart, kind=PRI)
        Call Mpi_AllReduce(tabtmp, tabnpart, procNB, MPI_PRI, Mpi_Sum, Mpi_Comm_World, mpierr)
        idmin = 0
        Do i = 0, procID - 1
            idmin = idmin + tabnpart(i)
        End Do
        Do i = 1, mynpart
            id(i) = idmin + int(i, kind=PRI)
        End Do
        Deallocate(tabtmp, tabnpart)

        Call Mpi_Barrier(MPICube, mpierr)
        Write(*,*) procID, ' dm   mynpart  ', mynpart
        If (mynpart > 0) Then
            Write(*,*) procID, ' dm   minmax(x) ', minval(x), maxval(x)
            Write(*,*) procID, ' dm   minmax(v) ', minval(v), maxval(v)
            Write(*,*) procID, ' dm   minmax(id)', minval(id), maxval(id)
            If (multimass_dm) Write(*,*) procID, ' dm   minmax(m) ', minval(mp), maxval(mp)
        End If

    End Subroutine dyablo_lecture


    !=======================================================================
    ! Write the cube of particles to be analysed by current process
    Subroutine outputcube()

        Use modparam
        Use modvariable
        Use modmpicom
        Implicit none

        Integer(kind=PRI) :: i, j
        Character(len=1200) :: filecube
        Character(len=5)  :: pid_char

        Write(pid_char(1:5),'(I5.5)') procID
        filecube = out_path(trim(root)//"_cube_"//pid_char)
        Open(Unit=Ucub, file=filecube, Form='Unformatted')

        Write(Ucub) int(mynpart, kind=4)
        Write(Ucub) procID
        Write(Ucub) xmin, xmax, ymin, ymax, zmin, zmax
        Write(Ucub) ((x(j,i), j=1,3), i=1, mynpart)
        Write(Ucub) ((v(j,i), j=1,3), i=1, mynpart)
        Write(Ucub) (id(i), i=1, mynpart)
        If (multimass_dm) Write(Ucub) (mp(i), i=1, mynpart)
        Close(Ucub)

    End Subroutine outputcube

    !=======================================================================
    Subroutine outputcubestar()

        Use modparam
        Use modvariable
        Use modmpicom
        Implicit none

        Integer(kind=PRI) :: i, j
        Character(len=1200) :: filecube
        Character(len=5)  :: pid_char

        If (nstar > 0) Then
            Write(pid_char(1:5),'(I5.5)') procID
            filecube = out_path(trim(root)//"_cubestar_"//pid_char)
            Open(Unit=Ucub, file=filecube, Form='Unformatted')

            Write(Ucub) int(nstar, kind=4)
            Write(Ucub) procID
            Write(Ucub) xmin, xmax, ymin, ymax, zmin, zmax
            Write(Ucub) ((xstar(j,i), j=1,3), i=1, nstar)
            Write(Ucub) ((vstar(j,i), j=1,3), i=1, nstar)
            Write(Ucub) (idstar(i), i=1, nstar)
            If (multimass_star) Write(Ucub) (mpstar(i), i=1, nstar)
            If (star) Then
                Write(Ucub) (tpstar(i), i=1, nstar)
                If (metal) Write(Ucub) (zpstar(i), i=1, nstar)
            End If
            Close(Ucub)
        End If

    End Subroutine outputcubestar

    !=======================================================================
    Subroutine outputmass(ns, smin, nbs, massamas, cdmamas)

        Use modparam
        Use modvariable
        Use modmpicom
        Implicit none

        Integer(kind=4),                   intent(in) :: ns
        Integer(kind=PRI),                 intent(in) :: smin
        Integer(kind=4)                               :: nbs
        Integer(kind=4),  dimension(nbs),  intent(in) :: massamas
        Real(kind=DP),    dimension(3,nbs),intent(in) :: cdmamas

        Character(len=1200) :: fileamas
        Character(len=5)  :: pid_char
        Integer(kind=4) :: i

        Write(pid_char(1:5),'(I5.5)') procID
        fileamas = out_path(trim(root)//"_masst_"//pid_char)

        Open(Umas, File=trim(fileamas), Status='Unknown', Form='Unformatted')
        Write(Umas) int(ns, kind=4)
        Do i = 1, nbs
            If (massamas(i) >= Mmin) Write(Umas) int(i+smin-1, kind=8), massamas(i), real(cdmamas(:,i), kind=SP)
        End Do
        Close(Umas)

    End Subroutine outputmass

    !=======================================================================
    Subroutine outputstruct(np, ns, nbs, idf, xf, vf, massamas, idf4)

        Use modparam
        Use modvariable
        Use modmpicom
        Implicit none

        Integer(kind=4),                  intent(in) :: np, ns, nbs
        Integer(kind=PRI),dimension(np),  intent(in) :: idf
        Integer(kind=PRI4),dimension(np), intent(in) :: idf4
        Integer(kind=4),  dimension(nbs), intent(in) :: massamas
        Real   (kind=SP), dimension(3,np),intent(in) :: xf, vf

        Integer(kind=4) :: i, j, k, b
        Character(len=1200) :: filestrct
        Character(len=5)  :: pid_char

        Write(pid_char(1:5),'(I5.5)') procID
        filestrct = out_path(trim(root)//"_strct_"//pid_char)

        Open(Ustr, File=trim(filestrct), Status='Unknown', Form='Unformatted')
        Write(Ustr) int(ns, kind=4)

        b = 1
        Do i = 1, nbs
            If (massamas(i) >= Mmin) Then
                Write(Ustr) massamas(i)
                Write(Ustr) ((xf(k,j), k=1,3), j=b, b+massamas(i)-1)
                Write(Ustr) ((vf(k,j), k=1,3), j=b, b+massamas(i)-1)
                Write(Ustr) (idf4(j), j=b, b+massamas(i)-1)
                b = b + massamas(i)
            Else
                b = b + massamas(i)
            End If
        End Do

        Close(Ustr)

    End Subroutine outputstruct

End Module modio
