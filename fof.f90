!------------------------------------------------------------------------
! Parallel Friend of Friend halo finder, adapted for dyablo HDF5 snapshots.
!
! Original sequential FoF: Edouard Audit (CEA).
! Parallel MPI implementation: Fabrice Roy (CNRS / LUTH, Observatoire de Paris).
! Dyablo IO + auto-DM-mass + ID handling: this branch.
!------------------------------------------------------------------------


Program friend

  Use modconst
  Use modvariable
  Use modparam
  Use modmpicom
  Use modfofpara
  Use modio
  Use modtiming

  Implicit none

  Integer :: errcode, ios
  Character(len=1200):: filelog

  Call Mpi_Init(mpierr)
  Call Mpi_Comm_size(Mpi_Comm_World, procNB, mpierr)

  dims = int(procNB**(1./3.))
  periods = .true.

  Call Mpi_Cart_create(Mpi_Comm_World, 3, dims, periods, .true., MPICube, mpierr)

  Call Mpi_Comm_rank  (MPICube, procID,    mpierr)
  Call Mpi_Cart_coords(MPICube, procID, 3, CubeCoord, mpierr)

  Call Mpi_Cart_shift(MPICube, 0, 1, voisin(1), voisin(2), mpierr)
  Call Mpi_Cart_shift(MPICube, 1, 1, voisin(3), voisin(4), mpierr)
  Call Mpi_Cart_shift(MPICube, 2, 1, voisin(5), voisin(6), mpierr)

  ! Packed-broadcast buffer:
  !   512*5  : 5 x Character(len=512) root, outdir, pathinput, namepart, nameinfo
  !   3      : 1 x Character(len=3)   code_index
  !   3*ks   : 3 x Integer(kind=4)    grpsize, Mmin, Mmax
  !   kp     : 1 x Real(kind=4)       perco
  !   7*4    : 7 x Logical            star, metal, outcube, dofof, readfromcube, dotimings, invent_ids
  h_length = 512*5 + 3 + 3 * bit_size(Mmin)/8 + kind(perco) + 7*4
  allocate (header(0:h_length-1))
  h_pos = 0

  If (procID == 0) Then

     Open(10, file='fof.in', iostat=errcode)
     If (errcode > 0) Then
        Print *, '** Error opening input file fof.in. Please check this file. **'
        Call Mpi_Abort(Mpi_Comm_World, errcode, mpierr)
     End If

     Read(10, '(A512)') root
     Read(10, *)       code_index
     Read(10, '(A512)') pathinput
     Read(10, '(A512)') namepart
     Read(10, '(A512)') nameinfo
     Read(10, *)       grpsize
     Read(10, *)       perco
     Read(10, *)       Mmin
     Read(10, *)       Mmax
     Read(10, *)       star
     Read(10, *)       metal
     Read(10, *)       outcube
     Read(10, *)       dofof
     Read(10, *)       readfromcube
     Read(10, *)       dotimings
     Read(10, *)       invent_ids
     Read(10, '(A512)', iostat=ios) outdir
     If (ios /= 0) outdir = '.'
     Close(10)

     ! Create output dir if needed (mkdir -p is idempotent). procID==0 only.
     If (len_trim(outdir) > 0 .and. trim(outdir) /= '.') Then
        Call execute_command_line('mkdir -p '//trim(outdir), wait=.true.)
     End If

     Print *, 'Parallel FoF (dyablo branch)'
     Print *, procNB, ' processes:'

     filelog = out_path(trim(root)//'.log')
     Open(Unit=Ulog, file=filelog)
     Write(Ulog, *) 'Parallel FoF (dyablo branch)'
     Write(Ulog, *) procNB, ' processes:'

     Call outputparameters()
     Write(*, *) "Output parameters written"

     Call Mpi_Pack(        root,512, Mpi_Character, header, h_length, h_pos, Mpi_Comm_World, mpierr)
     Call Mpi_Pack(      outdir,512, Mpi_Character, header, h_length, h_pos, Mpi_Comm_World, mpierr)
     Call Mpi_Pack(  code_index,  3, Mpi_Character, header, h_length, h_pos, Mpi_Comm_World, mpierr)
     Call Mpi_Pack(   pathinput,512, Mpi_Character, header, h_length, h_pos, Mpi_Comm_World, mpierr)
     Call Mpi_Pack(    namepart,512, Mpi_Character, header, h_length, h_pos, Mpi_Comm_World, mpierr)
     Call Mpi_Pack(    nameinfo,512, Mpi_Character, header, h_length, h_pos, Mpi_Comm_World, mpierr)
     Call Mpi_Pack(     grpsize,  1,   Mpi_Integer, header, h_length, h_pos, Mpi_Comm_World, mpierr)
     Call Mpi_Pack(       perco,  1,      Mpi_Real, header, h_length, h_pos, Mpi_Comm_World, mpierr)
     Call Mpi_Pack(        Mmin,  1,   Mpi_Integer, header, h_length, h_pos, Mpi_Comm_World, mpierr)
     Call Mpi_Pack(        Mmax,  1,   Mpi_Integer, header, h_length, h_pos, Mpi_Comm_World, mpierr)
     Call Mpi_Pack(        star,  1,   Mpi_Logical, header, h_length, h_pos, Mpi_Comm_World, mpierr)
     Call Mpi_Pack(       metal,  1,   Mpi_Logical, header, h_length, h_pos, Mpi_Comm_World, mpierr)
     Call Mpi_Pack(     outcube,  1,   Mpi_Logical, header, h_length, h_pos, Mpi_Comm_World, mpierr)
     Call Mpi_Pack(       dofof,  1,   Mpi_Logical, header, h_length, h_pos, Mpi_Comm_World, mpierr)
     Call Mpi_Pack(readfromcube,  1,   Mpi_Logical, header, h_length, h_pos, Mpi_Comm_World, mpierr)
     Call Mpi_Pack(   dotimings,  1,   Mpi_Logical, header, h_length, h_pos, Mpi_Comm_World, mpierr)
     Call Mpi_Pack(  invent_ids,  1,   Mpi_Logical, header, h_length, h_pos, Mpi_Comm_World, mpierr)

     Print *, ''
     Print *, "Root:                                           ", trim(root)
     Print *, "Output directory:                               ", trim(outdir)
     Print *, "Type of input files:                            ", code_index
     Print *, "Path to input files:                            ", trim(pathinput)
     Print *, "Particle file name (HDF5 for DY):               ", trim(namepart)
     Print *, "Info / .ini file name:                          ", trim(nameinfo)
     Print *, "Size of groups of inputs (RAMSES legacy):       ", grpsize
     Print *, "Percolation parameter:                          ", perco
     Print *, "Minimum mass of halo to be analyzed:            ", Mmin
     Print *, "Maximum mass of halo to be analyzed:            ", Mmax
     Print *, "Star particles present in the snapshot?         ", star
     Print *, "Metallicities present in the snapshot?          ", metal
     Print *, "Write cubes of particles:                       ", outcube
     Print *, "Perform friends of friends halo detection:      ", dofof
     Print *, "Read particles from cube files:                 ", readfromcube
     Print *, "Perform timings (imply extra synchronisations): ", dotimings
     Print *, "Force inventing new IDs (ignore snapshot /id):  ", invent_ids
     Print *, ''

     Write(Ulog,*) 'INPUT PARAMETERS fof.in'
     Write(Ulog,*) root
     Write(Ulog,*) outdir
     Write(Ulog,*) code_index
     Write(Ulog,*) pathinput
     Write(Ulog,*) namepart
     Write(Ulog,*) nameinfo
     Write(Ulog,*) grpsize
     Write(Ulog,*) perco
     Write(Ulog,*) Mmin
     Write(Ulog,*) Mmax
     Write(Ulog,*) star
     Write(Ulog,*) metal
     Write(Ulog,*) outcube
     Write(Ulog,*) dofof
     Write(Ulog,*) readfromcube
     Write(Ulog,*) dotimings
     Write(Ulog,*) invent_ids

  End If

  Call Mpi_Bcast(header, h_length, Mpi_Packed, 0, Mpi_Comm_World, mpierr)

  If (procID /= 0) Then
     Call Mpi_Unpack(header, h_length, h_pos,         root,512, Mpi_Character, Mpi_Comm_World, mpierr)
     Call Mpi_Unpack(header, h_length, h_pos,       outdir,512, Mpi_Character, Mpi_Comm_World, mpierr)
     Call Mpi_Unpack(header, h_length, h_pos,   code_index,  3, Mpi_Character, Mpi_Comm_World, mpierr)
     Call Mpi_Unpack(header, h_length, h_pos,    pathinput,512, Mpi_Character, Mpi_Comm_World, mpierr)
     Call Mpi_Unpack(header, h_length, h_pos,     namepart,512, Mpi_Character, Mpi_Comm_World, mpierr)
     Call Mpi_Unpack(header, h_length, h_pos,     nameinfo,512, Mpi_Character, Mpi_Comm_World, mpierr)
     Call Mpi_Unpack(header, h_length, h_pos,      grpsize,  1,   Mpi_Integer, Mpi_Comm_World, mpierr)
     Call Mpi_Unpack(header, h_length, h_pos,        perco,  1,      Mpi_Real, Mpi_Comm_World, mpierr)
     Call Mpi_Unpack(header, h_length, h_pos,         Mmin,  1,   Mpi_Integer, Mpi_Comm_World, mpierr)
     Call Mpi_Unpack(header, h_length, h_pos,         Mmax,  1,   Mpi_Integer, Mpi_Comm_World, mpierr)
     Call Mpi_Unpack(header, h_length, h_pos,         star,  1,   Mpi_Logical, Mpi_Comm_World, mpierr)
     Call Mpi_Unpack(header, h_length, h_pos,        metal,  1,   Mpi_Logical, Mpi_Comm_World, mpierr)
     Call Mpi_Unpack(header, h_length, h_pos,      outcube,  1,   Mpi_Logical, Mpi_Comm_World, mpierr)
     Call Mpi_Unpack(header, h_length, h_pos,        dofof,  1,   Mpi_Logical, Mpi_Comm_World, mpierr)
     Call Mpi_Unpack(header, h_length, h_pos, readfromcube,  1,   Mpi_Logical, Mpi_Comm_World, mpierr)
     Call Mpi_Unpack(header, h_length, h_pos,    dotimings,  1,   Mpi_Logical, Mpi_Comm_World, mpierr)
     Call Mpi_Unpack(header, h_length, h_pos,   invent_ids,  1,   Mpi_Logical, Mpi_Comm_World, mpierr)
  End If

  Deallocate(header)

  multimass_star = .true.
  multimass_dm   = .true.
  If ((star .and. multimass_star) .or. multimass_dm) Then
     multimass = .true.
  Else
     multimass = .false.
  End If

  ! 'DY ' = dyablo HDF5 snapshot. RAMSES support has been removed in this branch.
  If (code_index == 'DY ' .or. code_index == 'DY3' .or. code_index == 'DYA') Then
     Call dyablo_lecture()
  Else
     If (procID == 0) Print *, '** Wrong file type. This branch only supports DY (dyablo). **'
     errcode = 2
     Call Mpi_Abort(Mpi_Comm_World, errcode, mpierr)
     Stop
  End If

  If (procID == 0) Then
     Write(* ,*)   'nz = ', nres, 'ngrid = ', ngrid
     Write(Ulog,*) 'nz = ', nres, 'ngrid = ', ngrid
  End If

  If (outcube) Then
     If (procID == 0) Print *,      'Write particles distributed in a cartesian grid'
     If (procID == 0) Write(Ulog,*) 'Write particles distributed in a cartesian grid'
     Call outputcube()
     If (multimass_dm) deallocate(mp)
     If (star) Then
        If (nstar > 0) Then
           Call outputcubestar()
           Deallocate(xstar, vstar, idstar)
           Deallocate(tpstar)
           If (metal) deallocate(zpstar)
           If (multimass_star) deallocate(mpstar)
        End If
     End If
  End If

  If (procID == 0) Then
     Print *, ' '
     Print *, 'Friends of Friends halo detection'
     Print *, ' '
     Write(Ulog, *) 'Friends of Friends halo detection'
  End If

  If (dofof) Then
     Call parafof()
  Else
     tFoF     = 0.0
     tFoFinit = 0.0
     tFoFloc  = 0.0
     tRaccord = 0.0
     tObs     = 0.0
     tOut     = 0.0
  End If

  If (dotimings .and. procID == 0) Then
     Print *, 'Friend of Friend termine'
     Print *, ''
     Print *, 'Temps d''execution:'
     Print *, 'Lecture:', tReadRA, ' dont'
     Print *, '        initialisation        :', tInitRead
     Print *, '        lecture des fichiers  :', tReadFile
     Print *, '        partage des particules:', tTailPart
     Print *, ''
     Print *, 'Friend of Friend:', tFoF, 'dont'
     Print *, '        initialisation:', tFoFinit
     Print *, '        FoF local     :', tFoFloc
     Print *, '        raccordement  :', tRaccord
     Print *, '        calcul d''observables:', tObs
     Print *, '        sorties:', tOut

     Write(Ulog,*) 'Friend of Friend termine'
     Write(Ulog,*) ''
     Write(Ulog,*) 'Temps d''execution:'
     Write(Ulog,*) 'Lecture:', tReadRA, ' dont'
     Write(Ulog,*) '        initialisation        :', tInitRead
     Write(Ulog,*) '        lecture des fichiers  :', tReadFile
     Write(Ulog,*) '        partage des particules:', tTailPart
     Write(Ulog,*) ''
     Write(Ulog,*) 'Friend of Friend:', tFoF, 'dont'
     Write(Ulog,*) '        initialisation:', tFoFinit
     Write(Ulog,*) '        FoF local     :', tFoFloc
     Write(Ulog,*) '        raccordement  :', tRaccord
     Write(Ulog,*) '        calcul d''observables:', tObs
     Write(Ulog,*) '        sorties:', tOut
  End If

  Deallocate(x, v, id)
  If (allocated(id4)) Deallocate(id4)
  If (allocated(mp))  Deallocate(mp)
  If (allocated(tp))  Deallocate(tp)
  Close(Ulog)

  Call Mpi_Finalize(mpierr)

End Program friend
