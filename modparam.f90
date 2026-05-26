! Input parameters
Module modparam

  Use modconst
  Character(len=512) :: root               ! base name for output
  Character(len=512) :: outdir             ! directory for output files (created if missing; '.' or blank = current dir)
  Character(len=512) :: pathinput          ! path to the dyablo snapshot
  Character(len=512) :: nameinfo           ! dyablo .ini config file name (under pathinput)
  Character(len=512) :: namepart           ! dyablo particles HDF5 file name (under pathinput)
  Character(len=3)  :: code_index         ! input file type: 'DY ' = dyablo (only one supported here)
  Integer(kind=8)   :: Mmin, Mmax         ! min and max structure mass
  Integer(kind=4)   :: grpsize            ! kept for compatibility with the RAMSES tool; set to 0 for dyablo
  Real(kind=SP)     :: perco              ! percolation parameter for Friends of Friends halo detection
  Logical           :: outcube            ! should there be an output after the particles reading/sharing?
  Logical           :: star               ! are star particles mixed in with DM in the snapshot?
  Logical           :: metal              ! is metallicity stored? (always .false. for current dyablo)
  Logical           :: dofof              ! should the structures be detected?
  Logical           :: readfromcube       ! should the particles be read from cube files?
  Logical           :: dotimings          ! should there be timings?
  Logical           :: multimass          !do we need to read star or dm mass?
  Logical           :: multimass_dm, multimass_star ! do we want to write DM mass? Star mass?
  Logical           :: invent_ids         ! force inventing new IDs even when the snapshot has an /id dataset

End Module modparam
