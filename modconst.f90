Module modconst

  Use mpi

#ifdef LONGREAL
  Integer, parameter :: PR=8
#else
  Integer, parameter :: PR=4
#endif

! Define LONGINT in makefile to analyze simulations with more than 2**31 particles.
#ifdef LONGINT
  Integer, parameter :: PRI = 8
  Integer, parameter :: MPI_PRI = Mpi_Integer8
  Integer, parameter :: PRI4 = 4
  Integer, parameter :: MPI_PRI4 = Mpi_Integer
#else
  Integer, parameter :: PRI = 4
  Integer, parameter :: MPI_PRI = Mpi_Integer
  Integer, parameter :: PRI4 = 4
  Integer, parameter :: MPI_PRI4 = Mpi_Integer
#endif


  Integer, parameter :: DP = kind(1.d0)
  Integer, parameter :: SP = kind(1.e0)

  ! Output Units
  Integer, parameter :: Ulog=50, Ucub=51, Umas=52, Ustr=53, Uopa=54

  ! DM particle mass (single precision, in input file units).
  ! Auto-detected from the median mass at read time. Stays SP so that
  ! the original mass==massdmpart comparisons against tmpmp(SP) still work.
  Real(kind=SP) :: massdmpart = 0.0_SP

  Integer, parameter :: verb=0

End Module modconst
