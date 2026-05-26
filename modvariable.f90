! arrays of positions/velocities/ids
Module modvariable

  Use modconst
  Real(   kind=SP ), dimension(:,:), allocatable :: x   ! particles positions - Always simple precision
  Real(   kind=SP ), dimension(:,:), allocatable :: v   ! particles velocities - Always simple precision
  Integer(kind=PRI), dimension(:),   allocatable :: id  ! invented globally-unique sequential IDs (used by FoF)
  Integer(kind=PRI4), dimension(:),  allocatable :: id4 ! original IDs from the snapshot (if /id exists), else = id
  Real(   kind=SP ), dimension(:), allocatable :: tp  ! particle formation time (if star)
  Real(   kind=SP ), dimension(:), allocatable :: zp  ! particle metallicity (if star) - unused for dyablo
  Real(   kind=SP ), dimension(:), allocatable :: mp  ! particle mass (if multimass)

  Real(   kind=SP ), dimension(:,:), allocatable :: xstar  !same as above but for star non DM particles
  Real(   kind=SP ), dimension(:,:), allocatable :: vstar
  Integer(kind=PRI), dimension(:),   allocatable :: idstar
  Real(   kind=SP ), dimension(:), allocatable :: tpstar
  Real(   kind=SP ), dimension(:), allocatable :: zpstar
  Real(   kind=SP ), dimension(:), allocatable :: mpstar

End Module modvariable
