MODULE ed_ham_mod

    use Operator_mod, only: Operator
    use runtime_error_mod, only: ERROR_GENERIC, Terminate_on_error
    use ed_state_mod, only: ed_state, N_fermions
    use iso_fortran_env, only: output_unit, error_unit

    IMPLICIT NONE
    PRIVATE
    PUBLIC :: ed_ham

    integer, parameter :: dp=kind(0.d0)  ! double precision

    TYPE ed_ham_part
!--------------------------------------------------------------------
!> @author
!> ALF-project
!>
!> @brief
!> A block of the Hamiltonian with constant particle number, for exact diagonalisation.
      private
      INTEGER :: N_particles, N_states
      complex(dp), allocatable :: H(:,:)
      real(dp)   , allocatable :: eigenval(:)
    END TYPE ed_ham_part


    TYPE ed_ham
!--------------------------------------------------------------------
!> @author
!> ALF-project
!>
!> @brief
!> Defines a Hamiltonian for exact diagonalisation.
        private
        logical :: restrict_N_part
        INTEGER :: N_orbitals, N_SUN, N_FL, N_states, N_part
        type(ed_ham_part), allocatable :: H_part(:)
        real(dp), allocatable :: eigenval(:)
        integer, allocatable :: list(:,:) ! list(i,1) = N_particles
                                          ! list(i,2) = index in subhamiltonian
      CONTAINS
        private
        PROCEDURE :: init => ed_ham_init
        PROCEDURE :: add_matrixelement => ed_ham_add_matrixelement
        PROCEDURE :: add_op_t => ed_ham_add_op_t
        PROCEDURE :: add_op_v => ed_ham_add_op_v

        PROCEDURE, public :: build_ham => ed_ham_build_ham
        PROCEDURE, public :: energy => ed_ham_energy
        PROCEDURE, public :: get_eigenvalues => ed_ham_get_eigenvalues
    END TYPE ed_ham

CONTAINS

    subroutine ed_ham_build_ham(this, ndim, N_SUN, OP_T, OP_V, dtau, N_part)
!--------------------------------------------------------------------
!> @author
!> ALF-project
!>
!> @brief
!> Builds Hamiltonian. Inculdes initialisation and calculation of eigenvalues.
      IMPLICIT NONE
      class(ed_ham), intent(inout) :: this
      integer, intent(in) :: ndim, N_SUN
      Type (Operator), intent(in) :: Op_V(:,:)
      Type (Operator), intent(in) :: Op_T(:,:)
      real(dp), intent(in) :: dtau
      integer , intent(in) :: N_part

      INTEGER :: n, N_FL

      print*, "Building ED-Hamiltonian"
      this%N_part = N_part
      this%restrict_N_part = (N_part > -1)
      N_FL = Size(OP_T,2)
      call this%init(ndim, N_SUN, N_FL)
      do n=1, Size(OP_T,1)
        call this%add_op_t( OP_T(n,:), dtau )
      enddo
      do n=1, Size(OP_V,1)
        call this%add_op_v( OP_V(n,:), dtau )
      enddo
      print*, "done Building ED-Hamiltonian"
      call ed_ham_test_hermitian(this)
      call ed_ham_eigenvalues(this)
    end subroutine ed_ham_build_ham


    function ed_ham_energy(this, beta)
!--------------------------------------------------------------------
!> @author
!> ALF-project
!>
!> @brief
!> Calculates system energy at given reciprocal temperature beta
      IMPLICIT NONE
      class(ed_ham), intent(inout) :: this
      real(dp)     , intent(in)    :: beta

      INTEGER  :: i
      real(dp) :: Z, E, ed_ham_energy

      Z = 0.d0
      E = 0.d0

      do i=1, size(this%eigenval)
        Z = Z + exp(-beta * this%eigenval(i))
        E = E + exp(-beta * this%eigenval(i)) * this%eigenval(i)
      enddo

      ed_ham_energy = E / Z
    end function ed_ham_energy


    subroutine ed_ham_get_eigenvalues(this, eigenvalues)
!--------------------------------------------------------------------
!> @author
!> ALF-project
!>
!> @brief
!> Returns eigenvalues of the Hamiltonian.
      IMPLICIT NONE
      class(ed_ham), intent(inout) :: this

      real(dp), allocatable, intent(out) :: eigenvalues(:)

      allocate( eigenvalues(size(this%eigenval)) )

      eigenvalues = this%eigenval

    end subroutine ed_ham_get_eigenvalues


    subroutine ed_ham_init(this, N_orbitals, N_SUN, N_FL)
!--------------------------------------------------------------------
!> @author
!> ALF-project
!>
!> @brief
!> Initialises Hamiltonian
      IMPLICIT NONE
      class(ed_ham), INTENT(INOUT) :: this
      integer, intent(in) :: N_orbitals, N_SUN, N_FL

      integer              :: i, N_particles
      integer, allocatable :: sizes(:)

      this%N_orbitals = N_orbitals
      this%N_SUN      = N_SUN
      this%N_FL       = N_FL
      this%N_states   = 2**(N_orbitals*N_SUN*N_FL)

      allocate(this%H_part(0:N_orbitals*N_SUN*N_FL), this%list(0:this%N_states-1,2))

      do i=0, N_orbitals*N_SUN*N_FL
        this%H_part(i)%N_particles = i
        this%H_part(i)%N_states    = 0
      enddo

      do i=0, this%N_states-1
        N_particles = N_fermions(i)
        this%H_part(N_particles)%N_states = this%H_part(N_particles)%N_states + 1
        this%list(i,1) = N_particles
        this%list(i,2) = this%H_part(N_particles)%N_states
      enddo

      do i=0, N_orbitals*N_SUN*N_FL
        if( .not.this%restrict_N_part .or. i == this%N_part ) then
          allocate( this%H_part(i)%H(this%H_part(i)%N_states, this%H_part(i)%N_states) )
          this%H_part(i)%H(:,:) = cmplx(0.d0, 0.d0, dp)
        endif
      enddo
    end subroutine ed_ham_init


    subroutine ed_ham_add_matrixelement(this, i, j, val)
!--------------------------------------------------------------------
!> @author
!> ALF-project
!>
!> @brief
!> Adds value val to Hamiltonian matrix element with indices i and j
      IMPLICIT NONE
      class(ed_ham), intent(inout) :: this
      integer      , intent(in)    :: i, j
      complex(dp)  , intent(in)    :: val

      if( abs(val) < 1e-11 ) return

      if ( this%list(i,1) .ne. this%list(j,1) ) then
        write(error_unit,*) "Error in ed_ham_add_matrixelement", this%list(i,1), this%list(j,1), val
        CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
      endif

      this%H_part(this%list(i,1))%H(this%list(i,2), this%list(j,2)) = &
       & this%H_part(this%list(i,1))%H(this%list(i,2), this%list(j,2)) + val
    end subroutine ed_ham_add_matrixelement


    subroutine ed_ham_add_op_t(this, op_t, dtau)
!--------------------------------------------------------------------
!> @author
!> ALF-project
!>
!> @brief
!> Adds single-particle opertator op_t to Hamiltonian
      IMPLICIT NONE
      class(ed_ham) , intent(inout) :: this
      type(Operator), intent(in)    :: Op_T(:)
      real(dp)      , intent(in)    :: dtau

      integer :: i, n1, n2, sigma, s
      type(ed_state) :: state

      call state%init(this%N_orbitals, this%N_SUN, this%N_FL)

      do i=0, this%N_states-1
        if( .not.this%restrict_N_part .or. N_fermions(i) == this%N_part ) then
          do s=0, this%N_FL-1
            do sigma=0, this%N_SUN-1
              do n1=1, op_t(s+1)%N
                do n2=1, op_t(s+1)%N
                  call state%set(i)
                  call state%annihil_e( op_t(s+1)%p(n1)-1, sigma, s )
                  call state%create_e ( op_t(s+1)%p(n2)-1, sigma, s )

                  call this%add_matrixelement(state%get_i(), i, &
                      & state%get_factor() * op_t(s+1)%O(n2,n1) * op_t(s+1)%g / (-dtau) )
                enddo
              enddo
            enddo
          enddo
        endif
      enddo
    end subroutine ed_ham_add_op_t


    subroutine ed_ham_add_op_v(this, op_v, dtau)
!--------------------------------------------------------------------
!> @author
!> ALF-project
!>
!> @brief
!> Adds type 2 operator op_v to Hamiltonian
      IMPLICIT NONE
      class(ed_ham) , intent(inout) :: this
      type(Operator), intent(in)    :: op_v(:)
      real(dp)      , intent(in)    :: dtau

      integer :: i, n1, n2, m1, m2, s, s2, sigma, sigma2
      complex (dp) :: temp
      type(ed_state) :: state

      if( op_v(1)%type .ne. 2 ) then
          write(error_unit,*) "ED only implemented for OP_V type 2"
          CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
      endif

      call state%init(this%N_orbitals, this%N_SUN, this%N_FL)

      temp = cmplx(0.d0, 0.d0, dp)
      do s=1, this%N_FL
          temp = temp + op_v(s)%g * op_v(s)%alpha
      enddo
      temp = temp**2 * cmplx( real(this%N_SUN**2,dp) / (-dtau), 0.d0, dp )

      do i=0, this%N_states-1
        if( .not.this%restrict_N_part .or. N_fermions(i) == this%N_part ) then
          call this%add_matrixelement(i, i, temp )
          do s=0, this%N_FL-1
            do s2=0, this%N_FL-1
              do sigma=0, this%N_SUN-1
                do n1=1, op_v(s+1)%N
                  do n2=1, op_v(s+1)%N
                    call state%set(i)
                    call state%annihil_e( op_v(s+1)%p(n1)-1, sigma, s )
                    call state%create_e ( op_v(s+1)%p(n2)-1, sigma, s )

                    call this%add_matrixelement(state%get_i(), i, state%get_factor() * &
                    & op_v(s+1)%O(n2,n1) * op_v(s+1)%g * op_v(s2+1)%g * op_v(s2+1)%alpha * &
                    & cmplx( real(2*this%N_SUN,dp) / (-dtau), 0.d0, dp ) )

                    do sigma2=0, this%N_SUN-1
                      do m1=1, op_v(s+1)%N
                        do m2=1, op_v(s+1)%N
                          call state%set(i)
                          call state%annihil_e( op_v(s2+1)%p(m1)-1, sigma2, s2 )
                          call state%create_e ( op_v(s2+1)%p(m2)-1, sigma2, s2 )
                          call state%annihil_e( op_v(s +1)%p(n1)-1, sigma , s  )
                          call state%create_e ( op_v(s +1)%p(n2)-1, sigma , s  )

                          call this%add_matrixelement(state%get_i(), i, &
                            & state%get_factor() * op_v(s+1)%O(n2,n1) * op_v(s2+1)%O(m2,m1) * &
                            & op_v(s+1)%g * op_v(s2+1)%g / (-dtau) )
                        enddo
                      enddo
                    enddo
                  enddo
                enddo
              enddo
            enddo
          enddo
        endif
      enddo
    end subroutine ed_ham_add_op_v


    subroutine ed_ham_test_hermitian(this)
!--------------------------------------------------------------------
!> @author
!> ALF-project
!>
!> @brief
!> Test if H is hermitian
      IMPLICIT NONE
      class(ed_ham), intent(inout) :: this

      integer :: n, i, j
      real(dp), parameter :: zero = 1.d-10

      print*, "Test hermiticity"
      do n=0, this%N_orbitals*this%N_SUN
        if( .not.this%restrict_N_part .or. n == this%N_part ) then
          do i=1, size(this%H_part(n)%H, 1)
            do j=i, size(this%H_part(n)%H, 1)
              if( abs( this%H_part(n)%H(i,j)- conjg(this%H_part(n)%H(j,i)) ) > zero ) then
                write(error_unit,*) "H", n, i, j, "not hermitian",  this%H_part(n)%H(i,j)- conjg(this%H_part(n)%H(j,i))
                CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
              endif
            enddo
          enddo
        endif
      enddo
      print*, "Hermiticity test concluded"
    end subroutine ed_ham_test_hermitian


    subroutine ed_ham_part_eigenvalues(this)
!--------------------------------------------------------------------
!> @author
!> ALF-project
!>
!> @brief
!> Calculates eigenvalues of Hamiltonian block
      IMPLICIT NONE
      class(ed_ham_part), intent(inout) :: this

      INTEGER               :: INFO, LDA, LWORK
      real(dp), allocatable :: E(:)
      complex(dp), allocatable :: TAU(:), WORK(:)

      !TODO: LWORK can probably be chosen in a more intelligent way
      LWORK = 32 * this%N_states

      allocate( this%eigenval(this%N_states), E(this%N_states-1) )
      allocate( TAU(this%N_states-1), WORK(LWORK) )

      !ZHETRD - reduce a complex Hermitian matrix A to real symmetric tridiagonal
      !         form T by a unitary similarity transformation: Q**H * A * Q = T
      !SUBROUTINE ZHETRD( UPLO, N, A, LDA, D, E, TAU, WORK, LWORK, INFO )
      print*, 'N_states', this%N_states
      CALL ZHETRD( 'U', this%N_states, this%H, this%N_states, this%eigenval, E, TAU, WORK, LWORK, INFO )
      if ( INFO .ne. 0 ) then
        write(error_unit,*) "Error with ZHETRD in ed_ham_part_eigenvalues", INFO
        CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
      endif

      deallocate( this%H, TAU, WORK )

      !DSTERF -  computes all eigenvalues of a symmetric tridiagonal matrix using the
      !          Pal-Walker-Kahan variant of the QL or QR algorithm.
      !SUBROUTINE DSTERF( N, D, E, INFO )
      CALL DSTERF( this%N_states, this%eigenval, E, INFO )
      if ( INFO .ne. 0 ) then
        write(error_unit,*) "Error with DSTERF in ed_ham_part_eigenvalues", INFO
        CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
      endif

      deallocate( E )
    end subroutine ed_ham_part_eigenvalues


    subroutine ed_ham_eigenvalues(this)
!--------------------------------------------------------------------
!> @author
!> ALF-project
!>
!> @brief
!> Calculates eigenvalues of Hamiltonian
        IMPLICIT NONE
        class(ed_ham)  , intent(inout) :: this

        integer :: n, n_temp
        real(dp) :: buff

        print*, "Calculating eigenvalues"


        if( this%restrict_N_part ) then
          allocate( this%eigenval(this%H_part(this%N_part)%N_states) )
        else
          allocate( this%eigenval(this%N_states) )
        endif
        n_temp = 1

        do n=0, this%N_orbitals*this%N_SUN*this%N_FL
          if( .not.this%restrict_N_part .or. n == this%N_part ) then
            print "(I0,' out of ',I0)", n, this%N_orbitals*this%N_SUN*this%N_FL
            call ed_ham_part_eigenvalues(this%H_part(n))
            this%eigenval(n_temp:n_temp+this%H_part(n)%N_states-1) = this%H_part(n)%eigenval(:)
            deallocate( this%H_part(n)%eigenval )
            n_temp = n_temp + this%H_part(n)%N_states
          endif
        enddo

        !Sort this%eigenval
        if( .not.this%restrict_N_part .or. n == this%N_part ) then
          do n=1, this%N_states
            n_temp = minloc( this%eigenval(n:this%N_states), dim=1 ) + n - 1

            buff = this%eigenval(n)
            this%eigenval(n)      = this%eigenval(n_temp)
            this%eigenval(n_temp) = buff
          enddo
        endif

        print*, "Done calculating eigenvalues"
        print*, "Ground state energy:", this%eigenval(1)

        OPEN(Unit=50, file="ED_Eigenvalues", status="replace")
        do n=1, size(this%eigenval)
          write(50,*) this%eigenval(n)
        enddo
        close(50)
    end subroutine ed_ham_eigenvalues

end MODULE ed_ham_mod
