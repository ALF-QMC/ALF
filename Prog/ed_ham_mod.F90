MODULE ed_ham_mod

    use Operator_mod
    use ed_state_mod
    
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: ed_ham

    TYPE ed_ham
!--------------------------------------------------------------------
!> @author 
!> ALF-project
!>   
!> @brief 
!> Defines a basis spanning the Fock space, for exact diagonalisation. 
        private
        INTEGER :: N_orbitals, N_SUN, N_states
        complex(kind=kind(0.d0)), allocatable :: H(:,:)
        real(kind=kind(0.d0))   , allocatable :: eigenval(:)
        CONTAINS
            private
            PROCEDURE :: init => ed_ham_init
            PROCEDURE :: add_op_t => ed_ham_add_op_t
            PROCEDURE :: add_op_v => ed_ham_add_op_v
            
            PROCEDURE, public :: build_h => ed_ham_build_h
            PROCEDURE, public :: energy => ed_ham_energy
            PROCEDURE, public :: get_eigenvalues => ed_ham_get_eigenvalues
    END TYPE ed_ham

CONTAINS


    subroutine ed_ham_init(this, N_orbitals, N_SUN)
        IMPLICIT NONE
        class(ed_ham), INTENT(INOUT) :: this
        integer, intent(in) :: N_orbitals, N_SUN
        
        this%N_orbitals = N_orbitals
        this%N_SUN      = N_SUN
        this%N_states   = 2**(N_orbitals*N_SUN)
        allocate( this%H(this%N_states, this%N_states) )
        this%H(:,:) = cmplx(0.d0, 0.d0, kind=kind(0.d0))
    
    end subroutine ed_ham_init
    

    subroutine ed_ham_add_op_t(this, op_t, dtau)
        IMPLICIT NONE
        class(ed_ham)  , intent(inout) :: this
        type(Operator), intent(in)    :: Op_T
        real(Kind=Kind(0.d0)), intent(in) :: dtau
        
        integer :: i, n1, n2, s
        type(ed_state) :: state
        
        call state%init(this%N_orbitals)
        
        do i=0, this%N_states-1
            do s=0, this%N_SUN-1
                do n1=1, op_t%N
                    do n2=1, op_t%N
                        call state%set(i)
                        call state%annihil_e( op_t%p(n1)-1, s )
                        call state%create_e ( op_t%p(n2)-1, s )
                        
                        this%H(state%get_i()+1, i+1) = this%H(state%get_i()+1, i+1) &
                            & + state%get_factor() * op_t%O(n2,n1) * op_t%g / (-dtau)
                    enddo
                enddo
            enddo
        enddo
    
    end subroutine ed_ham_add_op_t
    

    subroutine ed_ham_add_op_v(this, op_v, dtau)
        IMPLICIT NONE
        class(ed_ham)  , intent(inout) :: this
        type(Operator), intent(in)    :: op_v
        real(Kind=Kind(0.d0)), intent(in) :: dtau
        
        integer :: i, n1, n2, m1, m2, s, s2
        type(ed_state) :: state
        
        if( op_v%type .ne. 2 ) then
            print*, "ED only implemented for OP_V type 2"
            stop
        endif
        
        call state%init(this%N_orbitals)
        
        do i=0, this%N_states-1
            this%H(i+1, i+1) = this%H(i+1, i+1) &
            & + op_v%g**2 * op_v%alpha**2 * &
            &   cmplx( - real(this%N_SUN**2,kind=kind(0.d0)) / dtau, 0.d0, kind=kind(0.d0) )
            do s=0, this%N_SUN-1
                do n1=1, op_v%N
                    do n2=1, op_v%N
                        call state%set(i)
                        call state%annihil_e( op_v%p(n1)-1, s )
                        call state%create_e ( op_v%p(n2)-1, s )
                        
                        this%H(state%get_i()+1, i+1) = this%H(state%get_i()+1, i+1) + state%get_factor() * &
                        & op_v%O(n2,n1) * op_v%g**2 * op_v%alpha * cmplx( - real(2*this%N_SUN,kind=kind(0.d0)) / dtau, 0.d0, kind=kind(0.d0) )
                        
                        do s2=0, this%N_SUN-1
                            do m1=1, op_v%N
                                do m2=1, op_v%N
                                    call state%set(i)
                                    call state%annihil_e( op_v%p(m1)-1, s2 )
                                    call state%create_e ( op_v%p(m2)-1, s2 )
                                    call state%annihil_e( op_v%p(n1)-1, s  )
                                    call state%create_e ( op_v%p(n2)-1, s  )
                                    
                                    this%H(state%get_i()+1, i+1) = this%H(state%get_i()+1, i+1) &
                                    & + state%get_factor() * op_v%O(n2,n1) * op_v%O(m2,m1) * op_v%g**2 / (-dtau)
                                enddo
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        enddo
    
    end subroutine ed_ham_add_op_v
    
          
    subroutine ed_ham_build_h(this, ndim, N_SUN, OP_T, OP_V, dtau)
        IMPLICIT NONE
        class(ed_ham)  , intent(inout) :: this
        integer, intent(in) :: ndim, N_SUN
        Type (Operator), intent(in) :: Op_V(:,:)
        Type (Operator), intent(in) :: Op_T(:,:)
        real(Kind=Kind(0.d0)), intent(in) :: dtau
        
        
        INTEGER :: nf, n
        print*, "Building ED-Hamiltonian"
        call this%init(ndim, N_SUN)
        do nf=1,1
            do n=1, Size(OP_T,1)
                call this%add_op_t( OP_T(n,nf), dtau )
            enddo
            do n=1, Size(OP_V,1)
                call this%add_op_v( OP_V(n,nf), dtau )
            enddo
        enddo
        print*, "done Building ED-Hamiltonian"
        call ed_ham_test_hermitian(this)
        call ed_ham_eigenvalues(this)
    
    end subroutine ed_ham_build_h
    
          
    subroutine ed_ham_test_hermitian(this)
        !Test if H is hermitian
        IMPLICIT NONE
        class(ed_ham)  , intent(inout) :: this
        
        integer :: i, j
        real(Kind=Kind(0.d0)) :: zero = 1.d-10
        
        print*, "Test hermiticity"
        !do i=1, size(this%H, 1)
        !    print*, this%H(:,i)
        !enddo
        do i=1, size(this%H, 1)
            do j=i, size(this%H, 1)
                if( abs( this%H(i,j)- conjg(this%H(j,i)) ) > zero ) then
                    print*, "H", i, j, "not hermitian",  this%H(i,j)- conjg(this%H(j,i))
                    stop 1
                endif
            enddo
        enddo
        print*, "Hermiticity test concluded"

    
    end subroutine ed_ham_test_hermitian
    
          
    subroutine ed_ham_eigenvalues(this)
        IMPLICIT NONE
        class(ed_ham)  , intent(inout) :: this
        
        INTEGER               :: INFO, LDA, LWORK, i
        real(kind=kind(0.d0)), allocatable :: E(:)
        complex(kind=kind(0.d0)), allocatable :: TAU(:), WORK(:)
        
        print*, "Calculating eigenvalues"
        
        !TODO: LWORK can probably be chosen in a more intelligent way
        LWORK = 32 * this%N_states 
        
        allocate( this%eigenval(this%N_states), E(this%N_states-1) )
        allocate( TAU(this%N_states-1), WORK(LWORK) )
        
        !ZHETRD - reduce a complex Hermitian matrix A to real symmetric tridiagonal
        !         form T by a unitary similarity transformation: Q**H * A * Q = T
        !SUBROUTINE ZHETRD( UPLO, N, A, LDA, D, E, TAU, WORK, LWORK, INFO )
        CALL ZHETRD( 'U', this%N_states, this%H, this%N_states, this%eigenval, E, TAU, WORK, LWORK, INFO )
        if ( INFO .ne. 0 ) then
            print*, "Error with ZHETRD in ed_ham_eigenvalues", INFO
            stop 
        endif
        
        !print*, LWORK, WORK(1)
        
        deallocate( this%H, TAU, WORK )
        
        !DSTERF -  computes all eigenvalues of a symmetric tridiagonal matrix using the
        !          Pal-Walker-Kahan variant of the QL or QR algorithm.
        !SUBROUTINE DSTERF( N, D, E, INFO )
        CALL DSTERF( this%N_states, this%eigenval, E, INFO )
        if ( INFO .ne. 0 ) then
            print*, "Error with DSTERF in ed_ham_eigenvalues", INFO
            stop 
        endif
        
        deallocate( E )
        
        print*, "Done calculating eigenvalues"
        print*, "Ground state energy:", this%eigenval(1)
          
        OPEN(Unit = 50,file="ED_Eigenvalues",status="replace")
        do i=1, this%N_states
            write(50,*) this%eigenval(i)
        enddo
        close(50)
        
    end subroutine ed_ham_eigenvalues
    
          
    function ed_ham_energy(this, beta)
        IMPLICIT NONE
        class(ed_ham)  , intent(inout) :: this
        real(kind=kind(0.d0)), intent(in) :: beta
        
        INTEGER               :: i
        real(kind=kind(0.d0)) :: Z, E, ed_ham_energy
        
        Z = 0.d0
        E = 0.d0
        
        do i=1, this%N_states
            Z = Z + exp(-beta * this%eigenval(i))
            E = E + exp(-beta * this%eigenval(i)) * this%eigenval(i)
        enddo
        
        ed_ham_energy = E / Z
        
    end function ed_ham_energy
    
          
    subroutine ed_ham_get_eigenvalues(this, eigenvalues)
        IMPLICIT NONE
        class(ed_ham)  , intent(inout) :: this
        
        real(kind=kind(0.d0)), allocatable, intent(out) :: eigenvalues(:)
        
        allocate( eigenvalues(this%N_states) )
        
        eigenvalues = this%eigenval
        
    end subroutine ed_ham_get_eigenvalues
    
end MODULE ed_ham_mod
