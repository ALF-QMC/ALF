

    module types

        implicit none
        integer, parameter :: dp = kind(1.d0) ! double precision

    end module

    module lattices_f2py
        use lattices_v3
        use types
        
        implicit none
        Integer :: la_N, la_Ns
        Integer, allocatable :: la_list(:,:),  la_invlist(:,:),  la_nnlist(:,:,:), &
                &               la_listk(:,:), la_invlistk(:,:), la_nnlistk(:,:,:), la_imj(:,:)
        Real(dp) , allocatable :: la_b1_p(:), la_b2_p(:), la_BZ1_p(:), la_BZ2_p(:), &
                &               la_b1_perp_p(:), la_b2_perp_p(:)
        
        contains
        
        subroutine Lattice_out(L1_p, L2_p, a1_p, a2_p)

            Implicit none

            Real(dp), dimension(:), intent(in) :: L1_p, L2_p, a1_p, a2_p
            
            Type (Lattice) :: Latt 
            Integer :: ndim, N, L
            
            call Make_lattice(L1_p, L2_p, a1_p, a2_p, Latt)
            
            ndim = size(L1_p)
            L  = (size(Latt%Invlist,1) - 1)/2
            
            N  = Latt%N
            
            allocate ( la_b1_p(ndim), la_b2_p(ndim), la_BZ1_p(ndim), la_BZ2_p(ndim) )
            allocate ( la_b1_perp_p(ndim), la_b2_perp_p(ndim) )
            Allocate ( la_List(N,ndim), la_Invlist(-L:L, -L:L ) )
            Allocate ( la_Listk(N,ndim), la_Invlistk(-L:L, -L:L) )
            Allocate ( la_nnlist(N,-1:1,-1:1), la_nnlistk(N,-1:1,-1:1) )
            Allocate ( la_imj(N,N) )
            
            la_N         = Latt%N
            la_Ns        = Latt%Ns
            la_list      = Latt%list
            la_invlist   = Latt%invlist
            la_nnlist    = Latt%nnlist
            la_listk     = Latt%listk
            la_invlistk  = Latt%invlistk
            la_nnlistk   = Latt%nnlistk
            la_imj       = Latt%imj
            la_b1_p      = Latt%b1_p
            la_b2_p      = Latt%b2_p
            la_BZ1_p     = Latt%BZ1_p
            la_BZ2_p     = Latt%BZ2_p
            la_b1_perp_p = Latt%b1_perp_p
            la_b2_perp_p = Latt%b2_perp_p
            
            call clear_Lattice(Latt)
            
        end subroutine Lattice_out
        
        subroutine Lattice_out_clean()

            Implicit none
            
            deallocate ( la_b1_p, la_b2_p, la_BZ1_p, la_BZ2_p )
            deallocate ( la_b1_perp_p, la_b2_perp_p )
            deallocate ( la_List, la_Invlist )
            deallocate ( la_Listk, la_Invlistk )
            deallocate ( la_nnlist, la_nnlistk )
            deallocate ( la_imj )
        
        end subroutine Lattice_out_clean
        
    end module lattices_f2py
