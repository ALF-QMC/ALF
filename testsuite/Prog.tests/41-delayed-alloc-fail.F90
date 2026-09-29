! A refused panel allocation must stop the run with a message naming the delay.
!
! delay_alloc allocates the panels with stat=, and on failure reports the depth,
! the size requested and ALF_DELAY_K before terminating. Without that guard the
! runtime aborts with a generic allocation error that does not point at the
! delay at all.

Program DelayedAllocFail

   Use delayed_update_mod

   Implicit None

   Integer, Parameter :: Ndim = 2**22, N_FL = 1, dmax = 2

   Call delay_alloc(Ndim, N_FL, dmax, Ndim)

   Write (*,*) "ERROR: a panel allocation of some 560 TB succeeded"
   Call delay_dealloc()
   Stop 2

End Program DelayedAllocFail