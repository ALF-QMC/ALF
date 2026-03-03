!  Copyright (C) 2016-2026 The ALF project
! 
!     The ALF project is free software: you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation, either version 3 of the License, or
!     (at your option) any later version.
! 
!     The ALF project is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.
! 
!     You should have received a copy of the GNU General Public License
!     along with Foobar.  If not, see http://www.gnu.org/licenses/.
!     
!     Under Section 7 of GPL version 3 we require you to fulfill the following additional terms:
!     
!     - It is our hope that this program makes a contribution to the scientific community. Being
!       part of that community we feel that it is reasonable to require you to give an attribution
!       back to the original authors if you have benefitted from this program.
!       Guidelines for a proper citation can be found on the project's homepage
!       http://alf.physik.uni-wuerzburg.de .
!       
!     - We require the preservation of the above copyright notice and this license in all original files.
!     
!     - We prohibit the misrepresentation of the origin of the original source files. To obtain 
!       the original source files please visit the homepage http://alf.physik.uni-wuerzburg.de .
! 
!     - If you make substantial changes to the program we require you to either consider contributing
!       to the ALF project or to mark your material in a reasonable way as different from the original version.


Module Random_Wrap
!--------------------------------------------------------------------
!> @author ALF-project
!> @brief Wrappers for random number generator.
!
!--------------------------------------------------------------------
   contains

!--------------------------------------------------------------------
!> @author ALF-project
!> @brief Query the length of the RNG seed vector.
!> @details Returns the number of integers required to completely specify
!> the state of the Fortran intrinsic random number generator.
!> @param[out] K Integer seed vector length (typically 33 for gfortran)
!--------------------------------------------------------------------
     
     Subroutine Get_seed_Len(K)
       Implicit none
       Integer :: K
       CALL RANDOM_SEED (SIZE=K) 
     end Subroutine Get_seed_Len

!--------------------------------------------------------------------
!> @author ALF-project
!> @brief Initialize RNG state from seed vector.
!> @details Sets up the Fortran intrinsic RNG using provided seed values.
!> If Iseed_vec(1)==0, initializes with default seed 8752143 and uses LCG
!> for remaining seed components. Otherwise distributes provided seeds.
!> @param[in] Iseed_vec Integer array of seed values (dimension K)
!--------------------------------------------------------------------
     
     Subroutine Ranset(Iseed_vec)
       Implicit none
       Integer, Dimension(:) :: Iseed_vec

       Integer :: K, N, i, Iseed
       Integer, allocatable :: Seed_start(:)
       Real (Kind=Kind(0.d0)) :: X

       N = size(Iseed_vec) 
       CALL RANDOM_SEED (SIZE=K)           
       Allocate         (SEED_start(K) ) 
       ! Setup SEED_start
       Iseed = Iseed_vec(1)
       If (Iseed == 0 )  then 
          Iseed = 8752143
          N = 0
       endif
       do i = 1,K
          if (i <= N) then 
             SEED_Start(i) = Iseed_vec(i)
          else
             X = lcg(Iseed)
             SEED_Start(i) = Iseed
          endif
       enddo
       CALL RANDOM_SEED (PUT = SEED_start(1:K)) 
       !Write(6,*) 'Starting seeds ', SEED_Start

     end Subroutine Ranset
       
!--------------------------------------------------------------------
!> @author ALF-project
!> @brief Retrieve current RNG state into seed vector.
!> @details Queries the Fortran intrinsic RNG and extracts current state
!> into the provided seed array for later restoration via Ranset.
!> @param[out] Iseed_vec Integer array to hold current RNG state
!--------------------------------------------------------------------

     Subroutine Ranget(Iseed_vec)
       Implicit none
       Integer, Dimension(:) :: Iseed_vec

       Integer :: K, N, i, Iseed
       Integer, allocatable :: Seed_end(:)
       Real (Kind=Kind(0.d0)) :: X

       N = size(Iseed_vec) 
       CALL RANDOM_SEED (SIZE=K)           
       Allocate         (SEED_end(K) ) 
       CALL RANDOM_SEED (GET = SEED_end(1:K)) 
       ! Setup SEED_start
       Iseed = Iseed_vec(1)
       do i = 1,N
          if (i <= K) then 
             Iseed_vec(i)  =  SEED_end(i)
          else
             X = lcg(Iseed)
             Iseed_vec(i) = Iseed
          endif
       enddo
       !Write(6,*) 'End seeds ', SEED_end

     end Subroutine Ranget

!--------------------------------------------------------------------
!> @author
!> Florian Goth
!
!> @brief 
!> This function returns a real Pseudo-Random Number using a Linear congrential
!> Random number generator. The range of the returned values is [0.0, 1.0).
!> Note that the value 1.0 will not be returned.
!
!> @param[in] seed An integer to seed the LCG.
!-------------------------------------------------------------------- 

     real (Kind=Kind(0.D0)) function lcg(seed)
       use, intrinsic :: iso_fortran_env
       implicit none
       integer :: seed
       integer(INT64) :: res, norm
       
       res = seed ! convert type
       res = 62089911*res + 4349
       norm = 2147483648_INT64 !specify 8 byte integer
       lcg = DBLE(MODULO(res, norm))/DBLE(norm)
       seed = Int(res,kind(0)) ! convert back
     end function lcg

!--------------------------------------------------------------------
!> @author ALF-project
!> @brief Return uniform random variate in [0.0, 1.0).
!> @details Wrapper around Fortran intrinsic RANDOM_NUMBER. Returns
!> uniformly distributed random number in half-open interval [0, 1).
!> @param[in] iq Optional index (not used, for interface consistency)
!> @return Real(8) uniform random number
!--------------------------------------------------------------------

     real (Kind=Kind(0.D0))  function  ranf_wrap(iq)
       implicit none
       integer, optional ::  iq
       Real (Kind=Kind(0.D0)) :: X
       Call Random_Number(X)
       ranf_wrap = X
     end function ranf_wrap
     
!--------------------------------------------------------------------
!> @author ALF-project
!> @brief Return standard normal (Gaussian) random variate.
!> @details Generates standard normal N(0,1) variate using Box-Muller
!> transformation: R = sqrt(-2*ln(u1)), θ = 2π*u2, x = R*cos(θ)
!> where u1, u2 are independent uniform [0,1) variates.
!> @param[in] iq Optional index (not used, for interface consistency)
!> @return Real(8) standard normal random variable
!--------------------------------------------------------------------

     real (kind=kind(0.D0))  function  rang_wrap(iq)
        
       ! Random variable according to the distribution:  exp(-x**2/2)/(sqrt(2*3.1415927))
       Implicit none
     
       integer, optional :: iq
       real (Kind=kind(0.D0)) ::  pi, ranmod, theta
       
       PI = 3.1415926536D0
       RANMOD = SQRT(-2.D0 * log(RANF_Wrap(iq)))
       THETA  = 2.D0 * PI * RANF_wrap(iq)
       rang_wrap = RANMOD * COS(THETA)
       
     end function rang_wrap
     
!--------------------------------------------------------------------
!> @author ALF-project
!> @brief Return uniform random integer in [1, N].
!> @details Generates uniformly distributed random integer in inclusive
!> range [1, N] by scaling and rounding a uniform [0,1) variate.
!> Values are clamped to [1, N] to handle rounding edge cases.
!> @param[in] N Upper bound (inclusive) for integer range
!> @return Integer in [1, N]
!--------------------------------------------------------------------

     integer function nranf(N)
       implicit none
       integer :: N
       
       nranf  = nint(ranf_wrap()*dble(N) + 0.5D0)
       
       if (nranf .lt. 1 ) nranf = 1
       if (nranf .gt. N ) nranf = N 
       
     end function nranf
     
   end Module Random_Wrap
