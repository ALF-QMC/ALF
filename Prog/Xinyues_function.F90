      subroutine Get_spin_fluctuation_configuration(Phi,Dtau,Beta,L1)


        ! This subroutine generates a random configuration of spin fluctuation
        ! Phi(L1,L1,3,Ltrot)   Ltrot = Beta/Dtau
        ! Phi(ix,iy,alpha,tau) is the field on time slice  tau  that couples to the sigma^{\alpha} Pauli matrix at 
        ! site \vec{i} = ix*a_1 + iy*a_2 

        implicit none
        real*8::  Dtau,Beta
        Type(Lattice), Intent(In)                         :: Latt
        integer:: Ltrot,L1,Latt%N
        integer, parameter:: nm=Latt%N*Ltrot*3   ! total number of degree of freedom
        integer, parameter:: nh=nm/2      ! number of gaussian pairs
        integer, parameter:: bzl=-L1/2+1  ! lower boundary of BZ
        integer, parameter:: bzh=L1/2     ! upper boundary of BZ
        integer, parameter:: bwl=-Ltrot/2+1  ! lower boundary of wBZ
        integer, parameter:: bwh=Ltrot/2     ! upper boundary of wBZ
       
        complex*16 pwl(L1,L1)
        complex*16 pww(Ltrot,Ltrot)
        real*8,  parameter:: pi=dacos(-1.d0)
        real*8,  parameter:: sq2=dsqrt(2.d0)
        real*8,  parameter:: chi0_inv=1.d0/beta ! set beta*chi0_inv=1
        real*8,  parameter:: xi=3.d0    ! correlatio length
        real*8,  parameter:: wsf=1.d0  ! measured in uint of k_BT
        real*8,  parameter:: factor=Beta*chi0_inv/xi**2/Ltrot
        integer  inl(bzl:bzh),inw(bwl:bwh)
        integer i,j,k,l,inq,inq1,ins
        real*8  gauss(nm)
        real*8  u,v,r,theta,x,y,rq
        real*8  cqx,cqy,wtv,phase,sumr
        integer qx,qy,wt,qx1,qx2
        integer qxb,qyb,wtb,dist
        integer occ(Latt%N*Ltrot)
        complex*16  phiq(Latt%N*Ltrot,3)
        complex*16  phi1(L1,L1*Ltrot)
        complex*16  phi2(Ltrot,Latt%N)
        real*8  Phi(Latt%N*Ltrot,3)

        do i=1,L1
        do qx=bzl,bzh
        inq=(qx-bzl)+1
        phase=(2.d0*pi*qx*i)/L1
        pwl(i,inq)=dcmplx(dcos(phase),dsin(phase))
        end do
        end do

        !generating factors for FT in advance
        !needs to be rearranged into the subroutine.
        do k=1,Ltrot
        do wt=bwl,bwh
        inq=(wt-bwl)+1
        phase=(2.d0*pi*wt*k)/Ltrot
        pww(k,inq)=dcmplx(dcos(phase),dsin(phase))
        end do
        end do


        k=1
        do i=1,nh
        call random_number(u)
        call random_number(v)
        r=dsqrt(-2.d0*dlog(1.d0-u))
        theta=v*2.d0*pi
        x=r*dcos(theta)
        y=r*dsin(theta)
        gauss(k)=x
        gauss(k+1)=y
        k=k+2
        end do

        do i=bzl,bzh
        inl(i)=-i
        end do
        inl(bzh)=bzh

        do i=bwl,bwh
        inw(i)=-i
        end do
        inw(bwh)=bwh

        k=1
        occ=0
        do qx=bzl,bzh
        cqx=dcos(qx*2.d0*pi/L1)
        do qy=bzl,bzh
        cqy=dcos(qy*2.d0*pi/L1)
        do wt=bwl,bwh
        wtv=wt*2.d0*pi

        qxb=inl(qx)
        qyb=inl(qy)
        wtb=inw(wt)
        inq=((qx-bzl)*L1+(qy-bzl))*Ltrot+(wt-bwl)+1
        inq1=((qxb-bzl)*L1+(qyb-bzl))*Ltrot+(wtb-bwl)+1

        dist=iabs(qx-qxb)+iabs(qy-qyb)+iabs(wt-wtb)

        rq=factor*(1.d0+(4.d0+2.d0*(cqx+cqy))*xi**2+dabs(wtv)/wsf)
        rq=dsqrt(rq)

        if(occ(inq).eq.0)then
        if(dist.eq.0)then
        phiq(inq,1)=dcmplx(gauss(k),0.d0)/rq
        phiq(inq,2)=dcmplx(gauss(k+1),0.d0)/rq
        phiq(inq,3)=dcmplx(gauss(k+2),0.d0)/rq
        k=k+3
        occ(inq)=1
        else
        phiq(inq,1)=dcmplx(gauss(k),gauss(k+1))/rq/sq2
        phiq(inq1,1)=dcmplx(gauss(k),-gauss(k+1))/rq/sq2
        phiq(inq,2)=dcmplx(gauss(k+2),gauss(k+3))/rq/sq2
        phiq(inq1,2)=dcmplx(gauss(k+2),-gauss(k+3))/rq/sq2
        phiq(inq,3)=dcmplx(gauss(k+4),gauss(k+5))/rq/sq2
        phiq(inq1,3)=dcmplx(gauss(k+4),-gauss(k+5))/rq/sq2
        k=k+6
        occ(inq)=1
        occ(inq1)=1
        end if
        end if  

        end do 
        end do 
        end do 
        
        do l=1,3
        do qx=bzl,bzh
        qx1=(qx-bzl)+1
        do qy=bzl,bzh
        do wt=bwl,bwh
        qx2=(qy-bzl)*Ltrot+(wt-bwl)+1
        inq=((qx-bzl)*L1+(qy-bzl))*Ltrot+(wt-bwl)+1
        phi1(qx1,qx2)=phiq(inq,l) 
        rq1(qx1,qx2)=rq0(inq)!
        end do
        end do
        end do

        phi1=matmul(pwl,phi1)

        do i=1,L1
        do qy=bzl,bzh
        do wt=bwl,bwh
        qx2=(qy-bzl)*Ltrot+(wt-bwl)+1
        inq=(i-1)*L1*Ltrot+(qy-bzl)*Ltrot+(wt-bwl)+1
        phiq(inq,l)=phi1(i,qx2)
        end do
        end do
        end do

        do qy=bzl,bzh
        qx1=(qy-bzl)+1
        do i=1,L1
        do wt=bwl,bwh
        qx2=(i-1)*Ltrot+(wt-bwl)+1
        inq=(i-1)*L1*Ltrot+(qy-bzl)*Ltrot+(wt-bwl)+1
        phi1(qx1,qx2)=phiq(inq,l) 
        end do
        end do
        end do

        phi1=matmul(pwl,phi1)

        do i=1,L1
        do j=1,L1
        do wt=bwl,bwh
        qx2=(i-1)*Ltrot+(wt-bwl)+1
        inq=(i-1)*L1*Ltrot+(j-1)*Ltrot+(wt-bwl)+1
        phiq(inq,l)=phi1(j,qx2)
        end do
        end do
        end do

        do wt=bwl,bwh
        qx1=(wt-bwl)+1
        do i=1,L1
        do j=1,L1
        qx2=(i-1)*L1+(j-1)+1
        inq=(i-1)*L1*Ltrot+(j-1)*Ltrot+(wt-bwl)+1
        phi2(qx1,qx2)=phiq(inq,l)
        end do
        end do
        end do

        phi2=matmul(pww,phi2)

        do i=1,L1
        do j=1,L1
        do k=1,Ltrot
        qx2=(i-1)*L1+(j-1)+1
        inq=(i-1)*L1*Ltrot+(j-1)*Ltrot+(k-1)+1
        phiq(inq,l)=phi2(k,qx2)
       
        end do
        end do
        end do

        end do

        Phi=dreal(phiq)/dsqrt(1.d0*Latt%N*Ltrot)
        
        
        return

        end


