#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module viscosity_base
!
!  This module should be renamed "global_deriv_mod.F90"
! 
!  It is a collection of derivative operators that must be applied to the field 
!  over the sphere (as opposed to derivative operators that can be applied element 
!  by element)
!
!
use thread_mod, only : omp_get_num_threads
use kinds, only : real_kind, iulog
use dimensions_mod, only : np, nlev,qsize,nelemd
use hybrid_mod, only : hybrid_t, hybrid_create
use parallel_mod, only : parallel_t, abortmp
use element_mod, only : element_t
use derivative_mod, only : derivative_t, laplace_sphere_wk, vlaplace_sphere_wk, vorticity_sphere, derivinit, divergence_sphere
use edgetype_mod, only : EdgeBuffer_t, EdgeDescriptor_t
use edge_mod, only : edgevpack, edgevunpack, edgevunpackmin, &
    edgevunpackmax, initEdgeBuffer, FreeEdgeBuffer, edgeSunpackmax, edgeSunpackmin,edgeSpack, &
    edgevpack_nlyr, edgevunpack_nlyr, edge_g

use bndry_mod, only : bndry_exchangev, bndry_exchangeS, bndry_exchangeS_start,bndry_exchangeS_finish
use control_mod, only : hypervis_scaling, nu, nu_div
use perf_mod, only: t_startf, t_stopf

implicit none
private
save

public :: neighbor_minmax
#ifdef _PRIM
public :: biharmonic_wk_scalar
public :: neighbor_minmax_start,neighbor_minmax_finish
public :: smooth_phis
#endif


!
! compute vorticity/divergence and then project to make continious
! high-level routines uses only for I/O
public :: compute_zeta_C0
public :: compute_div_C0
interface compute_zeta_C0
    module procedure compute_zeta_C0_hybrid       ! hybrid version
    module procedure compute_zeta_C0_par          ! single threaded
end interface
interface compute_div_C0
    module procedure compute_div_C0_hybrid
    module procedure compute_div_C0_par
end interface
interface make_c0
    module procedure make_c0_hybrid
    module procedure make_c0_hybrid1
    module procedure make_c0_par
end interface

public :: make_c0
public :: make_c0_vector
public :: compute_zeta_C0_contra    ! for older versions of sweq which carry
public :: compute_div_C0_contra     ! velocity around in contra-coordinates
public :: compute_eta_C0_contra

type (EdgeBuffer_t)          :: edge1

contains



#ifdef _PRIM

subroutine biharmonic_wk_scalar(elem,qtens,deriv,edgeq,hybrid,nets,nete)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! compute weak biharmonic operator
!    input:  qtens = Q
!    output: qtens = weak biharmonic of Q
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
type (hybrid_t)      , intent(in) :: hybrid
type (element_t)     , intent(inout), target :: elem(:)
integer :: nets,nete
real (kind=real_kind), dimension(np,np,nlev,qsize,nets:nete) :: qtens
type (EdgeBuffer_t)  , intent(inout) :: edgeq
type (derivative_t)  , intent(in) :: deriv

! local
integer :: k,kptr,i,j,ie,ic,q
real (kind=real_kind), dimension(np,np) :: lap_p
logical var_coef1

   !if tensor hyperviscosity with tensor V is used, then biharmonic operator is (\grad\cdot V\grad) (\grad \cdot \grad) 
   !so tensor is only used on second call to laplace_sphere_wk
   var_coef1 = .true.
   if(hypervis_scaling > 0)    var_coef1 = .false.



   do ie=nets,nete
#if (defined COLUMN_OPENMP)
!$omp parallel do private(k, q, lap_p)
#endif
      do q=1,qsize      
         do k=1,nlev    !  Potential loop inversion (AAM)
            lap_p(:,:)=qtens(:,:,k,q,ie)
! Original use of qtens on left and right hand sides caused OpenMP errors (AAM)
           qtens(:,:,k,q,ie)=laplace_sphere_wk(lap_p,deriv,elem(ie),var_coef=var_coef1)
         enddo
         call edgeVpack_nlyr(edgeq, elem(ie)%desc, qtens(:,:,:,q,ie),nlev,nlev*(q-1),qsize*nlev)
      enddo
   enddo

   call t_startf('biwksc_bexchV')
   call bndry_exchangeV(hybrid,edgeq)
   call t_stopf('biwksc_bexchV')
   
   do ie=nets,nete

      ! apply inverse mass matrix, then apply laplace again
#if (defined COLUMN_OPENMP)
!$omp parallel do private(k, q, lap_p)
#endif
      do q=1,qsize      
        call edgeVunpack_nlyr(edgeq,elem(ie)%desc,qtens(:,:,:,q,ie),nlev,nlev*(q-1),qsize*nlev)
        do k=1,nlev    !  Potential loop inversion (AAM)
           lap_p(:,:)=elem(ie)%rspheremp(:,:)*qtens(:,:,k,q,ie)
           qtens(:,:,k,q,ie)=laplace_sphere_wk(lap_p,deriv,elem(ie),var_coef=.true.)
        enddo
      enddo
   enddo
#ifdef DEBUGOMP
#if (defined HORIZ_OPENMP)
!$OMP BARRIER
#endif
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end subroutine

#endif



! wrapper functions
subroutine make_C0_par(zeta,elem,par)
type (parallel_t), intent(in) :: par
type (element_t)     , intent(in), target :: elem(:)
real (kind=real_kind), dimension(np,np,nlev,nelemd) :: zeta
! local
type (hybrid_t)   :: hybrid
hybrid = hybrid_create(par,0,1)
call make_c0_hybrid_klev(zeta,elem,hybrid,1,nelemd,nlev)
end subroutine

subroutine make_C0_hybrid(zeta,elem,hybrid,nets,nete)
type (hybrid_t)      , intent(in) :: hybrid
type (element_t)     , intent(in), target :: elem(:)
integer :: nets,nete
real (kind=real_kind), dimension(np,np,nlev,nets:nete) :: zeta
call make_c0_hybrid_klev(zeta,elem,hybrid,nets,nete,nlev)
end subroutine

subroutine make_C0_hybrid1(zeta,elem,hybrid,nets,nete)
type (hybrid_t)      , intent(in) :: hybrid
type (element_t)     , intent(in), target :: elem(:)
integer :: nets,nete
real (kind=real_kind), dimension(np,np,nets:nete) :: zeta
call make_c0_hybrid_klev(zeta,elem,hybrid,nets,nete,1)
end subroutine


subroutine make_C0_vector(v,elem,hybrid,nets,nete)
type (hybrid_t)      , intent(in) :: hybrid
type (element_t)     , intent(in), target :: elem(:)
integer :: nets,nete
real (kind=real_kind), dimension(np,np,2,nets:nete) :: v
call make_c0_hybrid_klev(v,elem,hybrid,nets,nete,2)
end subroutine





subroutine make_C0_hybrid_klev(zeta,elem,hybrid,nets,nete,klev)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! apply DSS (aka assembly procedure) to zeta.  
! this is a low-performance routine used for I/O and analysis.
! no need to optimize
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
type (hybrid_t)      , intent(in) :: hybrid
type (element_t)     , intent(in), target :: elem(:)
integer :: nets,nete,klev
real (kind=real_kind), dimension(np,np,klev,nets:nete) :: zeta

! local
integer :: k,i,j,ie,ic

do ie=nets,nete
   do k=1,klev
      zeta(:,:,k,ie)=zeta(:,:,k,ie)*elem(ie)%spheremp(:,:)
   enddo
   call edgeVpack_nlyr(edge_g,elem(ie)%desc,zeta(1,1,1,ie),klev,0,klev)
enddo
call bndry_exchangeV(hybrid,edge_g)
do ie=nets,nete
   call edgeVunpack_nlyr(edge_g,elem(ie)%desc,zeta(1,1,1,ie),klev,0,klev)

   do k=1,klev
      zeta(:,:,k,ie)=zeta(:,:,k,ie)*elem(ie)%rspheremp(:,:)
   enddo
enddo
end subroutine












subroutine compute_zeta_C0_contra(zeta,elem,par,nt)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! compute C0 vorticity.  That is, solve:  
!     < PHI, zeta > = <PHI, curl(elem%state%v >
!
!    input:  v (stored in elem()%, in contra-variant coordinates)
!    output: zeta(:,:,:,:)   
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

type (parallel_t)      , intent(in) :: par
type (element_t)     , intent(in), target :: elem(:)
integer :: nt
real (kind=real_kind), dimension(np,np,nlev,nelemd) :: zeta
real (kind=real_kind), dimension(np,np,2) :: ulatlon
real (kind=real_kind), dimension(np,np) :: v1,v2

! local
integer :: k,ie
type (derivative_t)          :: deriv

call derivinit(deriv)

do k=1,nlev
do ie=1,nelemd
    v1 = elem(ie)%state%v(:,:,1,k,nt)
    v2 = elem(ie)%state%v(:,:,2,k,nt)
    ulatlon(:,:,1) = elem(ie)%D(:,:,1,1)*v1 + elem(ie)%D(:,:,1,2)*v2
    ulatlon(:,:,2) = elem(ie)%D(:,:,2,1)*v1 + elem(ie)%D(:,:,2,2)*v2
   !    zeta(:,:,k,ie)=elem(ie)%state%zeta(:,:,k)
   zeta(:,:,k,ie)=vorticity_sphere(ulatlon,deriv,elem(ie))
enddo
enddo

call make_C0(zeta,elem,par)

end subroutine

subroutine compute_eta_C0_contra(eta,elem,par,nt)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! compute C0 absolute vorticity.  That is, solve:
!     < PHI, eta > = <PHI, curl(elem%state%v) + coriolis >
!
!    input:  v (stored in elem()%, in contra-variant coordinates)
!    output: zeta(:,:,:,:)   
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

type (parallel_t)      , intent(in) :: par
type (element_t)     , intent(in), target :: elem(:)
integer :: nt
real (kind=real_kind), dimension(np,np,nlev,nelemd) :: eta
real (kind=real_kind), dimension(np,np,2) :: ulatlon
real (kind=real_kind), dimension(np,np) :: v1,v2

! local
integer :: k,ie
type (derivative_t)          :: deriv

call derivinit(deriv)

do k=1,nlev
do ie=1,nelemd
    v1 = elem(ie)%state%v(:,:,1,k,nt)
    v2 = elem(ie)%state%v(:,:,2,k,nt)
    ulatlon(:,:,1) = elem(ie)%D(:,:,1,1)*v1 + elem(ie)%D(:,:,1,2)*v2
    ulatlon(:,:,2) = elem(ie)%D(:,:,2,1)*v1 + elem(ie)%D(:,:,2,2)*v2
   eta(:,:,k,ie)=vorticity_sphere(ulatlon,deriv,elem(ie)) + elem(ie)%fcor(:,:)
enddo
enddo

call make_C0(eta,elem,par)

end subroutine

subroutine compute_div_C0_contra(zeta,elem,par,nt)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! compute C0 divergence. That is, solve:  
!     < PHI, zeta > = <PHI, div(elem%state%v >
!
!    input:  v (stored in elem()%, in contra-variant coordinates)
!    output: zeta(:,:,:,:)   
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

type (parallel_t)      , intent(in) :: par
type (element_t)     , intent(in), target :: elem(:)
integer :: nt
real (kind=real_kind), dimension(np,np,nlev,nelemd) :: zeta
real (kind=real_kind), dimension(np,np,2) :: ulatlon
real (kind=real_kind), dimension(np,np) :: v1,v2

! local
integer :: k,ie
type (derivative_t)          :: deriv

call derivinit(deriv)

do k=1,nlev
do ie=1,nelemd
    v1 = elem(ie)%state%v(:,:,1,k,nt)
    v2 = elem(ie)%state%v(:,:,2,k,nt)
    ulatlon(:,:,1) = elem(ie)%D(:,:,1,1)*v1 + elem(ie)%D(:,:,1,2)*v2
    ulatlon(:,:,2) = elem(ie)%D(:,:,2,1)*v1 + elem(ie)%D(:,:,2,2)*v2
   !    zeta(:,:,k,ie)=elem(ie)%state%zeta(:,:,k)
   zeta(:,:,k,ie)=divergence_sphere(ulatlon,deriv,elem(ie))
enddo
enddo

call make_C0(zeta,elem,par)

end subroutine

subroutine compute_zeta_C0_par(zeta,elem,par,nt)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! compute C0 vorticity.  That is, solve:  
!     < PHI, zeta > = <PHI, curl(elem%state%v >
!
!    input:  v (stored in elem()%, in lat-lon coordinates)
!    output: zeta(:,:,:,:)   
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
type (parallel_t) :: par
type (element_t)     , intent(in), target :: elem(:)
real (kind=real_kind), dimension(np,np,nlev,nelemd) :: zeta
integer :: nt

! local
type (hybrid_t)              :: hybrid
integer :: k,i,j,ie,ic
type (derivative_t)          :: deriv

! single thread
hybrid = hybrid_create(par,0,1)

call compute_zeta_C0_hybrid(zeta,elem,hybrid,1,nelemd,nt)

end subroutine


subroutine compute_div_C0_par(zeta,elem,par,nt)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! compute C0 divergence. That is, solve:  
!     < PHI, zeta > = <PHI, div(elem%state%v >
!
!    input:  v (stored in elem()%, in lat-lon coordinates)
!    output: zeta(:,:,:,:)   
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

type (parallel_t) :: par
type (element_t)     , intent(in), target :: elem(:)
real (kind=real_kind), dimension(np,np,nlev,nelemd) :: zeta
integer :: nt

! local
type (hybrid_t)              :: hybrid
integer :: k,i,j,ie,ic
type (derivative_t)          :: deriv

! single thread
hybrid = hybrid_create(par,0,1)

call compute_div_C0_hybrid(zeta,elem,hybrid,1,nelemd,nt)

end subroutine



subroutine compute_zeta_C0_hybrid(zeta,elem,hybrid,nets,nete,nt)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! compute C0 vorticity.  That is, solve:  
!     < PHI, zeta > = <PHI, curl(elem%state%v >
!
!    input:  v (stored in elem()%, in lat-lon coordinates)
!    output: zeta(:,:,:,:)   
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

type (hybrid_t)      , intent(in) :: hybrid
type (element_t)     , intent(in), target :: elem(:)
integer :: nt,nets,nete
real (kind=real_kind), dimension(np,np,nlev,nets:nete) :: zeta

! local
integer :: k,i,j,ie,ic
type (derivative_t)          :: deriv

call derivinit(deriv)

do ie=nets,nete
#if (defined COLUMN_OPENMP)
!$omp parallel do private(k)
#endif
do k=1,nlev
   !    zeta(:,:,k,ie)=elem(ie)%state%zeta(:,:,k)
   zeta(:,:,k,ie)=vorticity_sphere(elem(ie)%state%v(:,:,:,k,nt),deriv,elem(ie))
enddo
enddo

call make_C0(zeta,elem,hybrid,nets,nete)

end subroutine


subroutine compute_div_C0_hybrid(zeta,elem,hybrid,nets,nete,nt)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! compute C0 divergence. That is, solve:  
!     < PHI, zeta > = <PHI, div(elem%state%v >
!
!    input:  v (stored in elem()%, in lat-lon coordinates)
!    output: zeta(:,:,:,:)   
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

type (hybrid_t)      , intent(in) :: hybrid
type (element_t)     , intent(in), target :: elem(:)
integer :: nt,nets,nete
real (kind=real_kind), dimension(np,np,nlev,nets:nete) :: zeta

! local
integer :: k,i,j,ie,ic
type (derivative_t)          :: deriv

call derivinit(deriv)

do ie=nets,nete
#if (defined COLUMN_OPENMP)
!$omp parallel do private(k)
#endif
do k=1,nlev
   !    zeta(:,:,k,ie)=elem(ie)%state%zeta(:,:,k)
   zeta(:,:,k,ie)=divergence_sphere(elem(ie)%state%v(:,:,:,k,nt),deriv,elem(ie))
enddo
enddo

call make_C0(zeta,elem,hybrid,nets,nete)

end subroutine








#ifdef _PRIM

subroutine neighbor_minmax(hybrid,edgeMinMax,nets,nete,min_neigh,max_neigh)
 
   type (hybrid_t)      , intent(in) :: hybrid
   type (EdgeBuffer_t)  , intent(inout) :: edgeMinMax
   integer :: nets,nete
   real (kind=real_kind) :: min_neigh(nlev,qsize,nets:nete)
   real (kind=real_kind) :: max_neigh(nlev,qsize,nets:nete)

   ! local 
   integer :: ie,q, k,kptr

   
   do ie=nets,nete
      kptr = 0
      call  edgeSpack(edgeMinMax,min_neigh(:,:,ie),qsize*nlev,kptr,2*qsize*nlev,ie)
      kptr = qsize*nlev
      call  edgeSpack(edgeMinMax,max_neigh(:,:,ie),qsize*nlev,kptr,2*qsize*nlev,ie)
   enddo
   
   call t_startf('nmm_bexchV')
   call bndry_exchangeS(hybrid,edgeMinMax)
   call t_stopf('nmm_bexchV')

   do ie=nets,nete
      kptr = 0
      call  edgeSunpackMIN(edgeMinMax,min_neigh(:,:,ie),qsize*nlev,kptr,2*qsize*nlev,ie)
      kptr = qsize*nlev
      call  edgeSunpackMAX(edgeMinMax,max_neigh(:,:,ie),qsize*nlev,kptr,2*qsize*nlev,ie)
      do q=1,qsize
      do k=1,nlev
          min_neigh(k,q,ie) = max(min_neigh(k,q,ie),0d0)
      enddo
      enddo
   enddo
  
end subroutine neighbor_minmax

subroutine neighbor_minmax_start(hybrid,edgeMinMax,nets,nete,min_neigh,max_neigh)

   type (hybrid_t)      , intent(in) :: hybrid
   type (EdgeBuffer_t)  , intent(inout) :: edgeMinMax
   integer :: nets,nete
   real (kind=real_kind) :: min_neigh(nlev,qsize,nets:nete)
   real (kind=real_kind) :: max_neigh(nlev,qsize,nets:nete)

   ! local 
   integer :: ie,q, k,kptr


   do ie=nets,nete
      kptr = 0
      call  edgeSpack(edgeMinMax,min_neigh(:,:,ie),qsize*nlev,kptr,2*qsize*nlev,ie)
      kptr = qsize*nlev
      call  edgeSpack(edgeMinMax,max_neigh(:,:,ie),qsize*nlev,kptr,2*qsize*nlev,ie)
   enddo

   call t_startf('nmm_bexchS_start')
   call bndry_exchangeS_start(hybrid,edgeMinMax)
   call t_stopf('nmm_bexchS_start')

end subroutine neighbor_minmax_start
subroutine neighbor_minmax_finish(hybrid,edgeMinMax,nets,nete,min_neigh,max_neigh)

   type (hybrid_t)      , intent(in) :: hybrid
   type (EdgeBuffer_t)  , intent(inout) :: edgeMinMax
   integer :: nets,nete
   real (kind=real_kind) :: min_neigh(nlev,qsize,nets:nete)
   real (kind=real_kind) :: max_neigh(nlev,qsize,nets:nete)

   ! local 
   integer :: ie,q, k,kptr

   call t_startf('nmm_bexchS_fini')
   call bndry_exchangeS_finish(hybrid,edgeMinMax)
   call t_stopf('nmm_bexchS_fini')

   do ie=nets,nete
      kptr = 0
      call  edgeSunpackMIN(edgeMinMax,min_neigh(:,:,ie),qsize*nlev,kptr,2*qsize*nlev,ie)
      kptr = qsize*nlev
      call  edgeSunpackMAX(edgeMinMax,max_neigh(:,:,ie),qsize*nlev,kptr,2*qsize*nlev,ie)
      do q=1,qsize
      do k=1,nlev
          min_neigh(k,q,ie) = max(min_neigh(k,q,ie),0d0)
      enddo
      enddo
   enddo

end subroutine neighbor_minmax_finish

subroutine smooth_phis(phis,elem,hybrid,deriv,nets,nete,minf,numcycle,p2filt,xgll)
  use control_mod, only : smooth_phis_nudt, hypervis_scaling
  implicit none

  integer :: nets,nete
  real (kind=real_kind), dimension(np,np,nets:nete), intent(inout)   :: phis
  type (hybrid_t)      , intent(in) :: hybrid
  type (element_t)     , intent(inout), target :: elem(:)
  type (derivative_t)  , intent(in) :: deriv
  real (kind=real_kind), intent(in)   :: minf
  integer,               intent(in) :: numcycle ! number of Laplace smoothing steps
  integer,               intent(in) :: p2filt   ! apply p2 filter 1=before  2=before and after each iteration
  real (kind=real_kind)             :: xgll(np)

  ! local
  type (EdgeBuffer_t)          :: edgebuf
  real (kind=real_kind), dimension(np,np,nets:nete) :: pstens
  real (kind=real_kind), dimension(nets:nete) :: pmin,pmax
  real (kind=real_kind) :: phis4(np)
  integer :: nt,ie,ic,i,j
  integer :: minmax_halo =-1   ! -1 = disabled.  
                               ! 0  = recompute each time

  if (p2filt>=1 .and. np/=4) then
     call abortmp('ERROR: topo smoothing p2 filter option only supported with np==4')
  endif

  ! create edge buffer for 1 field
  call initEdgeBuffer(hybrid%par,edgebuf,elem,1)


  ! compute local element neighbor min/max
  do ie=nets,nete
     pmin(ie)=minval(phis(:,:,ie))
     pmax(ie)=maxval(phis(:,:,ie))
  enddo

  do ic=1,minmax_halo    ! take the min/max over three element halo
  do ie=nets,nete
     pstens(:,:,ie)=pmin(ie)
     call edgeVpack(edgebuf,pstens(:,:,ie),1,0,ie)
  enddo
  call bndry_exchangeV(hybrid,edgebuf)
  do ie=nets,nete
     call edgeVunpackMin(edgebuf, pstens(:,:,ie), 1, 0, ie)
     pmin(ie)=minval(pstens(:,:,ie))
  enddo

  do ie=nets,nete
     pstens(:,:,ie)=pmax(ie)
     call edgeVpack(edgebuf,pstens(:,:,ie),1,0,ie)
  enddo
  call bndry_exchangeV(hybrid,edgebuf)
  do ie=nets,nete
     call edgeVunpackMax(edgebuf, pstens(:,:,ie), 1, 0, ie)
     pmax(ie)=maxval(pstens(:,:,ie))
  enddo
  enddo


  do ic=1,numcycle

     ! recompute halo each step?
     !if (minmax_halo==0) then
     !   do ie=nets,nete
     !      pmin(ie)=minval(phis(:,:,ie))
     !      pmax(ie)=maxval(phis(:,:,ie))
     !   enddo
     !endif
     if (p2filt>=1) then
        ! apply p2 filter before laplace
        do ie=nets,nete
           do i=1,np
              phis4=phis(i,:,ie)
              phis(i,2,ie)=(xgll(3)*phis4(1)+phis4(2)+phis4(3)+xgll(2)*phis4(4))/2
              phis(i,3,ie)=(xgll(2)*phis4(1)+phis4(2)+phis4(3)+xgll(3)*phis4(4))/2
           enddo
           do j=1,np
              phis4=phis(:,j,ie)
              phis(2,j,ie)=(xgll(3)*phis4(1)+phis4(2)+phis4(3)+xgll(2)*phis4(4))/2
              phis(3,j,ie)=(xgll(2)*phis4(1)+phis4(2)+phis4(3)+xgll(3)*phis4(4))/2
           enddo
        end do
     endif


     do ie=nets,nete
        pstens(:,:,ie)=laplace_sphere_wk(phis(:,:,ie),deriv,elem(ie),var_coef=.true.)
     enddo

     do ie=nets,nete
        !  ps(t+1) = ps(t) + Minv * DSS * M * RHS
        !  ps(t+1) = Minv * DSS * M [ ps(t) +  RHS ]
        ! but output of biharminc_wk is of the form M*RHS.  rewrite as:
        !  ps(t+1) = Minv * DSS * M [ ps(t) +  M*RHS/M ]
        ! so we can apply limiter to ps(t) +  (M*RHS)/M
        phis(:,:,ie)=phis(:,:,ie) + &
           smooth_phis_nudt*pstens(:,:,ie)/elem(ie)%spheremp(:,:)


        if (minmax_halo>=0) then
        ! remove new extrema.  could use conservative reconstruction from advection
        ! but no reason to conserve mean PHI.
        do i=1,np
        do j=1,np
           if (phis(i,j,ie)>pmax(ie)) phis(i,j,ie)=pmax(ie)
           if (phis(i,j,ie)<pmin(ie)) phis(i,j,ie)=pmin(ie)
        enddo
        enddo
        endif

        ! user specified minimum
        do i=1,np
        do j=1,np
           if (phis(i,j,ie)<minf) phis(i,j,ie)=minf
        enddo
        enddo

        phis(:,:,ie)=phis(:,:,ie)*elem(ie)%spheremp(:,:)
        call edgeVpack(edgebuf,phis(:,:,ie),1,0,ie)

     enddo

     call bndry_exchangeV(hybrid,edgebuf)

     do ie=nets,nete
        call edgeVunpack(edgebuf, phis(:,:,ie), 1, 0, ie)
        phis(:,:,ie)=phis(:,:,ie)*elem(ie)%rspheremp(:,:)
     enddo

  enddo

  if (p2filt==2) then
     ! apply final p2 filter 
     do ie=nets,nete
        do i=1,np
           phis4=phis(i,:,ie)
           phis(i,2,ie)=(xgll(3)*phis4(1)+phis4(2)+phis4(3)+xgll(2)*phis4(4))/2
           phis(i,3,ie)=(xgll(2)*phis4(1)+phis4(2)+phis4(3)+xgll(3)*phis4(4))/2
        enddo
        do j=1,np
           phis4=phis(:,j,ie)
           phis(2,j,ie)=(xgll(3)*phis4(1)+phis4(2)+phis4(3)+xgll(2)*phis4(4))/2
           phis(3,j,ie)=(xgll(2)*phis4(1)+phis4(2)+phis4(3)+xgll(3)*phis4(4))/2
        enddo
     end do
  endif
  



  call FreeEdgeBuffer(edgebuf) 

  end subroutine smooth_phis



#else


subroutine neighbor_minmax(elem,hybrid,edge3,nets,nete,nt,min_neigh,max_neigh,min_var,max_var,kmass)
!
! compute Q min&max over the element and all its neighbors
!
!
integer :: nets,nete,nt
type (hybrid_t)      , intent(in) :: hybrid
type (element_t)     , intent(inout) :: elem(:)
type (EdgeBuffer_t)  , intent(in) :: edge3
real (kind=real_kind) :: min_neigh(nlev,nets:nete)
real (kind=real_kind) :: max_neigh(nlev,nets:nete)
real (kind=real_kind),optional :: min_var(nlev,nets:nete)
real (kind=real_kind),optional :: max_var(nlev,nets:nete)
real (kind=real_kind) :: Qmin(np,np,nlev)
real (kind=real_kind) :: Qmax(np,np,nlev)
real (kind=real_kind) :: Qvar(np,np,nlev)
integer, optional :: kmass
type (EdgeDescriptor_t), allocatable :: desc(:)

! local
integer :: ie,k,q

  if(present(kmass))then
!the check if kmass is a valid number is done in sweq_mod
    do k=1,nlev
      if(k.ne.kmass)then
         do ie=nets,nete
            elem(ie)%state%p(:,:,k,nt)=elem(ie)%state%p(:,:,k,nt)/&
            elem(ie)%state%p(:,:,kmass,nt)
         enddo
      endif
    enddo
  endif

    ! compute p min, max
    do ie=nets,nete
#if (defined COLUMN_OPENMP)
!$omp parallel do private(k)
#endif
       do k=1,nlev
          Qmin(:,:,k)=minval(elem(ie)%state%p(:,:,k,nt))
          Qmax(:,:,k)=maxval(elem(ie)%state%p(:,:,k,nt))
          ! max - min - crude approximation to TV within the element:
          Qvar(:,:,k)=Qmax(1,1,k)-Qmin(1,1,k)
       enddo
       ! cant use more efficient edgeVpack_nlyr() interface because below we use the old
       ! edgeVunpackMin/Max interface
       call edgeVpack(edge3,Qmax,nlev,0,ie)
       call edgeVpack(edge3,Qmin,nlev,nlev,ie)
       call edgeVpack(edge3,Qvar,nlev,2*nlev,ie)
    enddo
    
    call t_startf('nmm_bexchV')
    call bndry_exchangeV(hybrid,edge3)
    call t_stopf('nmm_bexchV')
       
    do ie=nets,nete
#if (defined COLUMN_OPENMP)
!$omp parallel do private(k)
#endif
       do k=1,nlev
          Qmin(:,:,k)=minval(elem(ie)%state%p(:,:,k,nt))
          Qmax(:,:,k)=maxval(elem(ie)%state%p(:,:,k,nt))
       enddo

       ! now unpack the min
       if (present(min_var)) then
#if (defined COLUMN_OPENMP)
!$omp parallel do private(k)
#endif
          do k=1,nlev
             Qvar(:,:,k)=Qmax(1,1,k)-Qmin(1,1,k)
          enddo
! WARNING - edgeVunpackMin/Max take second argument as input/ouput
          call edgeVunpackMin(edge3,Qvar,nlev,2*nlev,ie)
#if (defined COLUMN_OPENMP)
!$omp parallel do private(k)
#endif
          do k=1,nlev
             min_var(k,ie)=minval(Qvar(:,:,k))
          enddo
       endif

       ! now unpack the max
       if (present(max_var)) then
#if (defined COLUMN_OPENMP)
!$omp parallel do private(k)
#endif
          do k=1,nlev
             Qvar(:,:,k)=Qmax(1,1,k)-Qmin(1,1,k)
          enddo
! WARNING - edgeVunpackMin/Max take second argument as input/ouput
          call edgeVunpackMax(edge3,Qvar,nlev,2*nlev,ie)
#if (defined COLUMN_OPENMP)
!$omp parallel do private(k)
#endif
          do k=1,nlev
             max_var(k,ie)=maxval(Qvar(:,:,k))
          enddo
       endif


! WARNING - edgeVunpackMin/Max take second argument as input/ouput
       call edgeVunpackMax(edge3,Qmax,nlev,0,ie)
       call edgeVunpackMin(edge3,Qmin,nlev,nlev,ie)
#if (defined COLUMN_OPENMP)
!$omp parallel do private(k)
#endif
       do k=1,nlev
          max_neigh(k,ie)=maxval(Qmax(:,:,k))
          min_neigh(k,ie)=minval(Qmin(:,:,k))
       enddo
       
    end do

#ifdef DEBUGOMP
#if (defined HORIZ_OPENMP)
!$OMP BARRIER
#endif
#endif

  if(present(kmass))then
    do k=1,nlev
       if(k.ne.kmass)then
          do ie=nets,nete
             elem(ie)%state%p(:,:,k,nt)=elem(ie)%state%p(:,:,k,nt)*&
             elem(ie)%state%p(:,:,kmass,nt)
          enddo
       endif
    enddo
  endif
end subroutine





#endif
end module viscosity_base
