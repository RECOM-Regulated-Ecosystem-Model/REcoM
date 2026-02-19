! ========================================================================
! Halo exchange routines + broadcast routines that collect information
! on the entire field (needed for output)
! The routines here are very similar, difference is the data type and
! exchange pattern.
! exchange_nod2D_i(arr(myDim_nod2D+eDim_nod2D))    INTEGER
! exchange_nod2D(arr(myDim_nod2D+eDim_nod2D))      WP
! exchange_nod3D(arr(nl-1,myDim_nod2D+eDim_nod2D)) WP
! exchange_nod3D_full(arr(nl,myDim_nod2D+eDim_nod2D)) WP
! exchange_edge2D(edge_array2D)     WP  not used currently  !!! no buffer!!!
! exchange_edge3D(edge_array3D)     WP  not used currently  !!! no buffer!!!
! exchange_elem3D(elem_array3D)     WP
! exchange_elem2d_full
! exchange_elem2d_full_i
! ========================================================================

module recom_g_comm

  use, intrinsic :: ISO_FORTRAN_ENV, only: int16, int32, real32, real64

  implicit none


contains

! ========================================================================
! General version of the communication routine for 2D nodal fields
subroutine recom_exchange_nod2D(nod_array2D, npes, sn, rn, MPI_COMM_FESOM, mype,            &
                                s_mpitype_nod2D, r_mpitype_nod2D, sPE, rPE, requests, nreq, &
                                luse_g2g)

IMPLICIT NONE

logical, intent(in),optional                  :: luse_g2g
integer, intent(in)                           :: sn, rn, npes, MPI_COMM_FESOM, mype
integer, intent(inout)                        :: nreq
integer, intent(in),    dimension(:)          :: sPE, rPE
integer, intent(inout), dimension(:)          :: requests
integer, intent(in),    dimension(:), pointer :: s_mpitype_nod2D, r_mpitype_nod2D
real(real64),   intent(inout)                 :: nod_array2D(:)

 if (npes > 1) then
    call recom_exchange_nod2D_begin(nod_array2D, npes, sn, rn, MPI_COMM_FESOM, mype, &
                                    s_mpitype_nod2D, r_mpitype_nod2D, sPE, rPE,      &
                                    requests, nreq, luse_g2g)
    call recom_exchange_nod_end(npes, nreq, requests)
 end if

END SUBROUTINE recom_exchange_nod2D

! ========================================================================
! General version of the communication routine for 2D nodal fields
subroutine recom_exchange_nod2D_begin(nod_array2D, npes, sn, rn, MPI_COMM_FESOM, mype,    &
                                      s_mpitype_nod2D, r_mpitype_nod2D, sPE, rPE, requests, nreq, &
                                      luse_g2g)
IMPLICIT NONE

logical, intent(in), optional                 :: luse_g2g
integer, intent(in)                           :: sn, rn, npes, MPI_COMM_FESOM, mype
integer, intent(inout)                        :: nreq
integer, intent(in),    dimension(:)          :: sPE, rPE
integer, intent(inout), dimension(:)          :: requests
integer, intent(in),    dimension(:), pointer :: s_mpitype_nod2D, r_mpitype_nod2D
real(real64), intent(inout)                   :: nod_array2D(:)

integer                               :: n, MPIerr
logical                               :: lg2g

if(present(luse_g2g)) then
   lg2g = luse_g2g
else
   lg2g = .false.
end if

  if (npes > 1) then
     !$ACC HOST_DATA USE_DEVICE(nod_array2D) IF(lg2g)

     DO n=1,rn
        call MPI_IRECV(nod_array2D, 1, r_mpitype_nod2D(n), rPE(n), &
             rPE(n), MPI_COMM_FESOM, requests(n), MPIerr)
     END DO
     DO n=1, sn
        call MPI_ISEND(nod_array2D, 1, s_mpitype_nod2D(n), sPE(n), &
             mype, MPI_COMM_FESOM, requests(rn+n), MPIerr)
     END DO

    !$ACC END HOST_DATA
     nreq = rn+sn

  end if

END SUBROUTINE recom_exchange_nod2D_begin

! ========================================================================
! General version of the communication routine for 3D nodal fields
! stored in (vertical, horizontal) format
subroutine recom_exchange_nod3D(nod_array3D, npes, sn, rn, MPI_COMM_FESOM, mype,    &
                                s_mpitype_nod3D, r_mpitype_nod3D, sPE, rPE, requests, nreq, &
                                luse_g2g)
IMPLICIT NONE

logical, intent(in), optional                       :: luse_g2g
integer, intent(in)                                 :: sn, rn, npes, MPI_COMM_FESOM, mype
integer, intent(inout)                              :: nreq
integer, intent(in),    dimension(:)                :: sPE, rPE
integer, intent(inout), dimension(:)                :: requests
integer, intent(in),    dimension(:, :, :), pointer :: s_mpitype_nod3D, r_mpitype_nod3D
real(real64), intent(inout)                         :: nod_array3D(:,:)

if (npes > 1) then
   call recom_exchange_nod3D_begin(nod_array3D, npes, sn, rn, MPI_COMM_FESOM, mype,       &
                                   s_mpitype_nod3D, r_mpitype_nod3D, sPE, rPE,  requests, &
                                   nreq, luse_g2g)
   call recom_exchange_nod_end(npes, nreq, requests)
endif

END SUBROUTINE recom_exchange_nod3D

! ========================================================================
! General version of the communication routine for 3D nodal fields
! stored in (vertical, horizontal) format
subroutine recom_exchange_nod3D_begin(nod_array3D, npes, sn, rn, MPI_COMM_FESOM, mype,    &
                                      s_mpitype_nod3D, r_mpitype_nod3D, sPE, rPE, requests, nreq, &
                                      luse_g2g)
IMPLICIT NONE

logical, intent(in), optional                       :: luse_g2g
integer, intent(in)                                 :: sn, rn, npes, MPI_COMM_FESOM, mype
integer, intent(inout)                              :: nreq
integer, intent(in),    dimension(:)                :: sPE, rPE
integer, intent(inout), dimension(:)                :: requests
integer, intent(in),    dimension(:, :, :), pointer :: s_mpitype_nod3D, r_mpitype_nod3D
real(real64), intent(inout)                         :: nod_array3D(:,:)

integer                               :: n, MPIerr
integer                               :: nz, nl1
logical                               :: lg2g

if(present(luse_g2g)) then
   lg2g = luse_g2g
else
   lg2g = .false.
end if

 if (npes > 1) then

    nl1=ubound(nod_array3D,1)

    if ((nl1<ubound(r_mpitype_nod3D, 2)-1) .or. (nl1>ubound(r_mpitype_nod3D, 2))) then
       if (mype==0) then
          print *,'Subroutine recom_exchange_nod3D not implemented for',nl1,'layers.'
          print *,'Adding the MPI datatypes is easy, see oce_modules.F90.'
       endif
       call MPI_Abort(MPI_COMM_FESOM, 1)
    endif

    !$ACC HOST_DATA USE_DEVICE(nod_array3D) IF(lg2g)

    DO n=1,rn
       call MPI_IRECV(nod_array3D, 1, r_mpitype_nod3D(n,nl1,1), rPE(n), &
                      rPE(n), MPI_COMM_FESOM, requests(n), MPIerr)
    END DO
    DO n=1, sn
       call MPI_ISEND(nod_array3D, 1, s_mpitype_nod3D(n,nl1,1), sPE(n), &
            mype, MPI_COMM_FESOM, requests(rn+n), MPIerr)
    END DO

    !$ACC END HOST_DATA
    nreq = rn+sn

 endif
END SUBROUTINE recom_exchange_nod3D_begin

!=======================================
! AND WAITING
!=======================================

SUBROUTINE recom_exchange_nod_end(npes, request_count, array_of_requests)
use mpi
IMPLICIT NONE

integer, intent(in) :: request_count, npes
integer, intent(inout), dimension(:) :: array_of_requests

integer :: MPIerr

if (npes > 1) &
  call MPI_WAITALL(request_count, array_of_requests, MPI_STATUSES_IGNORE, MPIerr)

END SUBROUTINE recom_exchange_nod_end

end module recom_g_comm

module recom_g_comm_auto
use recom_g_comm
implicit none
interface recom_exchange_nod
      module procedure recom_exchange_nod2D
      module procedure recom_exchange_nod3D
end interface recom_exchange_nod

interface recom_exchange_nod_begin
      module procedure recom_exchange_nod2D_begin
      module procedure recom_exchange_nod3D_begin
end interface recom_exchange_nod_begin

interface recom_exchange_nod_end
  module procedure recom_exchange_nod_end
end interface

!!$interface exchange_edge
!!$      module procedure exchange_edge2D
!!$!      module procedure exchange_edge3D  ! not available, not used
!!$end interface exchange_edge

private  ! hides items not listed on public statement
public :: recom_exchange_nod, recom_exchange_nod_begin, recom_exchange_nod_end
end module recom_g_comm_auto
