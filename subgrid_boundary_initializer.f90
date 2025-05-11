module subgrid_boundary_initializer
   implicit none
   private

   public :: initialize_subgrid_boundary
   public :: compute_subgrid_boundary_error

   public :: get_max_geom_progress_nummber
   public :: initialize_zero_subgrid
   public :: initialize_geom_progress_subgrid_boundary
   public :: compute_geom_progress_subgrid_boundary_error

   real*8, parameter :: R1 = 0.1d0, R2 = 1.0d0
   real*8, parameter :: x_min = 0.3d0, y_min = 0.0d0
   real*8, parameter :: len = 0.4d0

   real*8, parameter :: eps = 1.0d-6
   real*8, parameter :: q = 0.5d0

contains

   subroutine initialize_subgrid_boundary(u, u_size)
      implicit none
      integer, intent(in) :: u_size
      real*8,  intent(inout) :: u(u_size*u_size)

      integer :: l0, l1, i

      do l1 = 1, u_size
         do l0 = 1, u_size

            i = l0 + u_size*(l1 - 1)

            if (l0.eq.1 .or. l0.eq.u_size .or. l1.eq.1 .or. l1.eq.u_size) then
               u(i) = log(sqrt((x_min + len*(l1 - 1)/u_size)**2 + &
                  (y_min + len*(l0 - 1)/u_size)**2)*R2/(R1*R1))/log(R2/R1)
            else
               u(i) = 0.0d0
            end if

         end do
      end do

   end subroutine initialize_subgrid_boundary

   subroutine compute_subgrid_boundary_error(u, u_size, error)
      implicit none
      integer, intent(in) :: u_size
      real*8,  intent(inout) :: u(u_size*u_size)
      real*8,  intent(out) :: error

      integer :: l0, l1, i

      error = 0.0d0

      do l1 = 1, u_size
         do l0 = 1, u_size

            i = l0 + u_size*(l1 - 1)
            error = max(error, abs(u(i) - log(sqrt((x_min + len*(l1 - 1)/u_size)**2 + &
               (y_min + len*(l0 - 1)/u_size)**2)*R2/(R1*R1))/log(R2/R1)))

         end do
      end do

   end subroutine compute_subgrid_boundary_error

   subroutine get_max_geom_progress_nummber(k)
      implicit none
      integer,  intent(out) :: k

      k = ceiling(log(eps) / log(q))

   end subroutine get_max_geom_progress_nummber

   subroutine initialize_zero_subgrid(u, u_size)
      implicit none
      integer, intent(in) :: u_size
      real*8,  intent(inout) :: u(u_size*u_size)

      integer :: l0, l1, i

      do l1 = 1, u_size
         do l0 = 1, u_size
            i = l0 + u_size*(l1 - 1)
            u(i) = 0.0d0
         end do
      end do

   end subroutine initialize_zero_subgrid

   subroutine initialize_geom_progress_subgrid_boundary(u, u_size, k)
      implicit none
      integer, intent(in) :: u_size, k
      real*8,  intent(inout) :: u(u_size*u_size)

      integer :: l0, l1, i
      real*8 :: s

      s = 1.0d0 - q**k

      do l1 = 1, u_size
         do l0 = 1, u_size
            i = l0 + u_size*(l1 - 1)
            if (l0.eq.1 .or. l0.eq.u_size .or. l1.eq.1 .or. l1.eq.u_size) then
               u(i) = s
            end if
         end do
      end do

   end subroutine initialize_geom_progress_subgrid_boundary

   subroutine compute_geom_progress_subgrid_boundary_error(u, u_size, k, error)
      implicit none
      integer, intent(in) :: u_size, k
      real*8,  intent(inout) :: u(u_size*u_size)
      real*8,  intent(out) :: error

      integer :: l0, l1, i
      real*8 :: s

      s = 1.0d0 - q**k
      error = 0.0d0

      do l1 = 1, u_size
         do l0 = 1, u_size

            i = l0 + u_size*(l1 - 1)
            error = max(error, abs(u(i) - s))

         end do
      end do

   end subroutine compute_geom_progress_subgrid_boundary_error

end module subgrid_boundary_initializer
