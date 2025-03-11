module macrogrid_solvers
   implicit none

   public :: simple_iteration
   public :: conjugate_residuals

contains

   subroutine simple_iteration(subgrid_solver, macrogrid,   &
      macrogrid_size_x, macrogrid_size_y,                   &
      subgrid_size, omega,                                  &
      eps_subgrid, eps_interface,                           &
      max_iter_subgrid, max_iter_interface,                 &
      interface_error, time, iter)
      implicit none
      external :: subgrid_solver
      integer, intent(in) :: macrogrid_size_x, macrogrid_size_y, subgrid_size, max_iter_subgrid, max_iter_interface
      real*8,  intent(inout) :: macrogrid(macrogrid_size_x,macrogrid_size_y,subgrid_size,subgrid_size)
      real*8,  intent(in)   :: omega, eps_subgrid, eps_interface
      real*8,  intent(out)   :: interface_error, time
      integer,  intent(out)   :: iter

      integer :: iX, iY, k
      real*8  :: dx, dy, old_value, new_val, subgrid_error, start_time, end_time

      call compute_dxdy(macrogrid_size_x, macrogrid_size_y, subgrid_size, dx, dy)
      call cpu_time(start_time)

      do iter = 1, max_iter_interface

         interface_error = 0.0d0

         do iX = 1, macrogrid_size_x
            do iY = 1, macrogrid_size_y
               call subgrid_solver(macrogrid(iX,iY,:,:), subgrid_size, &
                  omega, eps_subgrid, max_iter_subgrid, subgrid_error, dx, dy)
            end do
         end do

         do iX = 1, macrogrid_size_x
            do iY = 1, macrogrid_size_y - 1
               do k = 2, subgrid_size - 1

                  old_value = macrogrid(iX, iY, k, subgrid_size)

                  new_val = (-macrogrid(iX, iY,   k, subgrid_size-2 ) - macrogrid(iX, iY+1,   k, 3 )  + &
                     4.0d0*macrogrid(iX, iY,   k,   subgrid_size-1 ) +4.0d0*macrogrid(iX, iY+1, k, 2 ))/6.0d0

                  macrogrid(iX, iY,   k, subgrid_size) = new_val
                  macrogrid(iX, iY+1, k, 1)            = new_val

                  interface_error = interface_error + dabs(new_val - old_value)
               end do
            end do
         end do

         do iX = 1, macrogrid_size_x - 1
            do iY = 1, macrogrid_size_y
               do k = 2, subgrid_size - 1

                  old_value = macrogrid(iX, iY, subgrid_size, k)

                  new_val = (-macrogrid(iX,iY, subgrid_size-2, k) + 4.0d0*macrogrid(iX+1, iY, 2, k ) + &
                     4.0d0*macrogrid(iX,   iY, subgrid_size-1, k) - macrogrid(iX+1, iY, 3, k) &
                     ) / 6.d0

                  macrogrid(iX,   iY, subgrid_size, k) = new_val
                  macrogrid(iX+1, iY, 1,            k) = new_val

                  interface_error = interface_error + dabs(new_val - old_value)
               end do
            end do
         end do

         do iX = 1, macrogrid_size_x - 1
            do iY = 1, macrogrid_size_y - 1

               old_value = macrogrid(iX, iY, subgrid_size, subgrid_size)
               new_val = ( &
                  macrogrid(iX, iY, subgrid_size-1, subgrid_size) +  &
                  macrogrid(iX+1, iY, 2, subgrid_size)  + &
                  macrogrid(iX, iY, subgrid_size, subgrid_size-1) + &
                  macrogrid(iX, iY+1, subgrid_size, 2 )  &
                  ) / 4.d0

               macrogrid(iX,   iY,   subgrid_size,   subgrid_size)   = new_val
               macrogrid(iX+1, iY,   1,              subgrid_size)   = new_val
               macrogrid(iX,   iY+1, subgrid_size,   1)              = new_val
               macrogrid(iX+1, iY+1, 1,              1)              = new_val

               interface_error = interface_error + dabs(new_val - old_value)
            end do
         end do

         if (interface_error < eps_interface) then
            exit
         end if

      end do

      call cpu_time(end_time)
      time = end_time - start_time

   end subroutine simple_iteration

   subroutine conjugate_residuals(subgrid_solver, macrogrid,&
      macrogrid_size_x, macrogrid_size_y,                   &
      subgrid_size, omega,                                  &
      eps_subgrid, eps_interface,                           &
      max_iter_subgrid, max_iter_interface,                 &
      interface_error, time, iter)
      implicit none
      external :: subgrid_solver
      integer, intent(in) :: macrogrid_size_x, macrogrid_size_y, subgrid_size, max_iter_subgrid, max_iter_interface
      real*8,  intent(inout) :: macrogrid(macrogrid_size_x,macrogrid_size_y,subgrid_size,subgrid_size)
      real*8,  intent(in)   :: omega, eps_subgrid, eps_interface
      real*8,  intent(out)   :: interface_error, time
      integer,  intent(out)   :: iter

      integer :: interface_size, i, iX, iY, k
      real(8), dimension((macrogrid_size_x*(macrogrid_size_y-1) + macrogrid_size_y*(macrogrid_size_x-1))*(subgrid_size-2)) :: &
         b, u_new, u_old, r_old, r_new, p_old, p_new, Ar_old, Ar_new, Ap_old, Ap_new

      real(8) :: alpha, beta, dx, dy, temp0, temp1, subgrid_error, old_value, new_val, start_time, end_time

      interface_size = (macrogrid_size_x*(macrogrid_size_y-1) + macrogrid_size_y*(macrogrid_size_x-1))*(subgrid_size-2)
      call compute_dxdy(macrogrid_size_x, macrogrid_size_y, subgrid_size, dx, dy)

      call cpu_time(start_time)

      do iX = 1, macrogrid_size_x
         do iY = 1, macrogrid_size_y
            call subgrid_solver(macrogrid(iX,iY,:,:), subgrid_size, &
               omega, eps_subgrid, max_iter_subgrid, subgrid_error, dx, dy)
         end do
      end do
      call s(macrogrid, dx, dy, macrogrid_size_x, macrogrid_size_y, subgrid_size, interface_size, b)

      u_old = 0.0d0

      call Cv(subgrid_solver, macrogrid, u_old, b, dx, dy, macrogrid_size_x, macrogrid_size_y, subgrid_size, interface_size, omega, eps_subgrid, max_iter_subgrid, Ar_old)
      r_new = - b - Ar_old

      p_new = r_new

      call Cv(subgrid_solver, macrogrid, r_new, b, dx, dy, macrogrid_size_x, macrogrid_size_y, subgrid_size, interface_size, omega, eps_subgrid, max_iter_subgrid, Ar_new)
      call dot_product(Ar_new, r_new, interface_size, temp0)
      call dot_product(Ar_new, Ar_new, interface_size, temp1)
      alpha = temp0 / temp1

      r_old = r_new - alpha * Ar_new

      call Cv(subgrid_solver, macrogrid, r_old, b, dx, dy, macrogrid_size_x, macrogrid_size_y, subgrid_size, interface_size, omega, eps_subgrid, max_iter_subgrid, Ar_old)
      call dot_product(Ar_old, r_old, interface_size, temp0)
      call dot_product(Ar_new, r_new, interface_size, temp1)
      beta = temp0 / temp1

      p_old = r_old + beta * p_new
      Ap_old = Ar_old + beta * Ar_new
      u_old = alpha * p_new

      do iter = 1, max_iter_interface

         call dot_product(Ar_old, r_old, interface_size, temp0)
         call dot_product(Ap_old, Ap_old, interface_size, temp1)
         alpha = temp0 / temp1

         u_new = u_old + alpha * p_old

         interface_error = 0.0d0
         do i = 1, interface_size
            interface_error = interface_error + dabs(u_new(i) - u_old(i))
         end do

         if (interface_error < eps_interface) then
            exit
         end if

         r_new = r_old - alpha * Ap_old

         call Cv(subgrid_solver, macrogrid, r_new, b, dx, dy, macrogrid_size_x, macrogrid_size_y, subgrid_size, interface_size, omega, eps_subgrid, max_iter_subgrid, Ar_new)
         call dot_product(Ar_new, r_new, interface_size, temp0)
         call dot_product(Ar_old, r_old, interface_size, temp1)
         beta = temp0 / temp1

         p_old = r_new + beta * p_old
         Ap_old = Ar_new + beta * Ap_old

         r_old = r_new
         Ar_old = Ar_new
         u_old = u_new

      end do

      call cpu_time(end_time)
      time = end_time - start_time

      call set_interface(macrogrid, u_new, macrogrid_size_x, macrogrid_size_y, subgrid_size, interface_size)
      do iX = 1, macrogrid_size_x
         do iY = 1, macrogrid_size_y
            call subgrid_solver(macrogrid(iX,iY,:,:), subgrid_size, &
               omega, eps_subgrid, max_iter_subgrid, subgrid_error, dx, dy)
         end do
      end do

      do iX = 1, macrogrid_size_x - 1
         do iY = 1, macrogrid_size_y - 1

            old_value = macrogrid(iX, iY, subgrid_size, subgrid_size)
            new_val = ( &
               macrogrid(iX, iY, subgrid_size-1, subgrid_size) +  &
               macrogrid(iX+1, iY, 2, subgrid_size)  + &
               macrogrid(iX, iY, subgrid_size, subgrid_size-1) + &
               macrogrid(iX, iY+1, subgrid_size, 2 )  &
               ) / 4.d0

            macrogrid(iX,   iY,   subgrid_size,   subgrid_size)   = new_val
            macrogrid(iX+1, iY,   1,              subgrid_size)   = new_val
            macrogrid(iX,   iY+1, subgrid_size,   1)              = new_val
            macrogrid(iX+1, iY+1, 1,              1)              = new_val
         end do
      end do

   end subroutine conjugate_residuals

   subroutine compute_dxdy(macrogrid_size_x, macrogrid_size_y, subgrid_size, dx, dy)
      implicit none
      integer, intent(in) :: macrogrid_size_x, macrogrid_size_y, subgrid_size
      real*8,  intent(out) :: dx, dy
      integer :: global_size_x
      integer :: global_size_y
      global_size_x = macrogrid_size_x*subgrid_size - (macrogrid_size_x - 1)
      global_size_y = macrogrid_size_y*subgrid_size - (macrogrid_size_y - 1)
      dx = 1.0d0 / dble(global_size_x - 1)
      dy = 1.0d0 / dble(global_size_y - 1)
   end subroutine compute_dxdy

   subroutine get_interface(macrogrid, macrogrid_size_x, macrogrid_size_y, subgrid_size, interface_size, u)
      implicit none
      integer, intent(in) :: macrogrid_size_x, macrogrid_size_y, subgrid_size, interface_size
      real*8,  intent(in) :: macrogrid(macrogrid_size_x,macrogrid_size_y,subgrid_size,subgrid_size)
      real*8,  intent(out) :: u(interface_size)

      integer :: iX, iY, i1, j1, k, i

      i = 1
      do iX = 1, macrogrid_size_x
         do iY = 1, macrogrid_size_y - 1
            do k = 2, subgrid_size - 1

               u(i) = macrogrid(iX, iY, k, subgrid_size)
               i = i + 1

            end do
         end do
      end do

      do iX = 1, macrogrid_size_x - 1
         do iY = 1, macrogrid_size_y
            do k = 2, subgrid_size - 1

               u(i) = macrogrid(iX, iY, subgrid_size, k)
               i = i + 1

            end do
         end do
      end do

   end subroutine get_interface

   subroutine set_interface(macrogrid, u, macrogrid_size_x, macrogrid_size_y, subgrid_size, interface_size)
      implicit none
      integer, intent(in) :: macrogrid_size_x, macrogrid_size_y, subgrid_size, interface_size
      real*8,  intent(in) :: u(interface_size)
      real*8,  intent(inout) :: macrogrid(macrogrid_size_x,macrogrid_size_y,subgrid_size,subgrid_size)

      integer :: iX, iY, i1, j1, k, i

      i = 1
      do iX = 1, macrogrid_size_x
         do iY = 1, macrogrid_size_y - 1
            do k = 2, subgrid_size - 1

               macrogrid(iX, iY, k, subgrid_size) = u(i)
               macrogrid(iX, iY+1, k, 1) = u(i)
               i = i + 1

            end do
         end do
      end do

      do iX = 1, macrogrid_size_x - 1
         do iY = 1, macrogrid_size_y
            do k = 2, subgrid_size - 1

               macrogrid(iX, iY, subgrid_size, k) = u(i)
               macrogrid(iX+1, iY, 1, k) = u(i)
               i = i + 1

            end do
         end do
      end do

   end subroutine set_interface

   subroutine s(macrogrid, dx, dy, macrogrid_size_x, macrogrid_size_y, subgrid_size, interface_size, s_vec)
      implicit none
      integer, intent(in) :: macrogrid_size_x, macrogrid_size_y, subgrid_size, interface_size
      real*8,  intent(in) :: dx, dy
      real*8,  intent(inout) :: macrogrid(macrogrid_size_x,macrogrid_size_y,subgrid_size,subgrid_size)
      real*8,  intent(out) :: s_vec(interface_size)

      integer :: iX, iY, i1, j1, k, i

      i = 1
      do iX = 1, macrogrid_size_x
         do iY = 1, macrogrid_size_y - 1
            do k = 2, subgrid_size - 1

               s_vec(i) = (3.0d0*macrogrid(iX, iY, k, subgrid_size) - 4.0d0*macrogrid(iX, iY, k, subgrid_size-1 ) + macrogrid(iX, iY, k, subgrid_size-2) - &
                  (- 3.0d0*macrogrid(iX, iY+1, k, 1) + 4.0d0*macrogrid(iX, iY+1, k, 2 ) - macrogrid(iX, iY+1, k, 3 ))) / (2.0d0 * dy)
               i = i + 1

            end do
         end do
      end do

      do iX = 1, macrogrid_size_x - 1
         do iY = 1, macrogrid_size_y
            do k = 2, subgrid_size - 1

               s_vec(i) = (3.0d0*macrogrid(iX, iY, subgrid_size, k) - 4.0d0*macrogrid(iX,iY, subgrid_size-1, k) + macrogrid(iX,iY, subgrid_size-2, k) - &
                  (- 3.0d0*macrogrid(iX+1, iY, 1, k) + 4.0d0*macrogrid(iX+1, iY, 2, k) - macrogrid(iX+1, iY, 3, k))) / (2.0d0 * dx)
               i = i + 1

            end do
         end do
      end do

   end subroutine s

   subroutine Cv(subgrid_solver, macrogrid, v, b, dx, dy, macrogrid_size_x, macrogrid_size_y, subgrid_size, interface_size, omega, eps_subgrid, max_iter_subgrid, Cv_vec)
      implicit none
      external :: subgrid_solver
      integer, intent(in) :: macrogrid_size_x, macrogrid_size_y, subgrid_size, interface_size, max_iter_subgrid
      real*8,  intent(in) :: v(interface_size), b(interface_size), dx, dy, omega, eps_subgrid
      real*8,  intent(inout) :: macrogrid(macrogrid_size_x,macrogrid_size_y,subgrid_size,subgrid_size)
      real*8,  intent(out) :: Cv_vec(interface_size)

      real*8 :: subgrid_error
      integer :: iX, iY, i1, j1, k, i

      call set_interface(macrogrid, v, macrogrid_size_x, macrogrid_size_y, subgrid_size, interface_size)

      do iX = 1, macrogrid_size_x
         do iY = 1, macrogrid_size_y
            call subgrid_solver(macrogrid(iX,iY,:,:), subgrid_size, &
               omega, eps_subgrid, max_iter_subgrid, subgrid_error, dx, dy)
         end do
      end do

      call s(macrogrid, dx, dy, macrogrid_size_x, macrogrid_size_y, subgrid_size, interface_size, Cv_vec)

      i = 1
      do iX = 1, macrogrid_size_x
         do iY = 1, macrogrid_size_y - 1
            do k = 2, subgrid_size - 1

               Cv_vec(i) = Cv_vec(i) - b(i)
               i = i + 1

            end do
         end do
      end do

      do iX = 1, macrogrid_size_x - 1
         do iY = 1, macrogrid_size_y
            do k = 2, subgrid_size - 1

               Cv_vec(i) = Cv_vec(i) - b(i)
               i = i + 1

            end do
         end do
      end do

   end subroutine Cv

   subroutine dot_product(a, b, size, result)
      implicit none
      integer, intent(in) :: size
      real*8, intent(in) :: a(size), b(size)
      real*8, intent(out) :: result
      integer :: i

      result = 0.0d0
      do i = 1, size
         result = result + a(i) * b(i)
      end do

   end subroutine dot_product

end module macrogrid_solvers
