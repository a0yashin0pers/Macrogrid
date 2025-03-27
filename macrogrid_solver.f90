module macrogrid_solver
   implicit none
   private

   interface
      subroutine interface_macrogrid_solver_method()
         implicit none
      end subroutine interface_macrogrid_solver_method

      subroutine interface_subgrid_solver_method(subgrid0, subgrid_size0, dx0, dy0)
         implicit none
         integer, intent(in) :: subgrid_size0
         real*8, intent(inout) :: subgrid0(subgrid_size0*subgrid_size0)
         real*8, intent(in) :: dx0, dy0
      end subroutine interface_subgrid_solver_method
   end interface

   public :: interface_macrogrid_solver_method, interface_subgrid_solver_method
   public :: run, configure_macrogrid_solver, get_results, get_omega

   public :: simple_iteration
   public :: sor, sor_fixed_omega
   public :: conjugate_residuals

   integer :: macrogrid_size_x, macrogrid_size_y, subgrid_size, max_iter
   real*8 :: omega, eps
   real*8, pointer :: macrogrid(:,:,:,:)

   real*8 :: dx, dy
   integer :: interface_size
   procedure(interface_subgrid_solver_method), pointer :: subgrid_solver_method => null()

   real*8 :: time
   integer :: iter

contains

   subroutine configure_macrogrid_solver(new_eps_interface, new_max_iter_interface, new_omega)
      implicit none
      real*8, intent(in) :: new_eps_interface, new_omega
      integer, intent(in) :: new_max_iter_interface

      eps = new_eps_interface
      max_iter = new_max_iter_interface
      omega = new_omega

   end subroutine configure_macrogrid_solver

   subroutine run(new_macrogrid, new_macrogrid_size_x, new_macrogrid_size_y, new_subgrid_size, &
      new_macrogrid_solver_method, new_subgrid_solver_method)

      implicit none
      integer, intent(in) :: new_macrogrid_size_x, new_macrogrid_size_y, new_subgrid_size
      real*8, intent(inout), target :: new_macrogrid(:,:,:,:)
      procedure(interface_macrogrid_solver_method) :: new_macrogrid_solver_method
      procedure(interface_subgrid_solver_method) :: new_subgrid_solver_method

      macrogrid_size_x = new_macrogrid_size_x
      macrogrid_size_y = new_macrogrid_size_y
      subgrid_size = new_subgrid_size

      interface_size = (macrogrid_size_x*(macrogrid_size_y-1) + macrogrid_size_y*(macrogrid_size_x-1))*(subgrid_size-2)

      dx = 1.0d0/dble(macrogrid_size_x*subgrid_size-(macrogrid_size_x-1)-1)
      dy = 1.0d0/dble(macrogrid_size_y*subgrid_size-(macrogrid_size_y-1)-1)

      macrogrid => new_macrogrid
      subgrid_solver_method => new_subgrid_solver_method
      call new_macrogrid_solver_method()

   end subroutine run

   subroutine get_results(res_time, res_iter)
      implicit none
      real*8, intent(out) :: res_time
      integer, intent(out) :: res_iter
      res_time = time
      res_iter = iter
   end subroutine get_results

   subroutine get_omega(res_omega)
      implicit none
      real*8, intent(out) :: res_omega
      res_omega = omega
   end subroutine get_omega

   subroutine compute_subgrids()
      implicit none
      integer :: iX, iY

      do iX = 1, macrogrid_size_x
         do iY = 1, macrogrid_size_y
            call subgrid_solver_method(macrogrid(iX,iY,:,:), subgrid_size, dx, dy)
         end do
      end do

   end subroutine compute_subgrids

   ! subroutine compute_subgrids_parallel()
   !    implicit none
   !    integer :: iX, iY

   !    !$OMP PARALLEL DO PRIVATE(iX, iY) COLLAPSE(2) SCHEDULE(static) &
   !    !$OMP DEFAULT(NONE) SHARED(macrogrid, macrogrid_size_x, macrogrid_size_y, subgrid_size, dx, dy)
   !    do iX = 1, macrogrid_size_x
   !       do iY = 1, macrogrid_size_y
   !          call subgrid_solver_method(macrogrid(iX,iY,:,:), subgrid_size, dx, dy)
   !       end do
   !    end do
   !    !$OMP END PARALLEL DO

   ! end subroutine compute_subgrids_parallel

   subroutine simple_iteration()
      implicit none
      real*8 :: old_vec(interface_size), new_vec(interface_size)
      real*8 :: interface_error, start_time, end_time
      integer :: iX, iY, k
      real*8 :: new_val, old_value

      call cpu_time(start_time)

      do iter = 1, max_iter

         interface_error = 0.0d0
         call compute_subgrids()

         do iX = 1, macrogrid_size_x
            do iY = 1, macrogrid_size_y - 1
               do k = 2, subgrid_size - 1

                  old_value = macrogrid(iX, iY, k, subgrid_size)

                  new_val = (4.0d0*(macrogrid(iX, iY, k, subgrid_size-1) + macrogrid(iX, iY+1, k, 2)) - &
                     macrogrid(iX, iY, k, subgrid_size-2) - macrogrid(iX, iY+1, k, 3))/6.0d0

                  macrogrid(iX, iY, k, subgrid_size) = new_val
                  macrogrid(iX, iY+1, k, 1) = new_val

                  interface_error = interface_error + dabs(new_val - old_value)

               end do
            end do
         end do

         do iX = 1, macrogrid_size_x - 1
            do iY = 1, macrogrid_size_y
               do k = 2, subgrid_size - 1

                  old_value = macrogrid(iX, iY, subgrid_size, k)

                  new_val = (4.0d0*(macrogrid(iX, iY, subgrid_size-1, k) + macrogrid(iX+1, iY, 2, k)) - &
                     macrogrid(iX,iY, subgrid_size-2, k) - macrogrid(iX+1, iY, 3, k))/6.d0

                  macrogrid(iX, iY, subgrid_size, k) = new_val
                  macrogrid(iX+1, iY, 1, k) = new_val

                  interface_error = interface_error + dabs(new_val - old_value)

               end do
            end do
         end do

         if (interface_error < eps) then
            exit
         end if

      end do

      call compute_intersection_nodes()

      call cpu_time(end_time)
      time = end_time - start_time

   end subroutine simple_iteration

   subroutine sor()
      implicit none
      real*8 :: old_vec(interface_size), new_vec(interface_size)
      real*8 :: interface_error, start_time, end_time
      integer :: iX, iY, k
      real*8 :: new_val, old_value, norm1, norm2, pz

      call cpu_time(start_time)

      omega = 1.0d0
      
      do iter = 1, max_iter

         interface_error = 0.0d0
         norm2 = 0.0d0
         call compute_subgrids()

         do iX = 1, macrogrid_size_x
            do iY = 1, macrogrid_size_y - 1
               do k = 2, subgrid_size - 1

                  old_value = macrogrid(iX, iY, k, subgrid_size)

                  new_val = (4.0d0*(macrogrid(iX, iY, k, subgrid_size-1) + macrogrid(iX, iY+1, k, 2)) - &
                     macrogrid(iX, iY, k, subgrid_size-2) - macrogrid(iX, iY+1, k, 3))/6.0d0

                  new_val = new_val * omega + (1.0d0-omega) * old_value

                  macrogrid(iX, iY, k, subgrid_size) = new_val
                  macrogrid(iX, iY+1, k, 1) = new_val

                  norm2 = norm2 + (new_val - old_value)*(new_val - old_value)
                  interface_error = interface_error + dabs(new_val - old_value)

               end do
            end do
         end do

         do iX = 1, macrogrid_size_x - 1
            do iY = 1, macrogrid_size_y
               do k = 2, subgrid_size - 1

                  old_value = macrogrid(iX, iY, subgrid_size, k)

                  new_val = (4.0d0*(macrogrid(iX, iY, subgrid_size-1, k) + macrogrid(iX+1, iY, 2, k)) - &
                     macrogrid(iX,iY, subgrid_size-2, k) - macrogrid(iX+1, iY, 3, k))/6.d0

                  new_val = new_val * omega + (1.0d0-omega) * old_value

                  macrogrid(iX, iY, subgrid_size, k) = new_val
                  macrogrid(iX+1, iY, 1, k) = new_val

                  norm2 = norm2 + (new_val - old_value)*(new_val - old_value)
                  interface_error = interface_error + dabs(new_val - old_value)

               end do
            end do
         end do

         if (iter > 1) then
            pz = dsqrt(norm2)/dsqrt(norm1)
            if (pz < 1.0d0) then
               omega = 2.0d0 / (1.0d0 + dsqrt(1.0d0 - pz))
            end if
         end if

         norm1 = norm2

         if (interface_error < eps) then
            exit
         end if

      end do

      call compute_intersection_nodes()

      call cpu_time(end_time)
      time = end_time - start_time

   end subroutine sor

   subroutine sor_fixed_omega()
      implicit none
      real*8 :: old_vec(interface_size), new_vec(interface_size)
      real*8 :: interface_error, start_time, end_time
      integer :: iX, iY, k
      real*8 :: new_val, old_value

      call cpu_time(start_time)

      do iter = 1, max_iter

         interface_error = 0.0d0
         call compute_subgrids()

         do iX = 1, macrogrid_size_x
            do iY = 1, macrogrid_size_y - 1
               do k = 2, subgrid_size - 1

                  old_value = macrogrid(iX, iY, k, subgrid_size)

                  new_val = (4.0d0*(macrogrid(iX, iY, k, subgrid_size-1) + macrogrid(iX, iY+1, k, 2)) - &
                     macrogrid(iX, iY, k, subgrid_size-2) - macrogrid(iX, iY+1, k, 3))/6.0d0

                  new_val = new_val * omega + (1.0d0-omega) * old_value

                  macrogrid(iX, iY, k, subgrid_size) = new_val
                  macrogrid(iX, iY+1, k, 1) = new_val

                  interface_error = interface_error + dabs(new_val - old_value)

               end do
            end do
         end do

         do iX = 1, macrogrid_size_x - 1
            do iY = 1, macrogrid_size_y
               do k = 2, subgrid_size - 1

                  old_value = macrogrid(iX, iY, subgrid_size, k)

                  new_val = (4.0d0*(macrogrid(iX, iY, subgrid_size-1, k) + macrogrid(iX+1, iY, 2, k)) - &
                     macrogrid(iX,iY, subgrid_size-2, k) - macrogrid(iX+1, iY, 3, k))/6.d0

                  new_val = new_val * omega + (1.0d0-omega) * old_value

                  macrogrid(iX, iY, subgrid_size, k) = new_val
                  macrogrid(iX+1, iY, 1, k) = new_val

                  interface_error = interface_error + dabs(new_val - old_value)

               end do
            end do
         end do

         if (interface_error < eps) then
            exit
         end if

      end do

      call compute_intersection_nodes()

      call cpu_time(end_time)
      time = end_time - start_time

   end subroutine sor_fixed_omega

   subroutine conjugate_residuals()
      implicit none
      real(8), dimension(interface_size) :: b, u_new, u_old, r_old, r_new, p_old, p_new
      real(8), dimension(interface_size) :: Ar_old, Ar_new, Ap_old, Ap_new
      real(8) :: alpha, beta, temp0, temp1, old_value, new_val, start_time, end_time, interface_error
      integer :: i

      call cpu_time(start_time)

      call compute_subgrids()
      call s(b)

      u_old = 0.0d0

      call Cv(u_old, b, Ar_old)
      r_new = - b - Ar_old

      p_new = r_new

      call Cv(r_new, b, Ar_new)
      temp0 = dot_product(Ar_new, r_new)
      temp1 = dot_product(Ar_new, Ar_new)
      alpha = temp0 / temp1

      r_old = r_new - alpha * Ar_new

      call Cv(r_old, b, Ar_old)
      temp0 = dot_product(Ar_old, r_old)
      temp1 = dot_product(Ar_new, r_new)
      beta = temp0 / temp1

      p_old = r_old + beta * p_new
      Ap_old = Ar_old + beta * Ar_new
      u_old = alpha * p_new

      do iter = 1, max_iter

         temp0 = dot_product(Ar_old, r_old)
         temp1 = dot_product(Ap_old, Ap_old)
         alpha = temp0 / temp1

         u_new = u_old + alpha * p_old

         interface_error = 0.0d0
         do i = 1, interface_size
            interface_error = interface_error + dabs(u_new(i) - u_old(i))
         end do
         if (interface_error < eps) then
            exit
         end if

         r_new = r_old - alpha * Ap_old

         call Cv(r_new, b, Ar_new)
         temp0 = dot_product(Ar_new, r_new)
         temp1 = dot_product(Ar_old, r_old)
         beta = temp0 / temp1

         p_old = r_new + beta * p_old
         Ap_old = Ar_new + beta * Ap_old

         r_old = r_new
         Ar_old = Ar_new
         u_old = u_new

      end do

      call cpu_time(end_time)
      time = end_time - start_time

      call set_interface(u_new)
      call compute_subgrids()
      call compute_intersection_nodes()

   end subroutine conjugate_residuals

   subroutine compute_intersection_nodes()
      implicit none
      integer :: iX, iY
      real*8 :: old_value, new_val

      do iX = 1, macrogrid_size_x - 1
         do iY = 1, macrogrid_size_y - 1

            old_value = macrogrid(iX, iY, subgrid_size, subgrid_size)

            new_val = (macrogrid(iX, iY, subgrid_size-1, subgrid_size) + &
               macrogrid(iX+1, iY, 2, subgrid_size) + &
               macrogrid(iX, iY, subgrid_size, subgrid_size-1) + &
               macrogrid(iX, iY+1, subgrid_size, 2 ))/4.d0

            macrogrid(iX, iY, subgrid_size, subgrid_size) = new_val
            macrogrid(iX+1, iY, 1, subgrid_size) = new_val
            macrogrid(iX, iY+1, subgrid_size, 1) = new_val
            macrogrid(iX+1, iY+1, 1, 1) = new_val

         end do
      end do

   end subroutine compute_intersection_nodes

   subroutine set_interface(vec)
      implicit none
      real*8, intent(in) :: vec(interface_size)
      integer :: iX, iY, k, i

      i = 1
      do iX = 1, macrogrid_size_x
         do iY = 1, macrogrid_size_y - 1
            do k = 2, subgrid_size - 1

               macrogrid(iX, iY, k, subgrid_size) = vec(i)
               macrogrid(iX, iY+1, k, 1) = vec(i)
               i = i + 1

            end do
         end do
      end do

      do iX = 1, macrogrid_size_x - 1
         do iY = 1, macrogrid_size_y
            do k = 2, subgrid_size - 1

               macrogrid(iX, iY, subgrid_size, k) = vec(i)
               macrogrid(iX+1, iY, 1, k) = vec(i)
               i = i + 1

            end do
         end do
      end do

   end subroutine set_interface

   subroutine s(s_vec)
      implicit none
      real*8,  intent(out) :: s_vec(interface_size)

      integer :: iX, iY, i1, j1, k, i

      i = 1
      do iX = 1, macrogrid_size_x
         do iY = 1, macrogrid_size_y - 1
            do k = 2, subgrid_size - 1

               s_vec(i) = (3.0d0*(macrogrid(iX, iY, k, subgrid_size) + macrogrid(iX, iY+1, k, 1)) - &
                  4.0d0*(macrogrid(iX, iY, k, subgrid_size-1 ) + macrogrid(iX, iY+1, k, 2 )) + &
                  macrogrid(iX, iY, k, subgrid_size-2) + macrogrid(iX, iY+1, k, 3 ))/(2.0d0*dy)
               i = i + 1

            end do
         end do
      end do

      do iX = 1, macrogrid_size_x - 1
         do iY = 1, macrogrid_size_y
            do k = 2, subgrid_size - 1

               s_vec(i) = (3.0d0*(macrogrid(iX, iY, subgrid_size, k) + macrogrid(iX+1, iY, 1, k)) - &
                  4.0d0*(macrogrid(iX,iY, subgrid_size-1, k) + macrogrid(iX+1, iY, 2, k)) + &
                  macrogrid(iX,iY, subgrid_size-2, k) + macrogrid(iX+1, iY, 3, k))/(2.0d0*dx)
               i = i + 1

            end do
         end do
      end do

   end subroutine s

   subroutine Cv(v, b, Cv_vec)
      implicit none
      real*8, intent(in) :: v(interface_size), b(interface_size)
      real*8, intent(out) :: Cv_vec(interface_size)

      real*8 :: subgrid_error

      call set_interface(v)
      call compute_subgrids()
      call s(Cv_vec)

      Cv_vec = Cv_vec - b

   end subroutine Cv

end module macrogrid_solver
