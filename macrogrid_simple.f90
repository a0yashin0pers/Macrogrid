program macrogrid_solver
   implicit none

   logical :: accumulate_subgrid_error = .false.

   integer, parameter :: SUBGRID_SIZES_COUNT = 8
   integer, dimension(SUBGRID_SIZES_COUNT), parameter :: subgrid_sizes = &
      [10, 18, 34, 66, 130, 258, 514, 1026]

   real*8, dimension(SUBGRID_SIZES_COUNT), parameter :: factors_const = &
      [1.52729d0, 1.69504d0, 1.82964d0, 1.90932d0, 1.95305d0, 1.97616d0, 1.98798d0, 1.99394d0]

   real*8, dimension(SUBGRID_SIZES_COUNT), parameter :: factors_log = &
      [1.52129d0, 1.69524d0, 1.82919d0, 1.90899d0, 1.95307d0, 1.97618d0, 1.98793d0, 1.99403d0]

   integer, parameter :: SUBGRID_CFGS_COUNT = 5
   integer, dimension(SUBGRID_CFGS_COUNT), parameter :: subgrid_configs = &
      [10, 18, 34, 66, 130]

   integer, parameter :: MACRIGRID_CFGS_COUNT = 3
   integer, dimension(2, MACRIGRID_CFGS_COUNT), parameter :: macrogrid_configs = reshape([ &
      1, 1, &
      2, 2, &
      4, 4  &
      ], shape=[2, MACRIGRID_CFGS_COUNT])

   integer, parameter :: max_iter_subgrid   = 100000
   integer, parameter :: max_iter_interface = 100000
   real*8,  parameter :: eps_subgrid        = 1.0d-5
   real*8,  parameter :: eps_interface      = 1.0d-5

   real*8 :: macrogrid(10, 10, 1026, 1026)

   integer :: i_cfg, i_size, sub_sz
   real*8  :: omega
   integer :: cfg_x, cfg_y

   print *, "-----------------------------------------------"
   print *, "accumulate_subgrid_error = ", accumulate_subgrid_error
   print *, "-----------------------------------------------"

   do i_cfg = 1, MACRIGRID_CFGS_COUNT
      cfg_x = macrogrid_configs(1, i_cfg)
      cfg_y = macrogrid_configs(2, i_cfg)

      print *, "=== Macrogrid size: ", cfg_x, " x ", cfg_y, " ==="

      do i_size = 1, SUBGRID_CFGS_COUNT
         sub_sz = subgrid_configs(i_size)
         print *, "  Subgrid size =", sub_sz

         omega = factors_log(i_size)
         call initialize_logarithmic_boundary(macrogrid, cfg_x, cfg_y, sub_sz)
         call measure_execution_time(accumulate_subgrid_error,             &
            macrogrid, cfg_x, cfg_y, sub_sz,      &
            omega,                                &
            max_iter_subgrid, max_iter_interface, &
            eps_subgrid, eps_interface)
         call compute_logarithmic_boundary_error(macrogrid, cfg_x, cfg_y, sub_sz)
      end do
   end do
end program macrogrid_solver

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

subroutine initialize_constant_boundary(macrogrid, macrogrid_size_x, macrogrid_size_y, subgrid_size)
   implicit none
   integer, intent(in) :: macrogrid_size_x, macrogrid_size_y, subgrid_size
   real*8,  intent(inout) :: macrogrid(macrogrid_size_x,macrogrid_size_y,subgrid_size,subgrid_size)

   integer :: global_size_x, global_size_y
   integer :: iX, iY, i1, j1, lX, lY

   global_size_x = macrogrid_size_x * subgrid_size - (macrogrid_size_x - 1)
   global_size_y = macrogrid_size_y * subgrid_size - (macrogrid_size_y - 1)

   do iX = 1, macrogrid_size_x
      do iY = 1, macrogrid_size_y
         do i1 = 1, subgrid_size
            do j1 = 1, subgrid_size
               lX = i1 + (iX-1)*subgrid_size - (iX-1)
               lY = j1 + (iY-1)*subgrid_size - (iY-1)

               if (lX == 1 .or. lX == global_size_x .or. &
                  lY == 1 .or. lY == global_size_y) then
                  macrogrid(iX, iY, i1, j1) = 1.0d0
               else
                  macrogrid(iX, iY, i1, j1) = 0.0d0
               end if
            end do
         end do
      end do
   end do
end subroutine initialize_constant_boundary

subroutine initialize_logarithmic_boundary(macrogrid, macrogrid_size_x, macrogrid_size_y, subgrid_size)
   implicit none
   real*8, parameter :: R1 = 0.1d0, R2 = 1.0d0
   real*8, parameter :: x_min = 0.3d0, y_min = 0.0d0
   real*8, parameter :: len = 0.4d0

   integer, intent(in) :: macrogrid_size_x, macrogrid_size_y, subgrid_size
   real*8,  intent(inout) :: macrogrid(macrogrid_size_x,macrogrid_size_y,subgrid_size,subgrid_size)

   integer :: global_size_x, global_size_y
   integer :: iX, iY, i1, j1, lX, lY
   real*8  :: boundary_val

   global_size_x = macrogrid_size_x * subgrid_size - (macrogrid_size_x - 1)
   global_size_y = macrogrid_size_y * subgrid_size - (macrogrid_size_y - 1)

   do iX = 1, macrogrid_size_x
      do iY = 1, macrogrid_size_y
         do i1 = 1, subgrid_size
            do j1 = 1, subgrid_size
               lX = i1 + (iX-1)*subgrid_size - (iX-1)
               lY = j1 + (iY-1)*subgrid_size - (iY-1)

               if (lX == 1 .or. lX == global_size_x .or. &
                  lY == 1 .or. lY == global_size_y) then
                  boundary_val = log( sqrt( (x_min + len*(lX-1)/(global_size_x-1))**2 + &
                     (y_min + len*(lY-1)/(global_size_y-1))**2 )*R2/(R1*R1) ) / &
                     log(R2/R1)

                  macrogrid(iX, iY, i1, j1) = boundary_val
               else
                  macrogrid(iX, iY, i1, j1) = 0.0d0
               end if
            end do
         end do
      end do
   end do
end subroutine initialize_logarithmic_boundary

subroutine solve_sor_subdomain(u, n, factor, eps, max_iter, iter_error, dx, dy)
   implicit none
   integer, intent(in) :: n, max_iter
   real*8,  intent(in) :: factor, eps, dx, dy
   real*8,  intent(inout) :: u(n*n)
   real*8,  intent(out)   :: iter_error

   integer :: iter, l0, l1, i
   real*8  :: sum_diff, u_old, denom, tmp
   real*8  :: invdx2, invdy2

   invdx2 = 1/(dx*dx)
   invdy2 = 1/(dy*dy)
   denom  = 2.0d0*invdx2 + 2.0d0*invdy2

   do iter = 1, max_iter
      sum_diff = 0.0d0
      i = n + 2

      do l1 = 3, n
         do l0 = 3, n
            u_old  = u(i)

            tmp = ( (u(i-1) + u(i+1))*invdx2 + &
               (u(i-n) + u(i+n))*invdy2 ) / denom

            u(i) = u_old*(1-factor) + factor*(tmp)

            sum_diff = sum_diff + dabs(u(i) - u_old)
            i = i + 1
         end do
         i = i + 2
      end do

      if (sum_diff < eps) then
         iter_error = sum_diff
         return
      end if
   end do

   iter_error = sum_diff
end subroutine solve_sor_subdomain


subroutine update_interfaces(macrogrid, macrogrid_size_x, macrogrid_size_y, subgrid_size, out_norm)
   implicit none
   integer, intent(in) :: macrogrid_size_x, macrogrid_size_y, subgrid_size
   real*8,  intent(inout) :: macrogrid(macrogrid_size_x,macrogrid_size_y,subgrid_size,subgrid_size)
   real*8,  intent(out)   :: out_norm

   integer :: iX, iY, k
   integer :: global_size_x, global_size_y
   real*8  :: dx, dy
   real*8  :: invdx2, invdy2, denom
   real*8  :: old_value, new_val

   global_size_x = macrogrid_size_x*subgrid_size - (macrogrid_size_x - 1)
   global_size_y = macrogrid_size_y*subgrid_size - (macrogrid_size_y - 1)
   dx = 1.0d0 / dble(global_size_x - 1)
   dy = 1.0d0 / dble(global_size_y - 1)

   out_norm = 0.0d0

   do iX = 1, macrogrid_size_x
      do iY = 1, macrogrid_size_y - 1
         do k = 2, subgrid_size - 1

            old_value = macrogrid(iX, iY, k, subgrid_size)

            new_val = (-macrogrid(iX, iY,   k, subgrid_size-2 ) - macrogrid(iX, iY+1,   k, 3 )  + &
               4.0d0*macrogrid(iX, iY,   k,   subgrid_size-1 ) +4.0d0*macrogrid(iX, iY+1, k, 2 ))/6.0d0


            macrogrid(iX, iY,   k, subgrid_size) = new_val
            macrogrid(iX, iY+1, k, 1)            = new_val

            out_norm = out_norm + dabs(new_val - old_value)
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

            out_norm = out_norm + dabs(new_val - old_value)
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

         out_norm = out_norm + dabs(new_val - old_value)
      end do
   end do
end subroutine update_interfaces


subroutine solve_macrogrid(accumulate_subgrid_error,                       &
   macrogrid, macrogrid_size_x, macrogrid_size_y, &
   subgrid_size, omega,                            &
   max_iter_sor, max_iter_interface,               &
   eps_sor, eps_interface,                         &
   out_iterations)
   implicit none
   logical, intent(in) :: accumulate_subgrid_error
   integer, intent(in) :: macrogrid_size_x, macrogrid_size_y
   integer, intent(in) :: subgrid_size, max_iter_sor, max_iter_interface
   real*8,  intent(in) :: omega, eps_sor, eps_interface
   real*8,  intent(inout) :: macrogrid(macrogrid_size_x,macrogrid_size_y,subgrid_size,subgrid_size)
   integer, intent(out) :: out_iterations

   integer :: iter, iX, iY
   real*8  :: norm_total, norm_sub_sor, norm_bound
   real*8  :: dx, dy

   call compute_dxdy(macrogrid_size_x, macrogrid_size_y, subgrid_size, dx, dy)

   out_iterations = 0

   do iter = 1, max_iter_interface
      norm_total = 0.0d0

      ! (1) Решение SOR в каждой подобласти
      do iX = 1, macrogrid_size_x
         do iY = 1, macrogrid_size_y
            call solve_sor_subdomain(macrogrid(iX,iY,:,:), subgrid_size, &
               omega, eps_sor, max_iter_sor, norm_sub_sor, dx, dy)
            if (accumulate_subgrid_error) then
               norm_total = norm_total + norm_sub_sor
            end if
         end do
      end do

      ! (2) Обновляем интерфейсы (учитывая dx, dy)
      call update_interfaces(macrogrid, macrogrid_size_x, macrogrid_size_y, subgrid_size, norm_bound)

      norm_total = norm_total + norm_bound

      ! Проверяем сходимость
      if (norm_total < eps_interface) then
         out_iterations = iter
         return
      end if
   end do

   out_iterations = max_iter_interface
end subroutine solve_macrogrid


 !-----------------------------------------------------------------------
 ! ПРОСТАЯ РАСЧЁТ ВРЕМЕНИ (EXECUTION TIME) И ВЫЗОВ solve_macrogrid
 !-----------------------------------------------------------------------
subroutine measure_execution_time(accumulate_subgrid_error,                     &
   macrogrid, macrogrid_size_x, macrogrid_size_y,&
   subgrid_size, omega,                          &
   max_iter_sor, max_iter_interface,             &
   eps_sor, eps_interface)
   implicit none
   logical, intent(in) :: accumulate_subgrid_error
   integer, intent(in) :: macrogrid_size_x, macrogrid_size_y
   integer, intent(in) :: subgrid_size, max_iter_sor, max_iter_interface
   real*8,  intent(in) :: omega, eps_sor, eps_interface
   real*8,  intent(inout) :: macrogrid(macrogrid_size_x,macrogrid_size_y,subgrid_size,subgrid_size)

   real*8 :: start_time, end_time, elapsed_time
   integer :: iterations

   call cpu_time(start_time)

   call solve_macrogrid(accumulate_subgrid_error,                &
      macrogrid, macrogrid_size_x, macrogrid_size_y, &
      subgrid_size, omega,                        &
      max_iter_sor, max_iter_interface,           &
      eps_sor, eps_interface,                     &
      iterations)

   call cpu_time(end_time)

   elapsed_time = end_time - start_time
   print *, '  Time to solve macrogrid = ', elapsed_time, ' seconds'
   print *, '  Number of interface iterations = ', iterations
end subroutine measure_execution_time


 !-----------------------------------------------------------------------
 ! ВЫЧИСЛЕНИЕ ОШИБКИ ДЛЯ "КОНСТАНТНОГО" ГРАНИЧНОГО УСЛОВИЯ
 !-----------------------------------------------------------------------
subroutine compute_constant_boundary_error(macrogrid, macrogrid_size_x, macrogrid_size_y, subgrid_size)
   implicit none
   integer, intent(in) :: macrogrid_size_x, macrogrid_size_y, subgrid_size
   real*8,  intent(in) :: macrogrid(macrogrid_size_x,macrogrid_size_y,subgrid_size,subgrid_size)

   integer :: iX, iY, i1, j1
   real*8  :: max_error, diff

   max_error = 0.0d0

   do iX = 1, macrogrid_size_x
      do iY = 1, macrogrid_size_y
         do i1 = 1, subgrid_size
            do j1 = 1, subgrid_size
               diff = abs(macrogrid(iX,iY,i1,j1) - 1.0d0)
               if (diff > max_error) max_error = diff
            end do
         end do
      end do
   end do

   print *, '  Max error (constant boundary) = ', max_error
end subroutine compute_constant_boundary_error


 !-----------------------------------------------------------------------
 ! ВЫЧИСЛЕНИЕ ОШИБКИ ДЛЯ "ЛОГАРИФМИЧЕСКОГО" ГРАНИЧНОГО УСЛОВИЯ
 !-----------------------------------------------------------------------
subroutine compute_logarithmic_boundary_error(macrogrid, macrogrid_size_x, macrogrid_size_y, subgrid_size)
   implicit none
   real*8, parameter :: R1 = 0.1d0, R2 = 1.0d0
   real*8, parameter :: x_min = 0.3d0, y_min = 0.0d0
   real*8, parameter :: len = 0.4d0

   integer, intent(in) :: macrogrid_size_x, macrogrid_size_y, subgrid_size
   real*8,  intent(in) :: macrogrid(macrogrid_size_x,macrogrid_size_y,subgrid_size,subgrid_size)

   integer :: global_size_x, global_size_y
   integer :: iX, iY, i1, j1, lX, lY
   real*8  :: max_error, val_exact, diff

   global_size_x = macrogrid_size_x * subgrid_size - (macrogrid_size_x - 1)
   global_size_y = macrogrid_size_y * subgrid_size - (macrogrid_size_y - 1)

   max_error = 0.0d0

   do iX = 1, macrogrid_size_x
      do iY = 1, macrogrid_size_y
         do i1 = 1, subgrid_size
            do j1 = 1, subgrid_size
               lX = i1 + (iX-1)*subgrid_size - (iX-1)
               lY = j1 + (iY-1)*subgrid_size - (iY-1)

               val_exact = log( sqrt( (x_min + len*(lX-1)/(global_size_x-1))**2 + &
                  (y_min + len*(lY-1)/(global_size_y-1))**2 )*R2/(R1*R1) ) / &
                  log(R2/R1)

               diff = abs(macrogrid(iX,iY,i1,j1) - val_exact)
               if (diff > max_error) max_error = diff
            end do
         end do
      end do
   end do

   print *, '  Max error (logarithmic boundary)  = ', max_error
end subroutine compute_logarithmic_boundary_error
