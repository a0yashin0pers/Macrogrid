program main
   use macrogrid_initializer
   use macrogrid_solvers
   use subgrid_solvers
   
   implicit none
   integer, parameter :: io = 10

   integer, parameter :: SUBGRID_SIZES_COUNT = 8
   integer, dimension(SUBGRID_SIZES_COUNT), parameter :: subgrid_sizes = &
      [10, 18, 34, 66, 130, 258, 514, 1026]

   real*8, dimension(SUBGRID_SIZES_COUNT), parameter :: factors_const = &
      [1.52729d0, 1.69504d0, 1.82964d0, 1.90932d0, 1.95305d0, 1.97616d0, 1.98798d0, 1.99394d0]

   real*8, dimension(SUBGRID_SIZES_COUNT), parameter :: factors_log = &
      [1.52129d0, 1.69524d0, 1.82919d0, 1.90899d0, 1.95307d0, 1.97618d0, 1.98793d0, 1.99403d0]

   integer, parameter :: SUBGRID_CFGS_COUNT = 5
   integer, dimension(SUBGRID_CFGS_COUNT), parameter :: subgrid_configs = &
      [1, 2, 3, 4, 5]

   integer, parameter :: MACRIGRID_CFGS_COUNT = 2
   integer, dimension(2, MACRIGRID_CFGS_COUNT), parameter :: macrogrid_configs = reshape([ &
      2, 2, &
      4, 4 &
      ], shape=[2, MACRIGRID_CFGS_COUNT])

   integer, parameter :: max_iter_subgrid   = 100000
   integer, parameter :: max_iter_interface = 100000
   real*8,  parameter :: eps_subgrid        = 1.0d-5
   real*8,  parameter :: eps_interface      = 1.0d-5

   real*8, allocatable :: macrogrid(:, :, :, :)
   integer :: i_cfg, i_size, cfg_x, cfg_y, sub_sz, iter
   real*8  :: omega, time, error

   allocate(macrogrid(10, 10, 1026, 1026))
   open(io, file='results.txt', status='unknown', action='write')

   write(io, *) "simple_iteration_one_iter with original_sor"
   do i_cfg = 1, MACRIGRID_CFGS_COUNT
      cfg_x = macrogrid_configs(1, i_cfg)
      cfg_y = macrogrid_configs(2, i_cfg)
      do i_size = 1, SUBGRID_CFGS_COUNT
         sub_sz = subgrid_sizes(subgrid_configs(i_size))
         omega = factors_log(subgrid_configs(i_size))
         call run_test(io, simple_iteration, original_sor, initialize_logarithmic_boundary, compute_logarithmic_boundary_error, &
            macrogrid, cfg_x, cfg_y, sub_sz, omega, eps_subgrid, eps_interface, 1, max_iter_interface)
      end do
   end do
   write(io, *) " "

   write(io, *) "simple_iteration_one_iter with tiling_sor"
   do i_cfg = 1, MACRIGRID_CFGS_COUNT
      cfg_x = macrogrid_configs(1, i_cfg)
      cfg_y = macrogrid_configs(2, i_cfg)
      do i_size = 1, SUBGRID_CFGS_COUNT
         sub_sz = subgrid_sizes(subgrid_configs(i_size))
         omega = factors_log(subgrid_configs(i_size))
         call run_test(io, simple_iteration, tiling_sor, initialize_logarithmic_boundary, compute_logarithmic_boundary_error, &
            macrogrid, cfg_x, cfg_y, sub_sz, omega, eps_subgrid, eps_interface, 1, max_iter_interface)
      end do
   end do
   write(io, *) " "

   write(io, *) "conjugate_residuals with original_sor"
   do i_cfg = 1, MACRIGRID_CFGS_COUNT
      cfg_x = macrogrid_configs(1, i_cfg)
      cfg_y = macrogrid_configs(2, i_cfg)
      do i_size = 1, SUBGRID_CFGS_COUNT
         sub_sz = subgrid_sizes(subgrid_configs(i_size))
         omega = factors_log(subgrid_configs(i_size))
         call run_test(io, conjugate_residuals, original_sor, initialize_logarithmic_boundary, compute_logarithmic_boundary_error, &
            macrogrid, cfg_x, cfg_y, sub_sz, omega, eps_subgrid, eps_interface, max_iter_subgrid, max_iter_interface)
      end do
   end do
   write(io, *) " "

   write(io, *) "conjugate_residuals with tiling_sor"
   do i_cfg = 1, MACRIGRID_CFGS_COUNT
      cfg_x = macrogrid_configs(1, i_cfg)
      cfg_y = macrogrid_configs(2, i_cfg)
      do i_size = 1, SUBGRID_CFGS_COUNT
         sub_sz = subgrid_sizes(subgrid_configs(i_size))
         omega = factors_log(subgrid_configs(i_size))
         call run_test(io, conjugate_residuals, tiling_sor, initialize_logarithmic_boundary, compute_logarithmic_boundary_error, &
            macrogrid, cfg_x, cfg_y, sub_sz, omega, eps_subgrid, eps_interface, max_iter_subgrid, max_iter_interface)
      end do
   end do
   write(io, *) " "

   write(io, *) "conjugate_residuals with subtiling_sor"
   do i_cfg = 1, MACRIGRID_CFGS_COUNT
      cfg_x = macrogrid_configs(1, i_cfg)
      cfg_y = macrogrid_configs(2, i_cfg)
      do i_size = 1, SUBGRID_CFGS_COUNT
         sub_sz = subgrid_sizes(subgrid_configs(i_size))
         omega = factors_log(subgrid_configs(i_size))
         call run_test(io, conjugate_residuals, subtiling_sor_4, initialize_logarithmic_boundary, compute_logarithmic_boundary_error, &
            macrogrid, cfg_x, cfg_y, sub_sz, omega, eps_subgrid, eps_interface, max_iter_subgrid, max_iter_interface)
      end do
   end do
   write(io, *) " "

   deallocate(macrogrid)
   close(io)

end program main

subroutine run_test(io, macrigrid_solver, subgrid_solver, &
   boundary_condition, compute_error,                 &
   macrogrid, macrogrid_size_x, macrogrid_size_y,     &
   subgrid_size, omega,                               &
   eps_subgrid, eps_interface,                        &
   max_iter_subgrid, max_iter_interface)
   implicit none
   external :: macrigrid_solver, subgrid_solver, boundary_condition, compute_error
   integer, intent(in) :: io, macrogrid_size_x, macrogrid_size_y, subgrid_size, max_iter_subgrid, max_iter_interface
   real*8,  intent(inout) :: macrogrid(macrogrid_size_x,macrogrid_size_y,subgrid_size,subgrid_size)
   real*8,  intent(in)   :: omega, eps_subgrid, eps_interface
   real*8   :: error, time
   integer  :: iter

   call boundary_condition(macrogrid, macrogrid_size_x, macrogrid_size_y, subgrid_size)
   call macrigrid_solver(subgrid_solver, macrogrid, macrogrid_size_x, macrogrid_size_y, subgrid_size, omega, eps_subgrid, eps_interface, &
      max_iter_subgrid, max_iter_interface, error, time, iter)
   call compute_error(macrogrid, macrogrid_size_x, macrogrid_size_y, subgrid_size, error)

   write(io,'(A,I5,A,I5,A,I5,A,F14.8,A,F14.8,A,I6)') &
      "Macrogrid_size_x&Macrogrid_size_y&Subgrid_size&Error&Run_time&Iterations#", &
      macrogrid_size_x, "&", macrogrid_size_y, "&", subgrid_size, "&", error, "&", time, "&", iter

end subroutine run_test

subroutine print_macrogrid(macrogrid, macrogrid_size, subgrid_size)
   implicit none
   integer, intent(in) :: macrogrid_size, subgrid_size
   real*8, intent(in) :: macrogrid(macrogrid_size, macrogrid_size, subgrid_size, subgrid_size)
   integer :: i, j, x, y

   do i = 1, macrogrid_size
      do y = 1, subgrid_size
         do j = 1, macrogrid_size
            do x = 1, subgrid_size
               write(*, '(10F8.3)', advance='no') macrogrid(j, i, x, y)
            end do
         end do
         print *, ' '
      end do
   end do
   print *, ' '

end subroutine print_macrogrid