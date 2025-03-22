program main
   use macrogrid_initializer
   use macrogrid_solver
   use subgrid_solver

   implicit none
   integer, parameter :: io = 10

   integer, dimension(8), parameter :: subgrid_sizes = &
      [10, 18, 34, 66, 130, 258, 514, 1026]
   real*8, dimension(8), parameter :: subgrid_omegas_const = &
      [1.52729d0, 1.69504d0, 1.82964d0, 1.90932d0, 1.95305d0, 1.97616d0, 1.98798d0, 1.99394d0]
   real*8, dimension(8), parameter :: subgrid_omegas_log = &
      [1.52129d0, 1.69524d0, 1.82919d0, 1.90899d0, 1.95307d0, 1.97618d0, 1.98793d0, 1.99403d0]

   integer, dimension(4), parameter :: subgrid_configs = &
      [1, 2, 3, 4]

   integer, dimension(8), parameter :: macrogrid_sizes = &
      [1, 2, 3, 4, 5, 6, 7, 8]
   real*8, dimension(8), parameter :: macrogrid_omegas_const = &
      [1.52729d0, 1.69504d0, 1.82964d0, 1.90932d0, 1.95305d0, 1.97616d0, 1.98798d0, 1.99394d0]
   real*8, dimension(8), parameter :: macrogrid_omegas_log = &
      [1.52129d0, 1.69524d0, 1.82919d0, 1.90899d0, 1.95307d0, 1.97618d0, 1.98793d0, 1.99403d0]

   integer, dimension(2), parameter :: macrogrid_configs = &
      [2, 4]

   integer :: mac_cfg, sub_cfg

   write(*, *) "Running tests"

   open(io, file='results.txt', status='unknown', action='write')

   write(io, *) "simple_iteration with tiling_sor"
   do mac_cfg = 1, size(macrogrid_configs)
      do sub_cfg = 1, size(subgrid_configs)

         call run_test(macrogrid_solver_method = simple_iteration, subgrid_solver_method = tiling_sor, &
            eps_subgrid = 1.0d-8, eps_interface = 1.0d-8, &
            max_iter_interface = 100000, max_iter_subgrid = 1, &
            macrogrid_size_idx = macrogrid_configs(mac_cfg), subgrid_size_idx = subgrid_configs(sub_cfg), &
            macrogrid_omega_idx = macrogrid_configs(mac_cfg), omega_subgrid_idx = subgrid_configs(sub_cfg), &
            tile_size = 4, subtile_level = 4)

      end do
   end do
   write(io, *) " "

   write(io, *) "sor with tiling_sor"
   do mac_cfg = 1, size(macrogrid_configs)
      do sub_cfg = 1, size(subgrid_configs)

         call run_test(macrogrid_solver_method = sor, subgrid_solver_method = tiling_sor, &
            eps_subgrid = 1.0d-8, eps_interface = 1.0d-8, &
            max_iter_interface = 100000, max_iter_subgrid = 1, &
            macrogrid_size_idx = macrogrid_configs(mac_cfg), subgrid_size_idx = subgrid_configs(sub_cfg), &
            macrogrid_omega_idx = macrogrid_configs(mac_cfg), omega_subgrid_idx = subgrid_configs(sub_cfg), &
            tile_size = 4, subtile_level = 4)

      end do
   end do
   write(io, *) " "

   write(io, *) "conjugate_residuals with tiling_sor"
   do mac_cfg = 1, size(macrogrid_configs)
      do sub_cfg = 1, size(subgrid_configs)

         call run_test(macrogrid_solver_method = conjugate_residuals, subgrid_solver_method = tiling_sor, &
            eps_subgrid = 1.0d-5, eps_interface = 1.0d-5, &
            max_iter_interface = 100000, max_iter_subgrid = 100000, &
            macrogrid_size_idx = macrogrid_configs(mac_cfg), subgrid_size_idx = subgrid_configs(sub_cfg), &
            macrogrid_omega_idx = macrogrid_configs(mac_cfg), omega_subgrid_idx = subgrid_configs(sub_cfg), &
            tile_size = 4, subtile_level = 4)

      end do
   end do
   write(io, *) " "

   write(io, *) "conjugate_residuals with subtiling_sor"
   do mac_cfg = 1, size(macrogrid_configs)
      do sub_cfg = 1, size(subgrid_configs)

         call run_test(macrogrid_solver_method = conjugate_residuals, subgrid_solver_method = subtiling_sor, &
            eps_subgrid = 1.0d-5, eps_interface = 1.0d-5, &
            max_iter_interface = 100000, max_iter_subgrid = 100000, &
            macrogrid_size_idx = macrogrid_configs(mac_cfg), subgrid_size_idx = subgrid_configs(sub_cfg), &
            macrogrid_omega_idx = macrogrid_configs(mac_cfg), omega_subgrid_idx = subgrid_configs(sub_cfg), &
            tile_size = 4, subtile_level = 4)

      end do
   end do
   write(io, *) " "

   close(io)

   write(*, *) "End of tests"

contains

   subroutine run_test (macrogrid_solver_method, subgrid_solver_method, &
      eps_subgrid, eps_interface, max_iter_interface, max_iter_subgrid, &
      macrogrid_size_idx, subgrid_size_idx, macrogrid_omega_idx, omega_subgrid_idx, tile_size, subtile_level)
      implicit none
      procedure(interface_macrogrid_solver_method) :: macrogrid_solver_method
      procedure(interface_subgrid_solver_method) :: subgrid_solver_method
      external :: macrogrid_solver_method, subgrid_solver_method
      integer, intent(in) :: macrogrid_size_idx, subgrid_size_idx, macrogrid_omega_idx, omega_subgrid_idx
      integer, intent(in) :: tile_size, subtile_level
      real*8, intent(in) :: eps_subgrid, eps_interface
      integer, intent(in) :: max_iter_interface, max_iter_subgrid

      real*8, allocatable :: macrogrid(:, :, :, :)

      real*8 :: error, time
      integer :: iter
      
      allocate(macrogrid(macrogrid_sizes(macrogrid_size_idx), macrogrid_sizes(macrogrid_size_idx), &
         subgrid_sizes(subgrid_size_idx), subgrid_sizes(subgrid_size_idx)))

      call initialize_logarithmic_boundary(macrogrid, &
         macrogrid_sizes(macrogrid_size_idx), macrogrid_sizes(macrogrid_size_idx), subgrid_sizes(subgrid_size_idx))

      call configure_subgrid_solver(eps_subgrid, max_iter_subgrid, &
         subgrid_omegas_log(omega_subgrid_idx), tile_size, subtile_level)

      call configure_macrogrid_solver(eps_interface, max_iter_interface, &
         macrogrid_omegas_log(macrogrid_omega_idx))

      call run(macrogrid, macrogrid_sizes(macrogrid_size_idx), macrogrid_sizes(macrogrid_size_idx), &
         subgrid_sizes(subgrid_size_idx), macrogrid_solver_method, subgrid_solver_method)

      call get_results(time, iter)

      call compute_logarithmic_boundary_error(macrogrid, &
         macrogrid_sizes(macrogrid_size_idx), macrogrid_sizes(macrogrid_size_idx), &
         subgrid_sizes(subgrid_size_idx), error)

      write(io,'(A,I5,A,I5,A,F14.8,A,F14.8,A,I6)') &
         "Macrogrid_size_x&Macrogrid_size_y&Subgrid_size&Error&Run_time&Iterations#", &
         macrogrid_sizes(macrogrid_size_idx), "&", subgrid_sizes(subgrid_size_idx), "&", error, "&", time, "&", iter

      deallocate(macrogrid)

   end subroutine run_test
end program main

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