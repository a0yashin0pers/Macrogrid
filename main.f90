program main
   use macrogrid_initializer
   use macrogrid_solver
   use subgrid_solver

   implicit none
   integer, parameter :: io = 10

   integer, dimension(8), parameter :: subgrid_sizes = &
      [10, 18, 34, 66, 130, 258, 514, 1026]
   real*8, dimension(8), parameter :: subgrid_omegas = &
      [1.52129d0, 1.69524d0, 1.82919d0, 1.90899d0, 1.95307d0, 1.97618d0, 1.98793d0, 1.99403d0]
   real*8, dimension(8), parameter :: subgrid_omegas_const = &
      [1.52729d0, 1.69504d0, 1.82964d0, 1.90932d0, 1.95305d0, 1.97616d0, 1.98798d0, 1.99394d0]

   integer, dimension(3), parameter :: subgrid_configs = &
      [1, 2, 3]

   integer, dimension(8), parameter :: macrogrid_sizes = &
      [1, 2, 3, 4, 5, 6, 7, 8]

   integer, dimension(2), parameter :: macrogrid_configs = &
      [2, 4]

   integer :: mac_cfg, sub_cfg
   real*8 :: omega

   omega = 1.0d0

   write(*, *) "Running tests"

   open(io, file='results.txt', status='unknown', action='write')

   write(io, *) "simple_iteration with tiling_sor"
   do mac_cfg = 1, size(macrogrid_configs)
      do sub_cfg = 1, size(subgrid_configs)

         call run_test(macrogrid_solver_method = simple_iteration, subgrid_solver_method = tiling_sor, &
            eps_subgrid = 1.0d-8, eps_interface = 1.0d-8, &
            max_iter_interface = 100000, max_iter_subgrid = 1, &
            macrogrid_size_idx = macrogrid_configs(mac_cfg), subgrid_size_idx = subgrid_configs(sub_cfg), &
            macrogrid_omega = omega, tile_size = 8, subtile_level = 8)

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
            macrogrid_omega = omega, tile_size = 8, subtile_level = 8)

      end do
   end do
   write(io, *) " "

   write(io, *) "sor_fixed_omega with tiling_sor"
   do mac_cfg = 1, size(macrogrid_configs)
      do sub_cfg = 1, size(subgrid_configs)

         call golden_section(macrogrid_configs(mac_cfg), subgrid_configs(sub_cfg), omega)

         call run_test(macrogrid_solver_method = sor_fixed_omega, subgrid_solver_method = tiling_sor, &
            eps_subgrid = 1.0d-8, eps_interface = 1.0d-8, &
            max_iter_interface = 100000, max_iter_subgrid = 1, &
            macrogrid_size_idx = macrogrid_configs(mac_cfg), subgrid_size_idx = subgrid_configs(sub_cfg), &
            macrogrid_omega = omega, tile_size = 8, subtile_level = 8)

      end do
   end do
   write(io, *) " "

   close(io)

   write(*, *) "End of tests"

contains

   subroutine run_test(macrogrid_solver_method, subgrid_solver_method, &
      eps_subgrid, eps_interface, max_iter_interface, max_iter_subgrid, &
      macrogrid_size_idx, subgrid_size_idx, macrogrid_omega, tile_size, subtile_level)
      implicit none
      procedure(interface_macrogrid_solver_method) :: macrogrid_solver_method
      procedure(interface_subgrid_solver_method) :: subgrid_solver_method
      integer, intent(in) :: macrogrid_size_idx, subgrid_size_idx
      integer, intent(in) :: tile_size, subtile_level
      real*8, intent(in) :: eps_subgrid, eps_interface
      real*8, intent(inout) :: macrogrid_omega
      integer, intent(in) :: max_iter_interface, max_iter_subgrid

      real*8, allocatable :: macrogrid(:, :, :, :)

      real*8 :: error, time
      integer :: iter
      
      allocate(macrogrid(macrogrid_sizes(macrogrid_size_idx), macrogrid_sizes(macrogrid_size_idx), &
         subgrid_sizes(subgrid_size_idx), subgrid_sizes(subgrid_size_idx)))

      call initialize_boundary(macrogrid, &
         macrogrid_sizes(macrogrid_size_idx), macrogrid_sizes(macrogrid_size_idx), subgrid_sizes(subgrid_size_idx))

      call configure_subgrid_solver(eps_subgrid, max_iter_subgrid, &
         subgrid_omegas(subgrid_size_idx), tile_size, subtile_level)

      call configure_macrogrid_solver(eps_interface, max_iter_interface, &
         macrogrid_omega)

      call run(macrogrid, macrogrid_sizes(macrogrid_size_idx), macrogrid_sizes(macrogrid_size_idx), &
         subgrid_sizes(subgrid_size_idx), macrogrid_solver_method, subgrid_solver_method)

      call get_results(time, iter)

      call compute_boundary_error(macrogrid, &
         macrogrid_sizes(macrogrid_size_idx), macrogrid_sizes(macrogrid_size_idx), &
         subgrid_sizes(subgrid_size_idx), error)

      call get_omega(macrogrid_omega)

      write(io,*) "Macrogrid_size&Subgrid_size&Error&Run_time&Iterations&Omega#", &
         macrogrid_sizes(macrogrid_size_idx), "&", subgrid_sizes(subgrid_size_idx), &
         "&", error, "&", time, "&", iter, "&", macrogrid_omega

      deallocate(macrogrid)

   end subroutine run_test

   subroutine golden_section(macrogrid_size_idx, subgrid_size_idx, macrogrid_omega)
      implicit none
      integer, intent(in) :: macrogrid_size_idx, subgrid_size_idx
      real*8, intent(out) :: macrogrid_omega
      
      integer :: tile_size = 8, subtile_level = 8
      real*8 :: eps_subgrid = 1.0d-8, eps_interface = 1.0d-8, eps_golden = 1.0d-5
      integer :: max_iter_interface = 100000, max_iter_subgrid = 1

      real*8, allocatable :: macrogrid(:, :, :, :)

      real*8 :: error, time, time1, time2
      integer :: iter,iter1,iter2

      real*8 :: k, a, b, wa, wb

      k = 0.6180339887d0
      a = 1.0d0
      b = 2.0d0
      wa = 1.0d0
      wb = 2.0d0

      allocate(macrogrid(macrogrid_sizes(macrogrid_size_idx), macrogrid_sizes(macrogrid_size_idx), &
         subgrid_sizes(subgrid_size_idx), subgrid_sizes(subgrid_size_idx)))

      macrogrid_omega = 1.5d0

      do 
         if (abs(wa-wb) < eps_golden) then 
         macrogrid_omega = (wb+wa)/2
         exit
         end if

         wa = a + (1-k)*(b-a)
         wb = a + k*(b-a)

         call initialize_boundary(macrogrid, &
            macrogrid_sizes(macrogrid_size_idx), macrogrid_sizes(macrogrid_size_idx), subgrid_sizes(subgrid_size_idx))
         call configure_subgrid_solver(eps_subgrid, max_iter_subgrid, &
            subgrid_omegas(subgrid_size_idx), tile_size, subtile_level)
         call configure_macrogrid_solver(eps_interface, max_iter_interface, &
            wa)
         call run(macrogrid, macrogrid_sizes(macrogrid_size_idx), macrogrid_sizes(macrogrid_size_idx), &
            subgrid_sizes(subgrid_size_idx), sor_fixed_omega, tiling_sor)
         call get_results(time, iter)
         iter1 = iter
         time1 = time
         time = 0.0d0
         iter = 0

         call initialize_boundary(macrogrid, &
            macrogrid_sizes(macrogrid_size_idx), macrogrid_sizes(macrogrid_size_idx), subgrid_sizes(subgrid_size_idx))
         call configure_subgrid_solver(eps_subgrid, max_iter_subgrid, &
            subgrid_omegas(subgrid_size_idx), tile_size, subtile_level)
         call configure_macrogrid_solver(eps_interface, max_iter_interface, &
            wb)
         call run(macrogrid, macrogrid_sizes(macrogrid_size_idx), macrogrid_sizes(macrogrid_size_idx), &
            subgrid_sizes(subgrid_size_idx), sor_fixed_omega, tiling_sor)
         call get_results(time, iter)
         iter2 = iter
         time2 = time
         time = 0.0d0
         iter = 0

         if (iter1 > iter2) then 
         a = wa
         else 
         b = wb
         end if

      end do

      call initialize_boundary(macrogrid, &
            macrogrid_sizes(macrogrid_size_idx), macrogrid_sizes(macrogrid_size_idx), subgrid_sizes(subgrid_size_idx))
      call configure_subgrid_solver(eps_subgrid, max_iter_subgrid, &
            subgrid_omegas(subgrid_size_idx), tile_size, subtile_level)
      call configure_macrogrid_solver(eps_interface, max_iter_interface, &
            macrogrid_omega)
      call run(macrogrid, macrogrid_sizes(macrogrid_size_idx), macrogrid_sizes(macrogrid_size_idx), &
            subgrid_sizes(subgrid_size_idx), sor_fixed_omega, tiling_sor)
      call get_results(time, iter)

      call compute_boundary_error(macrogrid, &
         macrogrid_sizes(macrogrid_size_idx), macrogrid_sizes(macrogrid_size_idx), &
         subgrid_sizes(subgrid_size_idx), error)

      deallocate(macrogrid)

   end subroutine golden_section

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


end program main