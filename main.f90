program main
   use macrogrid_boundary_initializer
   use subgrid_boundary_initializer
   use macrogrid_solver
   use subgrid_solver

   implicit none
   integer, parameter :: io = 10

   integer, dimension(8), parameter :: subgrid_sizes = &
      [10, 18, 34, 66, 130, 258, 514, 1026]
   real*8, dimension(8), parameter :: subgrid_omegas = &
      [1.52129d0, 1.69524d0, 1.82919d0, 1.90899d0, 1.95307d0, 1.97618d0, 1.98793d0, 1.99403d0]
   real*8, dimension(8), parameter :: subgrid_const_omegas = &
      [1.52729d0, 1.69504d0, 1.82964d0, 1.90932d0, 1.95305d0, 1.97616d0, 1.98798d0, 1.99394d0]
   integer, dimension(8), parameter :: subgrid_repeats_count = &
      [10000, 10000, 1000, 1000, 100, 100, 10, 10]

   integer, dimension(8), parameter :: macrogrid_sizes = &
      [1, 2, 3, 4, 5, 6, 7, 8]
   real*8, dimension(8), parameter :: macrogrid_omegas_1 = &
      [1.9d0, 1.99d0, 1.999d0, 1.999d0, 1.999d0, 1.999d0, 1.999d0, 1.999d0]
   real*8, dimension(8), parameter :: macrogrid_omegas_2 = &
      [1.97137756696708d0, 1.98302605323770d0, 1.99130449164847d0, &
      1.99618167072694d0, 1.99801933675286d0, 1.99950121735996d0, &
      1.99926204012712d0, 1.999d0]
   real*8, dimension(8), parameter :: macrogrid_omegas_4 = &
      [1.99283687348090d0, 1.99684870120073d0, 1.99836547729607d0, &
      1.99866111661454d0, 1.99955171858458d0, 1.999d0, &
      1.999d0, 1.999d0]
   real*8, dimension(8,8) :: macrogrid_omegas

   macrogrid_omegas(1,:) = macrogrid_omegas_1(:)
   macrogrid_omegas(2,:) = macrogrid_omegas_2(:)
   macrogrid_omegas(3,:) = macrogrid_omegas_1(:)
   macrogrid_omegas(4,:) = macrogrid_omegas_4(:)
   macrogrid_omegas(5,:) = macrogrid_omegas_1(:)
   macrogrid_omegas(6,:) = macrogrid_omegas_1(:)
   macrogrid_omegas(7,:) = macrogrid_omegas_1(:)
   macrogrid_omegas(8,:) = macrogrid_omegas_1(:)

   call macrogrid_test()

contains

   subroutine subgrid_test()
      implicit none

      integer, dimension(6), parameter :: subgrid_configs = &
         [3, 4, 5, 6, 7, 8]

      integer :: sub_cfg

      write(*, *) "Running subgrid tests"

      open(io, file='results.txt', status='unknown', action='write', recl=10000)

      write(io, *) "original_sor"
      do sub_cfg = 1, size(subgrid_configs)

         call run_subgrid_test(subgrid_solver_method = original_sor, &
            eps_subgrid = 1.0d-5, &
            max_iter_subgrid = 100000, &
            subgrid_size_idx = subgrid_configs(sub_cfg), &
            tile_size = 0, subtile_level = 0)

      end do
      write(io, *) " "

      write(io, *) "tiling_sor_1"
      do sub_cfg = 1, size(subgrid_configs)

         call run_subgrid_test(subgrid_solver_method = tiling_sor, &
            eps_subgrid = 1.0d-5, &
            max_iter_subgrid = 100000, &
            subgrid_size_idx = subgrid_configs(sub_cfg), &
            tile_size = 1, subtile_level = 0)

      end do
      write(io, *) " "

      write(io, *) "tiling_sor_2"
      do sub_cfg = 1, size(subgrid_configs)

         call run_subgrid_test(subgrid_solver_method = tiling_sor, &
            eps_subgrid = 1.0d-5, &
            max_iter_subgrid = 100000, &
            subgrid_size_idx = subgrid_configs(sub_cfg), &
            tile_size = 2, subtile_level = 0)

      end do
      write(io, *) " "

      write(io, *) "tiling_sor_4"
      do sub_cfg = 1, size(subgrid_configs)

         call run_subgrid_test(subgrid_solver_method = tiling_sor, &
            eps_subgrid = 1.0d-5, &
            max_iter_subgrid = 100000, &
            subgrid_size_idx = subgrid_configs(sub_cfg), &
            tile_size = 4, subtile_level = 0)

      end do
      write(io, *) " "

      write(io, *) "tiling_sor_8"
      do sub_cfg = 1, size(subgrid_configs)

         call run_subgrid_test(subgrid_solver_method = tiling_sor, &
            eps_subgrid = 1.0d-5, &
            max_iter_subgrid = 100000, &
            subgrid_size_idx = subgrid_configs(sub_cfg), &
            tile_size = 8, subtile_level = 0)

      end do
      write(io, *) " "

      write(io, *) "tiling_sor_16"
      do sub_cfg = 1, size(subgrid_configs)

         call run_subgrid_test(subgrid_solver_method = tiling_sor, &
            eps_subgrid = 1.0d-5, &
            max_iter_subgrid = 100000, &
            subgrid_size_idx = subgrid_configs(sub_cfg), &
            tile_size = 16, subtile_level = 0)

      end do
      write(io, *) " "

      write(io, *) "subtiling_sor_1"
      do sub_cfg = 1, size(subgrid_configs)

         call run_subgrid_test(subgrid_solver_method = subtiling_sor_1, &
            eps_subgrid = 1.0d-5, &
            max_iter_subgrid = 100000, &
            subgrid_size_idx = subgrid_configs(sub_cfg), &
            tile_size = 0, subtile_level = 0)

      end do
      write(io, *) " "

      write(io, *) "subtiling_sor_2"
      do sub_cfg = 1, size(subgrid_configs)

         call run_subgrid_test(subgrid_solver_method = subtiling_sor_2, &
            eps_subgrid = 1.0d-5, &
            max_iter_subgrid = 100000, &
            subgrid_size_idx = subgrid_configs(sub_cfg), &
            tile_size = 0, subtile_level = 0)

      end do
      write(io, *) " "

      write(io, *) "subtiling_sor_4"
      do sub_cfg = 1, size(subgrid_configs)

         call run_subgrid_test(subgrid_solver_method = subtiling_sor_4, &
            eps_subgrid = 1.0d-5, &
            max_iter_subgrid = 100000, &
            subgrid_size_idx = subgrid_configs(sub_cfg), &
            tile_size = 0, subtile_level = 0)

      end do
      write(io, *) " "

      write(io, *) "subtiling_sor_8"
      do sub_cfg = 1, size(subgrid_configs)

         call run_subgrid_test(subgrid_solver_method = subtiling_sor_8, &
            eps_subgrid = 1.0d-5, &
            max_iter_subgrid = 100000, &
            subgrid_size_idx = subgrid_configs(sub_cfg), &
            tile_size = 0, subtile_level = 0)

      end do
      write(io, *) " "

      write(io, *) "subtiling_sor_16"
      do sub_cfg = 1, size(subgrid_configs)

         call run_subgrid_test(subgrid_solver_method = subtiling_sor_16, &
            eps_subgrid = 1.0d-5, &
            max_iter_subgrid = 100000, &
            subgrid_size_idx = subgrid_configs(sub_cfg), &
            tile_size = 0, subtile_level = 0)

      end do
      write(io, *) " "

      close(io)

      write(*, *) "End of subgrid tests"
   end subroutine subgrid_test

   subroutine macrogrid_test()
      implicit none

      integer, dimension(7), parameter :: subgrid_configs = &
         [1, 2, 3, 4, 5, 6, 7]

      integer, dimension(2), parameter :: macrogrid_configs = &
         [2, 4]

      integer :: mac_cfg, sub_cfg

      write(*, *) "Running macrogrid tests"

      open(io, file='results.txt', status='unknown', action='write', recl=10000)

      write(io, *) "many_iter sor_fixed_omega with original_sor"
      do mac_cfg = 1, size(macrogrid_configs)
         do sub_cfg = 1, size(subgrid_configs)

            call run_macrogrid_test(macrogrid_solver_method = sor_fixed_omega, subgrid_solver_method = original_sor, &
               eps_subgrid = 1.0d-8, eps_interface = 1.0d-8, &
               max_iter_interface = 100000, max_iter_subgrid = 100000, &
               macrogrid_size_idx = macrogrid_configs(mac_cfg), subgrid_size_idx = subgrid_configs(sub_cfg), &
               tile_size = 0, subtile_level = 0)

         end do
      end do
      write(io, *) " "

      write(io, *) "sor_fixed_omega with original_sor"
      do mac_cfg = 1, size(macrogrid_configs)
         do sub_cfg = 1, size(subgrid_configs)

            call run_macrogrid_test(macrogrid_solver_method = sor_fixed_omega, subgrid_solver_method = original_sor, &
               eps_subgrid = 1.0d-8, eps_interface = 1.0d-8, &
               max_iter_interface = 100000, max_iter_subgrid = 1, &
               macrogrid_size_idx = macrogrid_configs(mac_cfg), subgrid_size_idx = subgrid_configs(sub_cfg), &
               tile_size = 0, subtile_level = 0)

         end do
      end do
      write(io, *) " "

      write(io, *) "sor_fixed_omega with tiling_sor"
      do mac_cfg = 1, size(macrogrid_configs)
         do sub_cfg = 1, size(subgrid_configs)

            call run_macrogrid_test(macrogrid_solver_method = sor_fixed_omega, subgrid_solver_method = tiling_sor, &
               eps_subgrid = 1.0d-8, eps_interface = 1.0d-8, &
               max_iter_interface = 100000, max_iter_subgrid = 1, &
               macrogrid_size_idx = macrogrid_configs(mac_cfg), subgrid_size_idx = subgrid_configs(sub_cfg), &
               tile_size = 8, subtile_level = 0)

         end do
      end do
      write(io, *) " "

      write(io, *) "conjugate_residuals with original_sor"
      do mac_cfg = 1, size(macrogrid_configs)
         do sub_cfg = 1, size(subgrid_configs)

            call run_macrogrid_test(macrogrid_solver_method = conjugate_residuals, subgrid_solver_method = original_sor, &
               eps_subgrid = 1.0d-8, eps_interface = 1.0d-8, &
               max_iter_interface = 100000, max_iter_subgrid = 100000, &
               macrogrid_size_idx = macrogrid_configs(mac_cfg), subgrid_size_idx = subgrid_configs(sub_cfg), &
               tile_size = 0, subtile_level = 0)

         end do
      end do
      write(io, *) " "

      write(io, *) "conjugate_residuals with tiling_sor"
      do mac_cfg = 1, size(macrogrid_configs)
         do sub_cfg = 1, size(subgrid_configs)

            call run_macrogrid_test(macrogrid_solver_method = conjugate_residuals, subgrid_solver_method = tiling_sor, &
               eps_subgrid = 1.0d-8, eps_interface = 1.0d-8, &
               max_iter_interface = 100000, max_iter_subgrid = 100000, &
               macrogrid_size_idx = macrogrid_configs(mac_cfg), subgrid_size_idx = subgrid_configs(sub_cfg), &
               tile_size = 8, subtile_level = 0)

         end do
      end do
      write(io, *) " "

      write(io, *) "conjugate_residuals with subtiling_sor"
      do mac_cfg = 1, size(macrogrid_configs)
         do sub_cfg = 1, size(subgrid_configs)

            call run_macrogrid_test(macrogrid_solver_method = conjugate_residuals, subgrid_solver_method = subtiling_sor_4, &
               eps_subgrid = 1.0d-8, eps_interface = 1.0d-8, &
               max_iter_interface = 100000, max_iter_subgrid = 100000, &
               macrogrid_size_idx = macrogrid_configs(mac_cfg), subgrid_size_idx = subgrid_configs(sub_cfg), &
               tile_size = 0, subtile_level = 0)

         end do
      end do
      write(io, *) " "

      close(io)

      write(*, *) "End of macrogrid tests"
   end subroutine macrogrid_test

   subroutine run_subgrid_test(subgrid_solver_method, &
      eps_subgrid, max_iter_subgrid, &
      subgrid_size_idx, tile_size, subtile_level)
      implicit none
      procedure(i_subgrid_solver_method) :: subgrid_solver_method
      integer, intent(in) :: subgrid_size_idx
      integer, intent(in) :: tile_size, subtile_level
      real*8, intent(in) :: eps_subgrid
      integer, intent(in) :: max_iter_subgrid

      real*8, allocatable :: subgrid(:)

      real*8 :: error, time, sum_time
      integer :: iter, i

      allocate(subgrid(subgrid_sizes(subgrid_size_idx)*subgrid_sizes(subgrid_size_idx)))

      sum_time = 0.0d0

      do i = 1, subgrid_repeats_count(subgrid_size_idx)

         call initialize_subgrid_boundary(subgrid, subgrid_sizes(subgrid_size_idx))

         call configure_subgrid_solver(eps_subgrid, max_iter_subgrid, &
            subgrid_omegas(subgrid_size_idx), tile_size, subtile_level)

         call run_subgrid_solver(subgrid, subgrid_sizes(subgrid_size_idx), subgrid_solver_method)

         call get_subgrid_solver_results(time, iter)

         call compute_subgrid_boundary_error(subgrid, subgrid_sizes(subgrid_size_idx), error)

         sum_time = sum_time + time

      end do

      time = sum_time / real(subgrid_repeats_count(subgrid_size_idx))

      write(io,*) "Subgrid_size&Error&Run_time&Iterations#", &
         subgrid_sizes(subgrid_size_idx), "&", error, "&", time, "&", iter

      deallocate(subgrid)

   end subroutine run_subgrid_test

   subroutine run_macrogrid_test(macrogrid_solver_method, subgrid_solver_method, &
      eps_subgrid, eps_interface, max_iter_interface, max_iter_subgrid, &
      macrogrid_size_idx, subgrid_size_idx, tile_size, subtile_level)
      implicit none
      procedure(interface_macrogrid_solver_method) :: macrogrid_solver_method
      procedure(interface_subgrid_solver_method) :: subgrid_solver_method
      integer, intent(in) :: macrogrid_size_idx, subgrid_size_idx
      integer, intent(in) :: tile_size, subtile_level
      real*8, intent(in) :: eps_subgrid, eps_interface
      integer, intent(in) :: max_iter_interface, max_iter_subgrid

      real*8, allocatable :: macrogrid(:, :, :, :)

      real*8 :: error, time
      integer :: iter

      allocate(macrogrid(macrogrid_sizes(macrogrid_size_idx), macrogrid_sizes(macrogrid_size_idx), &
         subgrid_sizes(subgrid_size_idx), subgrid_sizes(subgrid_size_idx)))

      call initialize_macrogrid_boundary(macrogrid, &
         macrogrid_sizes(macrogrid_size_idx), macrogrid_sizes(macrogrid_size_idx), subgrid_sizes(subgrid_size_idx))

      call configure_subgrid_solver(eps_subgrid, max_iter_subgrid, &
         subgrid_omegas(subgrid_size_idx), tile_size, subtile_level)

      call configure_macrogrid_solver(eps_interface, max_iter_interface, &
         macrogrid_omegas(macrogrid_size_idx, subgrid_size_idx))

      call run_macrogrid_solver(macrogrid, macrogrid_sizes(macrogrid_size_idx), macrogrid_sizes(macrogrid_size_idx), &
         subgrid_sizes(subgrid_size_idx), macrogrid_solver_method, subgrid_solver_method)

      call get_macrogrid_solver_results(time, iter)

      call compute_macrogrid_boundary_error(macrogrid, &
         macrogrid_sizes(macrogrid_size_idx), macrogrid_sizes(macrogrid_size_idx), &
         subgrid_sizes(subgrid_size_idx), error)

      write(io,*) "Macrogrid_size&Subgrid_size&Error&Run_time&Iterations#", &
         macrogrid_sizes(macrogrid_size_idx), "&", subgrid_sizes(subgrid_size_idx), &
         "&", error, "&", time, "&", iter

      deallocate(macrogrid)

   end subroutine run_macrogrid_test

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

         call initialize_macrogrid_boundary(macrogrid, &
            macrogrid_sizes(macrogrid_size_idx), macrogrid_sizes(macrogrid_size_idx), subgrid_sizes(subgrid_size_idx))
         call configure_subgrid_solver(eps_subgrid, max_iter_subgrid, &
            subgrid_omegas(subgrid_size_idx), tile_size, subtile_level)
         call configure_macrogrid_solver(eps_interface, max_iter_interface, &
            wa)
         call run_macrogrid_solver(macrogrid, macrogrid_sizes(macrogrid_size_idx), macrogrid_sizes(macrogrid_size_idx), &
            subgrid_sizes(subgrid_size_idx), sor_fixed_omega, tiling_sor)
         call get_macrogrid_solver_results(time, iter)
         iter1 = iter
         time1 = time
         time = 0.0d0
         iter = 0

         call initialize_macrogrid_boundary(macrogrid, &
            macrogrid_sizes(macrogrid_size_idx), macrogrid_sizes(macrogrid_size_idx), subgrid_sizes(subgrid_size_idx))
         call configure_subgrid_solver(eps_subgrid, max_iter_subgrid, &
            subgrid_omegas(subgrid_size_idx), tile_size, subtile_level)
         call configure_macrogrid_solver(eps_interface, max_iter_interface, &
            wb)
         call run_macrogrid_solver(macrogrid, macrogrid_sizes(macrogrid_size_idx), macrogrid_sizes(macrogrid_size_idx), &
            subgrid_sizes(subgrid_size_idx), sor_fixed_omega, tiling_sor)
         call get_macrogrid_solver_results(time, iter)
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
