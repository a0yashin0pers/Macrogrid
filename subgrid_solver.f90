module subgrid_solver
   implicit none
   private

   public :: configure_subgrid_solver

   public :: original_sor
   public :: tiling_sor
   public :: subtiling_sor
   public :: subtiling_sor_test_version
   public :: subtiling_sor_4
   public :: subtiling_sor_8
   public :: subtiling_sor_16

   real*8 :: eps, omega
   integer :: max_iter
   integer :: tile_size, subtile_level

contains

   subroutine configure_subgrid_solver(new_eps, new_max_iter, new_factor, new_tile_size, new_subtile_level)
      implicit none
      real*8, intent(in) :: new_eps, new_factor
      integer, intent(in) :: new_max_iter, new_tile_size, new_subtile_level

      eps = new_eps
      max_iter = new_max_iter
      omega = new_factor
      tile_size = new_tile_size
      subtile_level = new_subtile_level

   end subroutine configure_subgrid_solver

   subroutine original_sor(u, u_size, dx, dy)
      implicit none
      integer, intent(in) :: u_size
      real*8, intent(inout) :: u(u_size*u_size)
      real*8, intent(in) :: dx, dy

      integer :: iter, i, l0, l1, l2, l3
      real*8 :: invdx2, invdy2, f0, f1, u_old, error

      invdx2 = 1.0d0/(dx*dx)
      invdy2 = 1.0d0/(dy*dy)

      f0 = 1.0d0 - omega
      f1 = omega / (2.0d0*invdx2 + 2.0d0*invdy2)

      do iter = 1, max_iter

         error = 0.0d0
         i = u_size + 2

         do l1 = 3, u_size
            do l0 = 3, u_size

               u_old = u(i)
               u(i) = f0*u_old + &
                  f1*((u(i-1) + u(i+1))*invdx2 + &
                  (u(i-u_size) + u(i+u_size))*invdy2)

               error = error + abs(u(i) - u_old)
               i = i + 1

            end do
            i = i + 2
         end do

         if (error < eps) then
            return
         end if

      end do
   end subroutine original_sor

   subroutine tiling_sor(u, u_size, dx, dy)
      implicit none
      integer, intent(in) :: u_size
      real*8, intent(inout) :: u(u_size*u_size)
      real*8, intent(in) :: dx, dy

      integer :: tile_count, iter, i, l0, l1, l2, l3
      real*8 :: invdx2, invdy2, f0, f1, u_old, error

      tile_count = (u_size - 2)/tile_size

      invdx2 = 1.0d0/(dx*dx)
      invdy2 = 1.0d0/(dy*dy)

      f0 = 1.0d0 - omega
      f1 = omega / (2.0d0*invdx2 + 2.0d0*invdy2)

      do iter = 1, max_iter

         error = 0.0d0
         i = u_size + 2

         do l3 = 1, tile_count
            do l2 = 1, tile_count
               do l1 = 1, tile_size
                  do l0 = 1, tile_size

                     u_old = u(i)
                     u(i) = f0*u_old + &
                        f1*((u(i-1) + u(i+1))*invdx2 + &
                        (u(i-u_size) + u(i+u_size))*invdy2)

                     error = error + abs(u(i) - u_old)
                     i = i + 1

                  end do
                  i = i + u_size - tile_size
               end do
               i = i - tile_size*(u_size - 1)
            end do
            i = i + u_size*(tile_size - 1) + 2
         end do

         if (error < eps) then
            return
         end if

      end do
   end subroutine tiling_sor

   subroutine subtiling_sor(u, u_size, dx, dy)
      implicit none
      integer, intent(in) :: u_size
      real*8, intent(inout) :: u(u_size*u_size)
      real*8, intent(in) :: dx, dy

      integer :: tile_count, iter, i, l0, l1, l2, l3, l4
      real*8 :: invdx2, invdy2, f0, f1, u_old, error

      tile_count = (u_size - 2)/tile_size

      invdx2 = 1.0d0/(dx*dx)
      invdy2 = 1.0d0/(dy*dy)

      f0 = 1.0d0 - omega
      f1 = omega / (2.0d0*invdx2 + 2.0d0*invdy2)

      do iter = 1, max_iter

         error = 0.0d0
         i = u_size + 2

         do l2 = 0, subtile_level
            do l1 = 1, tile_size - l2
               do l0 = 1, tile_size - l2

                  u_old = u(i)
                  u(i) = f0*u_old + &
                     f1*((u(i-1) + u(i+1))*invdx2 + &
                     (u(i-u_size) + u(i+u_size))*invdy2)

                  if (l2.eq.subtile_level) then
                     error = error + abs(u(i) - u_old)
                  end if
                  i = i + 1

               end do
               i = i + u_size - l0 + 1
            end do
            i = i - u_size*(l1 - 1)
         end do
         i = i + tile_size

         do l3 = 2, tile_count - 1
            do l2 = 0, subtile_level
               do l1 = 1, tile_size - l2
                  do l0 = 1, tile_size

                     u_old = u(i)
                     u(i) = f0*u_old + &
                        f1*((u(i-1) + u(i+1))*invdx2 + &
                        (u(i-u_size) + u(i+u_size))*invdy2)

                     if (l2.eq.subtile_level) then
                        error = error + abs(u(i) - u_old)
                     end if
                     i = i + 1

                  end do
                  i = i + u_size - l0 + 1
               end do
               i = i - u_size*(l1 - 1) - 1
            end do
            i = i + tile_size + subtile_level + 1
         end do

         do l2 = 0, subtile_level
            do l1 = 1, tile_size - l2
               do l0 = 1, tile_size + l2

                  u_old = u(i)
                  u(i) = f0*u_old + &
                     f1*((u(i-1) + u(i+1))*invdx2 + &
                     (u(i-u_size) + u(i+u_size))*invdy2)

                  if (l2.eq.subtile_level) then
                     error = error + abs(u(i) - u_old)
                  end if
                  i = i + 1

               end do
               i = i + u_size - l0 + 1
            end do
            i = i - u_size*(l1 - 1) - 1
         end do
         i = i + tile_size + subtile_level + u_size*(tile_size - 1) + 3

         do l4 = 2, tile_count - 1
            do l2 = 0, subtile_level
               do l1 = 1, tile_size
                  do l0 = 1, tile_size - l2

                     u_old = u(i)
                     u(i) = f0*u_old + &
                        f1*((u(i-1) + u(i+1))*invdx2 + &
                        (u(i-u_size) + u(i+u_size))*invdy2)

                     if (l2.eq.subtile_level) then
                        error = error + abs(u(i) - u_old)
                     end if
                     i = i + 1

                  end do
                  i = i + u_size - l0 + 1
               end do
               i = i - u_size*l1
            end do
            i = i + u_size*(subtile_level + 1) + tile_size

            do l3 = 2, tile_count - 1
               do l2 = 0, subtile_level
                  do l1 = 1, tile_size
                     do l0 = 1, tile_size

                        u_old = u(i)
                        u(i) = f0*u_old + &
                           f1*((u(i-1) + u(i+1))*invdx2 + &
                           (u(i-u_size) + u(i+u_size))*invdy2)

                        if (l2.eq.subtile_level) then
                           error = error + abs(u(i) - u_old)
                        end if
                        i = i + 1

                     end do
                     i = i + u_size - l0 + 1
                  end do
                  i = i - u_size*l1 - 1
               end do
               i = i + u_size*(subtile_level + 1) + tile_size + subtile_level + 1
            end do

            do l2 = 0, subtile_level
               do l1 = 1, tile_size
                  do l0 = 1, tile_size + l2

                     u_old = u(i)
                     u(i) = f0*u_old + &
                        f1*((u(i-1) + u(i+1))*invdx2 + &
                        (u(i-u_size) + u(i+u_size))*invdy2)

                     if (l2.eq.subtile_level) then
                        error = error + abs(u(i) - u_old)
                     end if
                     i = i + 1

                  end do
                  i = i + u_size - l0 + 1
               end do
               i = i - u_size*l1 - 1
            end do
            i = i + u_size*(subtile_level + 1) + tile_size + subtile_level + u_size*(tile_size - 1) + 3
         end do

         do l2 = 0, subtile_level
            do l1 = 1, tile_size + l2
               do l0 = 1, tile_size - l2

                  u_old = u(i)
                  u(i) = f0*u_old + &
                     f1*((u(i-1) + u(i+1))*invdx2 + &
                     (u(i-u_size) + u(i+u_size))*invdy2)

                  if (l2.eq.subtile_level) then
                     error = error + abs(u(i) - u_old)
                  end if
                  i = i + 1

               end do
               i = i + u_size - l0 + 1
            end do
            i = i - u_size*l1
         end do
         i = i + u_size*(subtile_level + 1) + tile_size

         do l3 = 2, tile_count - 1
            do l2 = 0, subtile_level
               do l1 = 1, tile_size + l2
                  do l0 = 1, tile_size

                     u_old = u(i)
                     u(i) = f0*u_old + &
                        f1*((u(i-1) + u(i+1))*invdx2 + &
                        (u(i-u_size) + u(i+u_size))*invdy2)

                     if (l2.eq.subtile_level) then
                        error = error + abs(u(i) - u_old)
                     end if
                     i = i + 1

                  end do
                  i = i + u_size - l0 + 1
               end do
               i = i - u_size*l1 - 1
            end do
            i = i + u_size*(subtile_level + 1) + tile_size + subtile_level + 1
         end do

         do l2 = 0, subtile_level
            do l1 = 1, tile_size + l2
               do l0 = 1, tile_size + l2

                  u_old = u(i)
                  u(i) = f0*u_old + &
                     f1*((u(i-1) + u(i+1))*invdx2 + &
                     (u(i-u_size) + u(i+u_size))*invdy2)

                  if (l2.eq.subtile_level) then
                     error = error + abs(u(i) - u_old)
                  end if
                  i = i + 1

               end do
               i = i + u_size - l0 + 1
            end do
            i = i - u_size*l1 - 1
         end do

         if (error < eps) then
            return
         end if

      end do
   end subroutine subtiling_sor

   subroutine subtiling_sor_test_version(u, u_size, dx, dy)
      implicit none
      integer, intent(in) :: u_size
      real*8, intent(inout) :: u(u_size*u_size)
      real*8, intent(in) :: dx, dy

      integer, dimension(:,:), allocatable :: subtile_indices
      integer :: tile_count, subtile_count, subtile_idx, iter, i, l0, l1, l2, l3, l4, s0, s1, s2, s3
      real*8 :: invdx2, invdy2, f0, f1, u_old, error

      tile_count = (u_size - 2)/tile_size
      subtile_count = tile_count*tile_count*(subtile_level+1)
      allocate(subtile_indices(subtile_count, 4))
      subtile_idx = 1

      do l4 = 1, tile_count

         if (l4.eq.1) then
            s1 = -1
         else if (l4.eq.tile_count) then
            s1 = 1
         else
            s1 = 0
         end if

         if (l4.eq.1) then
            s3 = 0
         else
            s3 = 1
         end if

         do l3 = 1, tile_count

            if (l3.eq.1) then
               s0 = -1
            else if (l3.eq.tile_count) then
               s0 = 1
            else
               s0 = 0
            end if

            if (l3.eq.1) then
               s2 = 0
            else
               s2 = 1
            end if

            do l2 = 0, subtile_level
               subtile_indices(subtile_idx, 1) = (l4-1)*tile_size*u_size + (l3-1)*tile_size + u_size + 2 - l2*s3*u_size - l2*s2
               subtile_indices(subtile_idx, 2) = tile_size + l2*s1
               subtile_indices(subtile_idx, 3) = tile_size + l2*s0
               if (l2.eq.subtile_level) then
                  subtile_indices(subtile_idx, 4) = 1
               else
                  subtile_indices(subtile_idx, 4) = 0
               end if
               subtile_idx = subtile_idx + 1

            end do
         end do
      end do

      invdx2 = 1.0d0/(dx*dx)
      invdy2 = 1.0d0/(dy*dy)

      f0 = 1.0d0 - omega
      f1 = omega / (2.0d0*invdx2 + 2.0d0*invdy2)

      do iter = 1, max_iter

         error = 0.0d0
         i = u_size + 2

         do subtile_idx = 1, subtile_count
            i = subtile_indices(subtile_idx, 1)
            l2 = subtile_indices(subtile_idx, 2)
            l3 = subtile_indices(subtile_idx, 3)
            l4 = subtile_indices(subtile_idx, 4)
            if (l4.eq.1) then
               do l1 = 1, l2
                  do l0 = 1, l3

                     u_old = u(i)
                     u(i) = f0*u_old + &
                        f1*((u(i-1) + u(i+1))*invdx2 + &
                        (u(i-u_size) + u(i+u_size))*invdy2)
                     error = error + abs(u(i) - u_old)
                     i = i + 1

                  end do
                  i = i + u_size - l3
               end do
            else
               do l1 = 1, l2
                  do l0 = 1, l3

                     u_old = u(i)
                     u(i) = f0*u_old + &
                        f1*((u(i-1) + u(i+1))*invdx2 + &
                        (u(i-u_size) + u(i+u_size))*invdy2)
                     i = i + 1

                  end do
                  i = i + u_size - l3
               end do
            end if
         end do

         if (error < eps) then
            exit
         end if

      end do

      deallocate(subtile_indices)

   end subroutine subtiling_sor_test_version

   subroutine subtiling_sor_4(u, u_size, dx, dy)
      implicit none
      integer, intent(in) :: u_size
      real*8, intent(inout) :: u(u_size*u_size)
      real*8, intent(in) :: dx, dy

      integer, parameter :: subtile_level_param = 4
      integer :: tile_size_param = 4
      integer :: tile_count, iter, i, l0, l1, l2, l3, l4
      real*8 :: invdx2, invdy2, f0, f1, u_old, error
      tile_count = (u_size - 2)/tile_size_param
      invdx2 = 1.0d0/(dx*dx)
      invdy2 = 1.0d0/(dy*dy)
      f0 = 1.0d0 - omega
      f1 = omega / (2.0d0*invdx2 + 2.0d0*invdy2)
      do iter = 1, max_iter
         error = 0.0d0
         i = u_size + 2
         do l1 = 1, tile_size_param - 0
            do l0 = 1, tile_size_param - 0
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - u_size*(l1 - 1)
         do l1 = 1, tile_size_param - 1
            do l0 = 1, tile_size_param - 1
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - u_size*(l1 - 1)
         do l1 = 1, tile_size_param - 2
            do l0 = 1, tile_size_param - 2
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - u_size*(l1 - 1)
         do l1 = 1, tile_size_param - 3
            do l0 = 1, tile_size_param - 3
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - u_size*(l1 - 1)
         do l1 = 1, tile_size_param - 4
            do l0 = 1, tile_size_param - 4
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               error = error + Abs(u_old - u(i))
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i + tile_size_param - u_size*(l1 - 1)
         do l3 = 2, tile_count - 1
            do l1 = 1, tile_size_param - 0
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - u_size*(l1 - 1) - 1
            do l1 = 1, tile_size_param - 1
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - u_size*(l1 - 1) - 1
            do l1 = 1, tile_size_param - 2
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - u_size*(l1 - 1) - 1
            do l1 = 1, tile_size_param - 3
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - u_size*(l1 - 1) - 1
            do l1 = 1, tile_size_param - 4
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  error = error + Abs(u_old - u(i))
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i + subtile_level_param + tile_size_param - u_size*(l1 - 1)
         end do
         do l1 = 1, tile_size_param - 0
            do l0 = 1, tile_size_param + 0
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - u_size*(l1 - 1) - 1
         do l1 = 1, tile_size_param - 1
            do l0 = 1, tile_size_param + 1
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - u_size*(l1 - 1) - 1
         do l1 = 1, tile_size_param - 2
            do l0 = 1, tile_size_param + 2
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - u_size*(l1 - 1) - 1
         do l1 = 1, tile_size_param - 3
            do l0 = 1, tile_size_param + 3
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - u_size*(l1 - 1) - 1
         do l1 = 1, tile_size_param - 4
            do l0 = 1, tile_size_param + 4
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               error = error + Abs(u_old - u(i))
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - l1*u_size + subtile_level_param + tile_size_param*u_size + tile_size_param + 2
         do l4 = 2, tile_count - 1
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param - 0
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - l1*u_size
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param - 1
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - l1*u_size
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param - 2
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - l1*u_size
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param - 3
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - l1*u_size
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param - 4
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  error = error + Abs(u_old - u(i))
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - l1*u_size + tile_size_param + u_size*(subtile_level_param + 1)
            do l3 = 2, tile_count - 1
               do l1 = 1, tile_size_param
                  do l0 = 1, tile_size_param
                     u_old = u(i)
                     u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                     i = i + 1
                  end do
                  i = i - l0 + u_size + 1
               end do
               i = i - l1*u_size - 1
               do l1 = 1, tile_size_param
                  do l0 = 1, tile_size_param
                     u_old = u(i)
                     u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                     i = i + 1
                  end do
                  i = i - l0 + u_size + 1
               end do
               i = i - l1*u_size - 1
               do l1 = 1, tile_size_param
                  do l0 = 1, tile_size_param
                     u_old = u(i)
                     u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                     i = i + 1
                  end do
                  i = i - l0 + u_size + 1
               end do
               i = i - l1*u_size - 1
               do l1 = 1, tile_size_param
                  do l0 = 1, tile_size_param
                     u_old = u(i)
                     u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                     i = i + 1
                  end do
                  i = i - l0 + u_size + 1
               end do
               i = i - l1*u_size - 1
               do l1 = 1, tile_size_param
                  do l0 = 1, tile_size_param
                     u_old = u(i)
                     u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                     error = error + Abs(u_old - u(i))
                     i = i + 1
                  end do
                  i = i - l0 + u_size + 1
               end do
               i = i - l1*u_size + subtile_level_param + tile_size_param + u_size*(subtile_level_param + 1)
            end do
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param + 0
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - l1*u_size - 1
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param + 1
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - l1*u_size - 1
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param + 2
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - l1*u_size - 1
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param + 3
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - l1*u_size - 1
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param + 4
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  error = error + Abs(u_old - u(i))
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - l1*u_size + subtile_level_param*u_size + subtile_level_param + tile_size_param*u_size + tile_size_param + 2
         end do
         do l1 = 1, tile_size_param + 0
            do l0 = 1, tile_size_param - 0
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - l1*u_size
         do l1 = 1, tile_size_param + 1
            do l0 = 1, tile_size_param - 1
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - l1*u_size
         do l1 = 1, tile_size_param + 2
            do l0 = 1, tile_size_param - 2
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - l1*u_size
         do l1 = 1, tile_size_param + 3
            do l0 = 1, tile_size_param - 3
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - l1*u_size
         do l1 = 1, tile_size_param + 4
            do l0 = 1, tile_size_param - 4
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               error = error + Abs(u_old - u(i))
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - l1*u_size + tile_size_param + u_size*(subtile_level_param + 1)
         do l3 = 2, tile_count - 1
            do l1 = 1, tile_size_param + 0
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - l1*u_size - 1
            do l1 = 1, tile_size_param + 1
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - l1*u_size - 1
            do l1 = 1, tile_size_param + 2
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - l1*u_size - 1
            do l1 = 1, tile_size_param + 3
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - l1*u_size - 1
            do l1 = 1, tile_size_param + 4
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  error = error + Abs(u_old - u(i))
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - l1*u_size + subtile_level_param + tile_size_param + u_size*(subtile_level_param + 1)
         end do
         do l1 = 1, tile_size_param + 0
            do l0 = 1, tile_size_param + 0
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - l1*u_size - 1
         do l1 = 1, tile_size_param + 1
            do l0 = 1, tile_size_param + 1
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - l1*u_size - 1
         do l1 = 1, tile_size_param + 2
            do l0 = 1, tile_size_param + 2
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - l1*u_size - 1
         do l1 = 1, tile_size_param + 3
            do l0 = 1, tile_size_param + 3
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - l1*u_size - 1
         do l1 = 1, tile_size_param + 4
            do l0 = 1, tile_size_param + 4
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               error = error + Abs(u_old - u(i))
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - l1*u_size - 1
         if (error < eps) then
            return
         end if
      end do
   end subroutine subtiling_sor_4

   subroutine subtiling_sor_8(u, u_size, dx, dy)
      implicit none
      integer, intent(in) :: u_size
      real*8, intent(inout) :: u(u_size*u_size)
      real*8, intent(in) :: dx, dy

      integer, parameter :: subtile_level_param = 8
      integer :: tile_size_param = 8
      integer :: tile_count, iter, i, l0, l1, l2, l3, l4
      real*8 :: invdx2, invdy2, f0, f1, u_old, error
      tile_count = (u_size - 2)/tile_size_param
      invdx2 = 1.0d0/(dx*dx)
      invdy2 = 1.0d0/(dy*dy)
      f0 = 1.0d0 - omega
      f1 = omega / (2.0d0*invdx2 + 2.0d0*invdy2)
      do iter = 1, max_iter
         error = 0.0d0
         i = u_size + 2
         do l1 = 1, tile_size_param - 0
            do l0 = 1, tile_size_param - 0
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - u_size*(l1 - 1)
         do l1 = 1, tile_size_param - 1
            do l0 = 1, tile_size_param - 1
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - u_size*(l1 - 1)
         do l1 = 1, tile_size_param - 2
            do l0 = 1, tile_size_param - 2
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - u_size*(l1 - 1)
         do l1 = 1, tile_size_param - 3
            do l0 = 1, tile_size_param - 3
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - u_size*(l1 - 1)
         do l1 = 1, tile_size_param - 4
            do l0 = 1, tile_size_param - 4
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - u_size*(l1 - 1)
         do l1 = 1, tile_size_param - 5
            do l0 = 1, tile_size_param - 5
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - u_size*(l1 - 1)
         do l1 = 1, tile_size_param - 6
            do l0 = 1, tile_size_param - 6
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - u_size*(l1 - 1)
         do l1 = 1, tile_size_param - 7
            do l0 = 1, tile_size_param - 7
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - u_size*(l1 - 1)
         do l1 = 1, tile_size_param - 8
            do l0 = 1, tile_size_param - 8
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               error = error + Abs(u_old - u(i))
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i + tile_size_param - u_size*(l1 - 1)
         do l3 = 2, tile_count - 1
            do l1 = 1, tile_size_param - 0
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - u_size*(l1 - 1) - 1
            do l1 = 1, tile_size_param - 1
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - u_size*(l1 - 1) - 1
            do l1 = 1, tile_size_param - 2
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - u_size*(l1 - 1) - 1
            do l1 = 1, tile_size_param - 3
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - u_size*(l1 - 1) - 1
            do l1 = 1, tile_size_param - 4
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - u_size*(l1 - 1) - 1
            do l1 = 1, tile_size_param - 5
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - u_size*(l1 - 1) - 1
            do l1 = 1, tile_size_param - 6
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - u_size*(l1 - 1) - 1
            do l1 = 1, tile_size_param - 7
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - u_size*(l1 - 1) - 1
            do l1 = 1, tile_size_param - 8
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  error = error + Abs(u_old - u(i))
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i + subtile_level_param + tile_size_param - u_size*(l1 - 1)
         end do
         do l1 = 1, tile_size_param - 0
            do l0 = 1, tile_size_param + 0
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - u_size*(l1 - 1) - 1
         do l1 = 1, tile_size_param - 1
            do l0 = 1, tile_size_param + 1
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - u_size*(l1 - 1) - 1
         do l1 = 1, tile_size_param - 2
            do l0 = 1, tile_size_param + 2
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - u_size*(l1 - 1) - 1
         do l1 = 1, tile_size_param - 3
            do l0 = 1, tile_size_param + 3
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - u_size*(l1 - 1) - 1
         do l1 = 1, tile_size_param - 4
            do l0 = 1, tile_size_param + 4
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - u_size*(l1 - 1) - 1
         do l1 = 1, tile_size_param - 5
            do l0 = 1, tile_size_param + 5
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - u_size*(l1 - 1) - 1
         do l1 = 1, tile_size_param - 6
            do l0 = 1, tile_size_param + 6
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - u_size*(l1 - 1) - 1
         do l1 = 1, tile_size_param - 7
            do l0 = 1, tile_size_param + 7
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - u_size*(l1 - 1) - 1
         do l1 = 1, tile_size_param - 8
            do l0 = 1, tile_size_param + 8
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               error = error + Abs(u_old - u(i))
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - l1*u_size + subtile_level_param + tile_size_param*u_size + tile_size_param + 2
         do l4 = 2, tile_count - 1
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param - 0
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - l1*u_size
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param - 1
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - l1*u_size
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param - 2
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - l1*u_size
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param - 3
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - l1*u_size
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param - 4
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - l1*u_size
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param - 5
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - l1*u_size
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param - 6
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - l1*u_size
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param - 7
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - l1*u_size
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param - 8
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  error = error + Abs(u_old - u(i))
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - l1*u_size + tile_size_param + u_size*(subtile_level_param + 1)
            do l3 = 2, tile_count - 1
               do l1 = 1, tile_size_param
                  do l0 = 1, tile_size_param
                     u_old = u(i)
                     u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                     i = i + 1
                  end do
                  i = i - l0 + u_size + 1
               end do
               i = i - l1*u_size - 1
               do l1 = 1, tile_size_param
                  do l0 = 1, tile_size_param
                     u_old = u(i)
                     u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                     i = i + 1
                  end do
                  i = i - l0 + u_size + 1
               end do
               i = i - l1*u_size - 1
               do l1 = 1, tile_size_param
                  do l0 = 1, tile_size_param
                     u_old = u(i)
                     u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                     i = i + 1
                  end do
                  i = i - l0 + u_size + 1
               end do
               i = i - l1*u_size - 1
               do l1 = 1, tile_size_param
                  do l0 = 1, tile_size_param
                     u_old = u(i)
                     u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                     i = i + 1
                  end do
                  i = i - l0 + u_size + 1
               end do
               i = i - l1*u_size - 1
               do l1 = 1, tile_size_param
                  do l0 = 1, tile_size_param
                     u_old = u(i)
                     u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                     i = i + 1
                  end do
                  i = i - l0 + u_size + 1
               end do
               i = i - l1*u_size - 1
               do l1 = 1, tile_size_param
                  do l0 = 1, tile_size_param
                     u_old = u(i)
                     u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                     i = i + 1
                  end do
                  i = i - l0 + u_size + 1
               end do
               i = i - l1*u_size - 1
               do l1 = 1, tile_size_param
                  do l0 = 1, tile_size_param
                     u_old = u(i)
                     u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                     i = i + 1
                  end do
                  i = i - l0 + u_size + 1
               end do
               i = i - l1*u_size - 1
               do l1 = 1, tile_size_param
                  do l0 = 1, tile_size_param
                     u_old = u(i)
                     u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                     i = i + 1
                  end do
                  i = i - l0 + u_size + 1
               end do
               i = i - l1*u_size - 1
               do l1 = 1, tile_size_param
                  do l0 = 1, tile_size_param
                     u_old = u(i)
                     u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                     error = error + Abs(u_old - u(i))
                     i = i + 1
                  end do
                  i = i - l0 + u_size + 1
               end do
               i = i - l1*u_size + subtile_level_param + tile_size_param + u_size*(subtile_level_param + 1)
            end do
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param + 0
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - l1*u_size - 1
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param + 1
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - l1*u_size - 1
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param + 2
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - l1*u_size - 1
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param + 3
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - l1*u_size - 1
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param + 4
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - l1*u_size - 1
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param + 5
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - l1*u_size - 1
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param + 6
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - l1*u_size - 1
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param + 7
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - l1*u_size - 1
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param + 8
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  error = error + Abs(u_old - u(i))
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - l1*u_size + subtile_level_param*u_size + subtile_level_param + tile_size_param*u_size + tile_size_param + 2
         end do
         do l1 = 1, tile_size_param + 0
            do l0 = 1, tile_size_param - 0
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - l1*u_size
         do l1 = 1, tile_size_param + 1
            do l0 = 1, tile_size_param - 1
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - l1*u_size
         do l1 = 1, tile_size_param + 2
            do l0 = 1, tile_size_param - 2
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - l1*u_size
         do l1 = 1, tile_size_param + 3
            do l0 = 1, tile_size_param - 3
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - l1*u_size
         do l1 = 1, tile_size_param + 4
            do l0 = 1, tile_size_param - 4
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - l1*u_size
         do l1 = 1, tile_size_param + 5
            do l0 = 1, tile_size_param - 5
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - l1*u_size
         do l1 = 1, tile_size_param + 6
            do l0 = 1, tile_size_param - 6
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - l1*u_size
         do l1 = 1, tile_size_param + 7
            do l0 = 1, tile_size_param - 7
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - l1*u_size
         do l1 = 1, tile_size_param + 8
            do l0 = 1, tile_size_param - 8
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               error = error + Abs(u_old - u(i))
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - l1*u_size + tile_size_param + u_size*(subtile_level_param + 1)
         do l3 = 2, tile_count - 1
            do l1 = 1, tile_size_param + 0
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - l1*u_size - 1
            do l1 = 1, tile_size_param + 1
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - l1*u_size - 1
            do l1 = 1, tile_size_param + 2
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - l1*u_size - 1
            do l1 = 1, tile_size_param + 3
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - l1*u_size - 1
            do l1 = 1, tile_size_param + 4
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - l1*u_size - 1
            do l1 = 1, tile_size_param + 5
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - l1*u_size - 1
            do l1 = 1, tile_size_param + 6
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - l1*u_size - 1
            do l1 = 1, tile_size_param + 7
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - l1*u_size - 1
            do l1 = 1, tile_size_param + 8
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  error = error + Abs(u_old - u(i))
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - l1*u_size + subtile_level_param + tile_size_param + u_size*(subtile_level_param + 1)
         end do
         do l1 = 1, tile_size_param + 0
            do l0 = 1, tile_size_param + 0
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - l1*u_size - 1
         do l1 = 1, tile_size_param + 1
            do l0 = 1, tile_size_param + 1
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - l1*u_size - 1
         do l1 = 1, tile_size_param + 2
            do l0 = 1, tile_size_param + 2
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - l1*u_size - 1
         do l1 = 1, tile_size_param + 3
            do l0 = 1, tile_size_param + 3
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - l1*u_size - 1
         do l1 = 1, tile_size_param + 4
            do l0 = 1, tile_size_param + 4
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - l1*u_size - 1
         do l1 = 1, tile_size_param + 5
            do l0 = 1, tile_size_param + 5
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - l1*u_size - 1
         do l1 = 1, tile_size_param + 6
            do l0 = 1, tile_size_param + 6
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - l1*u_size - 1
         do l1 = 1, tile_size_param + 7
            do l0 = 1, tile_size_param + 7
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - l1*u_size - 1
         do l1 = 1, tile_size_param + 8
            do l0 = 1, tile_size_param + 8
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               error = error + Abs(u_old - u(i))
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - l1*u_size - 1
         if (error < eps) then
            return
         end if
      end do
   end subroutine subtiling_sor_8

   subroutine subtiling_sor_16(u, u_size, dx, dy)
      implicit none
      integer, intent(in) :: u_size
      real*8, intent(inout) :: u(u_size*u_size)
      real*8, intent(in) :: dx, dy

      integer, parameter :: subtile_level_param = 16
      integer :: tile_size_param = 16
      integer :: tile_count, iter, i, l0, l1, l2, l3, l4
      real*8 :: invdx2, invdy2, f0, f1, u_old, error
      tile_count = (u_size - 2)/tile_size_param
      invdx2 = 1.0d0/(dx*dx)
      invdy2 = 1.0d0/(dy*dy)
      f0 = 1.0d0 - omega
      f1 = omega / (2.0d0*invdx2 + 2.0d0*invdy2)
      do iter = 1, max_iter
         error = 0.0d0
         i = u_size + 2
         do l1 = 1, tile_size_param - 0
            do l0 = 1, tile_size_param - 0
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - u_size*(l1 - 1)
         do l1 = 1, tile_size_param - 1
            do l0 = 1, tile_size_param - 1
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - u_size*(l1 - 1)
         do l1 = 1, tile_size_param - 2
            do l0 = 1, tile_size_param - 2
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - u_size*(l1 - 1)
         do l1 = 1, tile_size_param - 3
            do l0 = 1, tile_size_param - 3
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - u_size*(l1 - 1)
         do l1 = 1, tile_size_param - 4
            do l0 = 1, tile_size_param - 4
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - u_size*(l1 - 1)
         do l1 = 1, tile_size_param - 5
            do l0 = 1, tile_size_param - 5
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - u_size*(l1 - 1)
         do l1 = 1, tile_size_param - 6
            do l0 = 1, tile_size_param - 6
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - u_size*(l1 - 1)
         do l1 = 1, tile_size_param - 7
            do l0 = 1, tile_size_param - 7
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - u_size*(l1 - 1)
         do l1 = 1, tile_size_param - 8
            do l0 = 1, tile_size_param - 8
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - u_size*(l1 - 1)
         do l1 = 1, tile_size_param - 9
            do l0 = 1, tile_size_param - 9
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - u_size*(l1 - 1)
         do l1 = 1, tile_size_param - 10
            do l0 = 1, tile_size_param - 10
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - u_size*(l1 - 1)
         do l1 = 1, tile_size_param - 11
            do l0 = 1, tile_size_param - 11
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - u_size*(l1 - 1)
         do l1 = 1, tile_size_param - 12
            do l0 = 1, tile_size_param - 12
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - u_size*(l1 - 1)
         do l1 = 1, tile_size_param - 13
            do l0 = 1, tile_size_param - 13
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - u_size*(l1 - 1)
         do l1 = 1, tile_size_param - 14
            do l0 = 1, tile_size_param - 14
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - u_size*(l1 - 1)
         do l1 = 1, tile_size_param - 15
            do l0 = 1, tile_size_param - 15
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - u_size*(l1 - 1)
         do l1 = 1, tile_size_param - 16
            do l0 = 1, tile_size_param - 16
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               error = error + Abs(u_old - u(i))
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i + tile_size_param - u_size*(l1 - 1)
         do l3 = 2, tile_count - 1
            do l1 = 1, tile_size_param - 0
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - u_size*(l1 - 1) - 1
            do l1 = 1, tile_size_param - 1
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - u_size*(l1 - 1) - 1
            do l1 = 1, tile_size_param - 2
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - u_size*(l1 - 1) - 1
            do l1 = 1, tile_size_param - 3
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - u_size*(l1 - 1) - 1
            do l1 = 1, tile_size_param - 4
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - u_size*(l1 - 1) - 1
            do l1 = 1, tile_size_param - 5
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - u_size*(l1 - 1) - 1
            do l1 = 1, tile_size_param - 6
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - u_size*(l1 - 1) - 1
            do l1 = 1, tile_size_param - 7
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - u_size*(l1 - 1) - 1
            do l1 = 1, tile_size_param - 8
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - u_size*(l1 - 1) - 1
            do l1 = 1, tile_size_param - 9
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - u_size*(l1 - 1) - 1
            do l1 = 1, tile_size_param - 10
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - u_size*(l1 - 1) - 1
            do l1 = 1, tile_size_param - 11
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - u_size*(l1 - 1) - 1
            do l1 = 1, tile_size_param - 12
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - u_size*(l1 - 1) - 1
            do l1 = 1, tile_size_param - 13
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - u_size*(l1 - 1) - 1
            do l1 = 1, tile_size_param - 14
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - u_size*(l1 - 1) - 1
            do l1 = 1, tile_size_param - 15
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - u_size*(l1 - 1) - 1
            do l1 = 1, tile_size_param - 16
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  error = error + Abs(u_old - u(i))
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i + subtile_level_param + tile_size_param - u_size*(l1 - 1)
         end do
         do l1 = 1, tile_size_param - 0
            do l0 = 1, tile_size_param + 0
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - u_size*(l1 - 1) - 1
         do l1 = 1, tile_size_param - 1
            do l0 = 1, tile_size_param + 1
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - u_size*(l1 - 1) - 1
         do l1 = 1, tile_size_param - 2
            do l0 = 1, tile_size_param + 2
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - u_size*(l1 - 1) - 1
         do l1 = 1, tile_size_param - 3
            do l0 = 1, tile_size_param + 3
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - u_size*(l1 - 1) - 1
         do l1 = 1, tile_size_param - 4
            do l0 = 1, tile_size_param + 4
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - u_size*(l1 - 1) - 1
         do l1 = 1, tile_size_param - 5
            do l0 = 1, tile_size_param + 5
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - u_size*(l1 - 1) - 1
         do l1 = 1, tile_size_param - 6
            do l0 = 1, tile_size_param + 6
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - u_size*(l1 - 1) - 1
         do l1 = 1, tile_size_param - 7
            do l0 = 1, tile_size_param + 7
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - u_size*(l1 - 1) - 1
         do l1 = 1, tile_size_param - 8
            do l0 = 1, tile_size_param + 8
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - u_size*(l1 - 1) - 1
         do l1 = 1, tile_size_param - 9
            do l0 = 1, tile_size_param + 9
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - u_size*(l1 - 1) - 1
         do l1 = 1, tile_size_param - 10
            do l0 = 1, tile_size_param + 10
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - u_size*(l1 - 1) - 1
         do l1 = 1, tile_size_param - 11
            do l0 = 1, tile_size_param + 11
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - u_size*(l1 - 1) - 1
         do l1 = 1, tile_size_param - 12
            do l0 = 1, tile_size_param + 12
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - u_size*(l1 - 1) - 1
         do l1 = 1, tile_size_param - 13
            do l0 = 1, tile_size_param + 13
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - u_size*(l1 - 1) - 1
         do l1 = 1, tile_size_param - 14
            do l0 = 1, tile_size_param + 14
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - u_size*(l1 - 1) - 1
         do l1 = 1, tile_size_param - 15
            do l0 = 1, tile_size_param + 15
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - u_size*(l1 - 1) - 1
         do l1 = 1, tile_size_param - 16
            do l0 = 1, tile_size_param + 16
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               error = error + Abs(u_old - u(i))
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - l1*u_size + subtile_level_param + tile_size_param*u_size + tile_size_param + 2
         do l4 = 2, tile_count - 1
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param - 0
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - l1*u_size
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param - 1
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - l1*u_size
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param - 2
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - l1*u_size
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param - 3
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - l1*u_size
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param - 4
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - l1*u_size
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param - 5
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - l1*u_size
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param - 6
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - l1*u_size
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param - 7
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - l1*u_size
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param - 8
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - l1*u_size
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param - 9
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - l1*u_size
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param - 10
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - l1*u_size
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param - 11
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - l1*u_size
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param - 12
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - l1*u_size
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param - 13
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - l1*u_size
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param - 14
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - l1*u_size
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param - 15
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - l1*u_size
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param - 16
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  error = error + Abs(u_old - u(i))
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - l1*u_size + tile_size_param + u_size*(subtile_level_param + 1)
            do l3 = 2, tile_count - 1
               do l1 = 1, tile_size_param
                  do l0 = 1, tile_size_param
                     u_old = u(i)
                     u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                     i = i + 1
                  end do
                  i = i - l0 + u_size + 1
               end do
               i = i - l1*u_size - 1
               do l1 = 1, tile_size_param
                  do l0 = 1, tile_size_param
                     u_old = u(i)
                     u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                     i = i + 1
                  end do
                  i = i - l0 + u_size + 1
               end do
               i = i - l1*u_size - 1
               do l1 = 1, tile_size_param
                  do l0 = 1, tile_size_param
                     u_old = u(i)
                     u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                     i = i + 1
                  end do
                  i = i - l0 + u_size + 1
               end do
               i = i - l1*u_size - 1
               do l1 = 1, tile_size_param
                  do l0 = 1, tile_size_param
                     u_old = u(i)
                     u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                     i = i + 1
                  end do
                  i = i - l0 + u_size + 1
               end do
               i = i - l1*u_size - 1
               do l1 = 1, tile_size_param
                  do l0 = 1, tile_size_param
                     u_old = u(i)
                     u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                     i = i + 1
                  end do
                  i = i - l0 + u_size + 1
               end do
               i = i - l1*u_size - 1
               do l1 = 1, tile_size_param
                  do l0 = 1, tile_size_param
                     u_old = u(i)
                     u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                     i = i + 1
                  end do
                  i = i - l0 + u_size + 1
               end do
               i = i - l1*u_size - 1
               do l1 = 1, tile_size_param
                  do l0 = 1, tile_size_param
                     u_old = u(i)
                     u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                     i = i + 1
                  end do
                  i = i - l0 + u_size + 1
               end do
               i = i - l1*u_size - 1
               do l1 = 1, tile_size_param
                  do l0 = 1, tile_size_param
                     u_old = u(i)
                     u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                     i = i + 1
                  end do
                  i = i - l0 + u_size + 1
               end do
               i = i - l1*u_size - 1
               do l1 = 1, tile_size_param
                  do l0 = 1, tile_size_param
                     u_old = u(i)
                     u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                     i = i + 1
                  end do
                  i = i - l0 + u_size + 1
               end do
               i = i - l1*u_size - 1
               do l1 = 1, tile_size_param
                  do l0 = 1, tile_size_param
                     u_old = u(i)
                     u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                     i = i + 1
                  end do
                  i = i - l0 + u_size + 1
               end do
               i = i - l1*u_size - 1
               do l1 = 1, tile_size_param
                  do l0 = 1, tile_size_param
                     u_old = u(i)
                     u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                     i = i + 1
                  end do
                  i = i - l0 + u_size + 1
               end do
               i = i - l1*u_size - 1
               do l1 = 1, tile_size_param
                  do l0 = 1, tile_size_param
                     u_old = u(i)
                     u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                     i = i + 1
                  end do
                  i = i - l0 + u_size + 1
               end do
               i = i - l1*u_size - 1
               do l1 = 1, tile_size_param
                  do l0 = 1, tile_size_param
                     u_old = u(i)
                     u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                     i = i + 1
                  end do
                  i = i - l0 + u_size + 1
               end do
               i = i - l1*u_size - 1
               do l1 = 1, tile_size_param
                  do l0 = 1, tile_size_param
                     u_old = u(i)
                     u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                     i = i + 1
                  end do
                  i = i - l0 + u_size + 1
               end do
               i = i - l1*u_size - 1
               do l1 = 1, tile_size_param
                  do l0 = 1, tile_size_param
                     u_old = u(i)
                     u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                     i = i + 1
                  end do
                  i = i - l0 + u_size + 1
               end do
               i = i - l1*u_size - 1
               do l1 = 1, tile_size_param
                  do l0 = 1, tile_size_param
                     u_old = u(i)
                     u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                     i = i + 1
                  end do
                  i = i - l0 + u_size + 1
               end do
               i = i - l1*u_size - 1
               do l1 = 1, tile_size_param
                  do l0 = 1, tile_size_param
                     u_old = u(i)
                     u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                     error = error + Abs(u_old - u(i))
                     i = i + 1
                  end do
                  i = i - l0 + u_size + 1
               end do
               i = i - l1*u_size + subtile_level_param + tile_size_param + u_size*(subtile_level_param + 1)
            end do
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param + 0
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - l1*u_size - 1
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param + 1
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - l1*u_size - 1
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param + 2
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - l1*u_size - 1
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param + 3
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - l1*u_size - 1
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param + 4
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - l1*u_size - 1
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param + 5
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - l1*u_size - 1
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param + 6
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - l1*u_size - 1
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param + 7
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - l1*u_size - 1
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param + 8
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - l1*u_size - 1
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param + 9
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - l1*u_size - 1
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param + 10
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - l1*u_size - 1
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param + 11
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - l1*u_size - 1
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param + 12
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - l1*u_size - 1
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param + 13
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - l1*u_size - 1
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param + 14
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - l1*u_size - 1
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param + 15
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - l1*u_size - 1
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param + 16
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  error = error + Abs(u_old - u(i))
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - l1*u_size + subtile_level_param*u_size + subtile_level_param + tile_size_param*u_size + tile_size_param + 2
         end do
         do l1 = 1, tile_size_param + 0
            do l0 = 1, tile_size_param - 0
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - l1*u_size
         do l1 = 1, tile_size_param + 1
            do l0 = 1, tile_size_param - 1
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - l1*u_size
         do l1 = 1, tile_size_param + 2
            do l0 = 1, tile_size_param - 2
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - l1*u_size
         do l1 = 1, tile_size_param + 3
            do l0 = 1, tile_size_param - 3
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - l1*u_size
         do l1 = 1, tile_size_param + 4
            do l0 = 1, tile_size_param - 4
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - l1*u_size
         do l1 = 1, tile_size_param + 5
            do l0 = 1, tile_size_param - 5
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - l1*u_size
         do l1 = 1, tile_size_param + 6
            do l0 = 1, tile_size_param - 6
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - l1*u_size
         do l1 = 1, tile_size_param + 7
            do l0 = 1, tile_size_param - 7
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - l1*u_size
         do l1 = 1, tile_size_param + 8
            do l0 = 1, tile_size_param - 8
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - l1*u_size
         do l1 = 1, tile_size_param + 9
            do l0 = 1, tile_size_param - 9
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - l1*u_size
         do l1 = 1, tile_size_param + 10
            do l0 = 1, tile_size_param - 10
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - l1*u_size
         do l1 = 1, tile_size_param + 11
            do l0 = 1, tile_size_param - 11
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - l1*u_size
         do l1 = 1, tile_size_param + 12
            do l0 = 1, tile_size_param - 12
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - l1*u_size
         do l1 = 1, tile_size_param + 13
            do l0 = 1, tile_size_param - 13
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - l1*u_size
         do l1 = 1, tile_size_param + 14
            do l0 = 1, tile_size_param - 14
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - l1*u_size
         do l1 = 1, tile_size_param + 15
            do l0 = 1, tile_size_param - 15
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - l1*u_size
         do l1 = 1, tile_size_param + 16
            do l0 = 1, tile_size_param - 16
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               error = error + Abs(u_old - u(i))
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - l1*u_size + tile_size_param + u_size*(subtile_level_param + 1)
         do l3 = 2, tile_count - 1
            do l1 = 1, tile_size_param + 0
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - l1*u_size - 1
            do l1 = 1, tile_size_param + 1
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - l1*u_size - 1
            do l1 = 1, tile_size_param + 2
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - l1*u_size - 1
            do l1 = 1, tile_size_param + 3
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - l1*u_size - 1
            do l1 = 1, tile_size_param + 4
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - l1*u_size - 1
            do l1 = 1, tile_size_param + 5
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - l1*u_size - 1
            do l1 = 1, tile_size_param + 6
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - l1*u_size - 1
            do l1 = 1, tile_size_param + 7
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - l1*u_size - 1
            do l1 = 1, tile_size_param + 8
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - l1*u_size - 1
            do l1 = 1, tile_size_param + 9
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - l1*u_size - 1
            do l1 = 1, tile_size_param + 10
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - l1*u_size - 1
            do l1 = 1, tile_size_param + 11
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - l1*u_size - 1
            do l1 = 1, tile_size_param + 12
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - l1*u_size - 1
            do l1 = 1, tile_size_param + 13
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - l1*u_size - 1
            do l1 = 1, tile_size_param + 14
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - l1*u_size - 1
            do l1 = 1, tile_size_param + 15
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - l1*u_size - 1
            do l1 = 1, tile_size_param + 16
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  error = error + Abs(u_old - u(i))
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - l1*u_size + subtile_level_param + tile_size_param + u_size*(subtile_level_param + 1)
         end do
         do l1 = 1, tile_size_param + 0
            do l0 = 1, tile_size_param + 0
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - l1*u_size - 1
         do l1 = 1, tile_size_param + 1
            do l0 = 1, tile_size_param + 1
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - l1*u_size - 1
         do l1 = 1, tile_size_param + 2
            do l0 = 1, tile_size_param + 2
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - l1*u_size - 1
         do l1 = 1, tile_size_param + 3
            do l0 = 1, tile_size_param + 3
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - l1*u_size - 1
         do l1 = 1, tile_size_param + 4
            do l0 = 1, tile_size_param + 4
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - l1*u_size - 1
         do l1 = 1, tile_size_param + 5
            do l0 = 1, tile_size_param + 5
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - l1*u_size - 1
         do l1 = 1, tile_size_param + 6
            do l0 = 1, tile_size_param + 6
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - l1*u_size - 1
         do l1 = 1, tile_size_param + 7
            do l0 = 1, tile_size_param + 7
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - l1*u_size - 1
         do l1 = 1, tile_size_param + 8
            do l0 = 1, tile_size_param + 8
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - l1*u_size - 1
         do l1 = 1, tile_size_param + 9
            do l0 = 1, tile_size_param + 9
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - l1*u_size - 1
         do l1 = 1, tile_size_param + 10
            do l0 = 1, tile_size_param + 10
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - l1*u_size - 1
         do l1 = 1, tile_size_param + 11
            do l0 = 1, tile_size_param + 11
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - l1*u_size - 1
         do l1 = 1, tile_size_param + 12
            do l0 = 1, tile_size_param + 12
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - l1*u_size - 1
         do l1 = 1, tile_size_param + 13
            do l0 = 1, tile_size_param + 13
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - l1*u_size - 1
         do l1 = 1, tile_size_param + 14
            do l0 = 1, tile_size_param + 14
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - l1*u_size - 1
         do l1 = 1, tile_size_param + 15
            do l0 = 1, tile_size_param + 15
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - l1*u_size - 1
         do l1 = 1, tile_size_param + 16
            do l0 = 1, tile_size_param + 16
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               error = error + Abs(u_old - u(i))
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - l1*u_size - 1
         if (error < eps) then
            return
         end if
      end do
   end subroutine subtiling_sor_16

end module subgrid_solver
