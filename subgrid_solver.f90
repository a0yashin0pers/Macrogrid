subroutine original_sor (u, u_size, factor, eps, max_iter, error, dx, dy)
   implicit none

   integer, intent(in) :: u_size, max_iter
   real*8,  intent(in) :: factor, eps, dx, dy
   real*8,  intent(inout) :: u(u_size*u_size)
   real*8,  intent(out)   :: error

   integer :: iter, i, l0, l1, l2, l3
   real*8 :: invdx2, invdy2, f0, f1, u_old

   invdx2 = 1.0d0/(dx*dx)
   invdy2 = 1.0d0/(dy*dy)

   f0 = 1.0d0 - factor
   f1 = factor / (2.0d0*invdx2 + 2.0d0*invdy2)

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

subroutine tiling_sor (u, u_size, factor, eps, max_iter, error, dx, dy)
   implicit none

   integer, parameter :: tile_size = 8

   integer, intent(in) :: u_size, max_iter
   real*8,  intent(in) :: factor, eps, dx, dy
   real*8,  intent(inout) :: u(u_size*u_size)
   real*8,  intent(out)   :: error

   integer :: tile_count, iter, i, l0, l1, l2, l3
   real*8 :: invdx2, invdy2, f0, f1, u_old

   tile_count = (u_size - 2)/tile_size

   invdx2 = 1.0d0/(dx*dx)
   invdy2 = 1.0d0/(dy*dy)

   f0 = 1.0d0 - factor
   f1 = factor / (2.0d0*invdx2 + 2.0d0*invdy2)

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

subroutine subtiling_sor (u, u_size, factor, eps, max_iter, error, dx, dy)
   implicit none

   integer, parameter :: tile_size = 8
   integer, parameter :: tile_level = tile_size

   integer, intent(in) :: u_size, max_iter
   real*8,  intent(in) :: factor, eps, dx, dy
   real*8,  intent(inout) :: u(u_size*u_size)
   real*8,  intent(out)   :: error

   integer :: tile_count, iter, i, l0, l1, l2, l3, l4
   real*8 :: invdx2, invdy2, f0, f1, u_old

   tile_count = (u_size - 2)/tile_size

   invdx2 = 1.0d0/(dx*dx)
   invdy2 = 1.0d0/(dy*dy)

   f0 = 1.0d0 - factor
   f1 = factor / (2.0d0*invdx2 + 2.0d0*invdy2)

   do iter = 1, max_iter

      error = 0.0d0
      i = u_size + 2

      do l2 = 0, tile_level
         do l1 = 1, tile_size - l2
            do l0 = 1, tile_size - l2

               u_old = u(i)
               u(i) = f0*u_old + &
                  f1*((u(i-1) + u(i+1))*invdx2 + &
                  (u(i-u_size) + u(i+u_size))*invdy2)

               if (l2.eq.tile_level) then
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
         do l2 = 0, tile_level
            do l1 = 1, tile_size - l2
               do l0 = 1, tile_size

                  u_old = u(i)
                  u(i) = f0*u_old + &
                     f1*((u(i-1) + u(i+1))*invdx2 + &
                     (u(i-u_size) + u(i+u_size))*invdy2)

                  if (l2.eq.tile_level) then
                     error = error + abs(u(i) - u_old)
                  end if
                  i = i + 1

               end do
               i = i + u_size - l0 + 1
            end do
            i = i - u_size*(l1 - 1) - 1
         end do
         i = i + tile_size + tile_level + 1
      end do

      do l2 = 0, tile_level
         do l1 = 1, tile_size - l2
            do l0 = 1, tile_size + l2

               u_old = u(i)
               u(i) = f0*u_old + &
                  f1*((u(i-1) + u(i+1))*invdx2 + &
                  (u(i-u_size) + u(i+u_size))*invdy2)

               if (l2.eq.tile_level) then
                  error = error + abs(u(i) - u_old)
               end if
               i = i + 1

            end do
            i = i + u_size - l0 + 1
         end do
         i = i - u_size*(l1 - 1) - 1
      end do
      i = i + tile_size + tile_level + u_size*(tile_size - 1) + 3

      do l4 = 2, tile_count - 1
         do l2 = 0, tile_level
            do l1 = 1, tile_size
               do l0 = 1, tile_size - l2

                  u_old = u(i)
                  u(i) = f0*u_old + &
                     f1*((u(i-1) + u(i+1))*invdx2 + &
                     (u(i-u_size) + u(i+u_size))*invdy2)

                  if (l2.eq.tile_level) then
                     error = error + abs(u(i) - u_old)
                  end if
                  i = i + 1

               end do
               i = i + u_size - l0 + 1
            end do
            i = i - u_size*l1
         end do
         i = i + u_size*(tile_level + 1) + tile_size

         do l3 = 2, tile_count - 1
            do l2 = 0, tile_level
               do l1 = 1, tile_size
                  do l0 = 1, tile_size

                     u_old = u(i)
                     u(i) = f0*u_old + &
                        f1*((u(i-1) + u(i+1))*invdx2 + &
                        (u(i-u_size) + u(i+u_size))*invdy2)

                     if (l2.eq.tile_level) then
                        error = error + abs(u(i) - u_old)
                     end if
                     i = i + 1

                  end do
                  i = i + u_size - l0 + 1
               end do
               i = i - u_size*l1 - 1
            end do
            i = i + u_size*(tile_level + 1) + tile_size + tile_level + 1
         end do

         do l2 = 0, tile_level
            do l1 = 1, tile_size
               do l0 = 1, tile_size + l2

                  u_old = u(i)
                  u(i) = f0*u_old + &
                     f1*((u(i-1) + u(i+1))*invdx2 + &
                     (u(i-u_size) + u(i+u_size))*invdy2)

                  if (l2.eq.tile_level) then
                     error = error + abs(u(i) - u_old)
                  end if
                  i = i + 1

               end do
               i = i + u_size - l0 + 1
            end do
            i = i - u_size*l1 - 1
         end do
         i = i + u_size*(tile_level + 1) + tile_size + tile_level + u_size*(tile_size - 1) + 3
      end do

      do l2 = 0, tile_level
         do l1 = 1, tile_size + l2
            do l0 = 1, tile_size - l2

               u_old = u(i)
               u(i) = f0*u_old + &
                  f1*((u(i-1) + u(i+1))*invdx2 + &
                  (u(i-u_size) + u(i+u_size))*invdy2)

               if (l2.eq.tile_level) then
                  error = error + abs(u(i) - u_old)
               end if
               i = i + 1

            end do
            i = i + u_size - l0 + 1
         end do
         i = i - u_size*l1
      end do
      i = i + u_size*(tile_level + 1) + tile_size

      do l3 = 2, tile_count - 1
         do l2 = 0, tile_level
            do l1 = 1, tile_size + l2
               do l0 = 1, tile_size

                  u_old = u(i)
                  u(i) = f0*u_old + &
                     f1*((u(i-1) + u(i+1))*invdx2 + &
                     (u(i-u_size) + u(i+u_size))*invdy2)

                  if (l2.eq.tile_level) then
                     error = error + abs(u(i) - u_old)
                  end if
                  i = i + 1

               end do
               i = i + u_size - l0 + 1
            end do
            i = i - u_size*l1 - 1
         end do
         i = i + u_size*(tile_level + 1) + tile_size + tile_level + 1
      end do

      do l2 = 0, tile_level
         do l1 = 1, tile_size + l2
            do l0 = 1, tile_size + l2

               u_old = u(i)
               u(i) = f0*u_old + &
                  f1*((u(i-1) + u(i+1))*invdx2 + &
                  (u(i-u_size) + u(i+u_size))*invdy2)

               if (l2.eq.tile_level) then
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

subroutine subtiling_sor_test_version(u, u_size, factor, eps, max_iter, error, dx, dy)
   implicit none

   integer, parameter :: tile_size = 4
   integer, parameter :: tile_level = tile_size

   integer, intent(in) :: u_size, max_iter
   real*8,  intent(in) :: factor, eps, dx, dy
   real*8,  intent(inout) :: u(u_size*u_size)
   real*8,  intent(out)   :: error

   integer, dimension(:,:), allocatable :: subtile_indices
   integer :: tile_count, subtile_count, subtile_idx, iter, i, l0, l1, l2, l3, l4, s0, s1, s2, s3
   real*8 :: invdx2, invdy2, f0, f1, u_old

   tile_count = (u_size - 2)/tile_size
   subtile_count = tile_count*tile_count*(tile_level+1)
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

         do l2 = 0, tile_level
            subtile_indices(subtile_idx, 1) = (l4-1)*tile_size*u_size + (l3-1)*tile_size + u_size + 2 - l2*s3*u_size - l2*s2
            subtile_indices(subtile_idx, 2) = tile_size + l2*s1
            subtile_indices(subtile_idx, 3) = tile_size + l2*s0
            if (l2.eq.tile_level) then
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

   f0 = 1.0d0 - factor
   f1 = factor / (2.0d0*invdx2 + 2.0d0*invdy2)

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