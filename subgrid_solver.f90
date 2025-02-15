subroutine original (u, u_size, factor, eps, max_iter, error, dx, dy)
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
end subroutine original

subroutine tiling (u, u_size, factor, eps, max_iter, error, dx, dy)
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
end subroutine tiling

subroutine subtiling (u, u_size, factor, eps, max_iter, error, dx, dy)
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
end subroutine subtiling