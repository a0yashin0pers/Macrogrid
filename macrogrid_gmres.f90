program macrogrid_solver
    implicit none
    logical :: one_sor = .true.
    integer, parameter :: max_iter = 1000000
    real*8, parameter :: epsilon = 1.0d-5
    integer, dimension(8) :: u_sizes = [10, 18, 34, 66, 130, 258, 514, 1026]
    real*8, dimension(8) :: factors = [1.52729d0, 1.69504d0, 1.82964d0, 1.90932d0, 1.95305d0, 1.97616d0, 1.98798d0, 1.99394d0]
    integer :: i, subgrid_size, num_subgrids, iterations
    real*8 :: omega
    real*8 :: macrogrid(10, 10, 1026, 1026)

    i = 1
    one_sor = .false.
    subgrid_size = u_sizes(i)
    omega = factors(i)
    num_subgrids = 2

    call initialize_macrogrid(macrogrid, num_subgrids, subgrid_size)
    call solve_macrogrid_with_gmres(one_sor, macrogrid, num_subgrids, subgrid_size, omega, max_iter, epsilon, iterations)
    call compute_error(macrogrid, num_subgrids, subgrid_size)

end program macrogrid_solver

subroutine print_macrogrid(macrogrid, num_subgrids, subgrid_size)
    implicit none
    integer, intent(in) :: num_subgrids, subgrid_size
    real*8, intent(in) :: macrogrid(num_subgrids, num_subgrids, subgrid_size, subgrid_size)
    integer :: i, j, x, y

    do i = 1, num_subgrids
        do y = 1, subgrid_size
            do j = 1, num_subgrids
                do x = 1, subgrid_size
                    write(*, '(10F8.3)', advance='no') macrogrid(j, i, x, y)
                end do
            end do
            print *, ' '
        end do
    end do
end subroutine print_macrogrid

subroutine initialize_macrogrid(macrogrid, num_subgrids, subgrid_size)
    implicit none
    integer, intent(in) :: num_subgrids, subgrid_size
    real*8, intent(inout) :: macrogrid(num_subgrids, num_subgrids, subgrid_size, subgrid_size)
    integer :: u_size, i0, j0, i1, j1, l0, l1

    u_size = num_subgrids*subgrid_size - num_subgrids + 1

    do i0 = 1, num_subgrids
        do j0 = 1, num_subgrids
            do i1 = 1, subgrid_size
                do j1 = 1, subgrid_size

                    l1 = i1 + (i0-1)*subgrid_size - i0 + 1
                    l0 = j1 + (j0-1)*subgrid_size - j0 + 1

                    if (l0.eq.1 .or. l0.eq.u_size .or. l1.eq.1 .or. l1.eq.u_size) then
                        macrogrid(i0,j0,i1,j1) = 1.0d0
                    else
                        macrogrid(i0,j0,i1,j1) = 0.0d0
                    end if
        
                end do
            end do
        end do
    end do
end subroutine initialize_macrogrid

subroutine solve_sor(u, u_size, factor, eps, max_iter, iter_error)
    implicit none
    integer, intent(in) :: u_size, max_iter
    real*8, intent(in) :: factor, eps
    real*8, intent(inout) :: u(u_size*u_size)
    real*8, intent(out) :: iter_error
    integer :: i, l0, l1, iter
    real*8 :: f0, f1, u_old

    f0 = 1.0d0 - factor
    f1 = factor * 0.25d0

    do iter = 1, max_iter
        iter_error = 0.0d0
        i = u_size + 2

        do l1 = 3, u_size
            do l0 = 3, u_size

                u_old = u(i)
                u(i) = f0 * u_old + f1 * (u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))

                iter_error = iter_error + abs(u(i) - u_old)
                i = i + 1

            end do
            i = i + 2
        end do

        if (iter_error < eps) exit
    end do

end subroutine solve_sor

subroutine scalar_product(a, b, size, result)
    implicit none
    integer, intent(in) :: size
    real*8, intent(in) :: a(size)
    real*8, intent(in) :: b(size)
    real*8, intent(out) :: result
    integer :: i

    result = 0.0d0
    do i = 1, size
        result = result + a(i) * b(i)
    end do

end subroutine scalar_product

subroutine get_border(macrogrid, num_subgrids, subgrid_size, border_size, border)
    implicit none
    integer, intent(in) :: num_subgrids, subgrid_size, border_size
    real*8, intent(in) :: macrogrid(num_subgrids, num_subgrids, subgrid_size, subgrid_size)
    real*8, intent(out) :: border(border_size)
    integer :: i, j, k, l

    l = 1

    do i = 1, num_subgrids
        do j = 1, num_subgrids-1
            do k = 2, subgrid_size-1
                border(l) = macrogrid(i,j,k,subgrid_size)
                l = l+1
            end do
        end do
    end do

    do i = 1, num_subgrids-1
        do j = 1, num_subgrids
            do k = 2, subgrid_size-1
                border(l) = macrogrid(i,j,subgrid_size,k)
                l = l+1
            end do
        end do
    end do

end subroutine get_border

subroutine initialize_b(macrogrid, num_subgrids, subgrid_size, border_size, b)
    implicit none
    integer, intent(in) :: num_subgrids, subgrid_size, border_size
    real*8, intent(in) :: macrogrid(num_subgrids, num_subgrids, subgrid_size, subgrid_size)
    real*8, intent(out) :: b(border_size)
    integer :: i, j, k, l

    l = 1

    do i = 1, num_subgrids
        do j = 1, num_subgrids-1
            do k = 2, subgrid_size-1
                b(l) = - 4.0d0*macrogrid(i,j,k,subgrid_size-1) - 4.0d0*macrogrid(i,j+1,k,2) + &
                macrogrid(i,j,k,subgrid_size-2) + macrogrid(i,j+1,k,3)
                l = l+1
            end do
        end do
    end do

    do i = 1, num_subgrids-1
        do j = 1, num_subgrids
            do k = 2, subgrid_size-1
                b(l) = - 4.0d0*macrogrid(i,j,subgrid_size-1,k) - 4.0d0*macrogrid(i+1,j,2,k) + &
                macrogrid(i,j,subgrid_size-2,k) + macrogrid(i+1,j,3,k)
                l = l+1
            end do
        end do
    end do

end subroutine initialize_b

subroutine Cp(macrogrid, p, b, num_subgrids, subgrid_size, border_size, result)
    implicit none
    integer, intent(in) :: num_subgrids, subgrid_size, border_size
    real*8, intent(inout) :: macrogrid(num_subgrids, num_subgrids, subgrid_size, subgrid_size)
    real*8, intent(in) :: b(border_size)
    real*8, intent(in) :: p(border_size)
    real*8, intent(out) :: result(border_size)
    integer :: i, j, k, l

    l = 1

    do i = 1, num_subgrids
        do j = 1, num_subgrids-1
            do k = 2, subgrid_size-1
                result(l) = - 4.0d0*macrogrid(i,j,k,subgrid_size-1) - 4.0d0*macrogrid(i,j+1,k,2) + &
                macrogrid(i,j,k,subgrid_size-2) + macrogrid(i,j+1,k,3) + &
                6.0d0 * p(l) - b(l)
                l = l + 1
            end do
        end do
    end do

    do i = 1, num_subgrids-1
        do j = 1, num_subgrids
            do k = 2, subgrid_size-1
                result(l) = - 4.0d0*macrogrid(i,j,subgrid_size-1,k) - 4.0d0*macrogrid(i+1,j,2,k) + &
                macrogrid(i,j,subgrid_size-2,k) + macrogrid(i+1,j,3,k) + &
                6.0d0 * p(l) - b(l)
                l = l + 1
            end do
        end do
    end do

end subroutine Cp

subroutine update_boundaries_gmres(macrogrid, b, num_subgrids, subgrid_size, border_size, max_iter, epsilon, norm)
    implicit none
    integer, intent(in) :: num_subgrids, subgrid_size, border_size, max_iter
    real*8, intent(in) :: epsilon
    real*8, intent(in) :: b(border_size)
    real*8, intent(inout) :: macrogrid(num_subgrids, num_subgrids, subgrid_size, subgrid_size)
    real*8, intent(out) :: norm
    real*8 :: r(border_size)
    real*8 :: r_1(border_size)
    real*8 :: p(border_size)
    real*8 :: p_1(border_size)
    real*8 :: u(border_size)
    real*8 :: u_old(border_size)
    real*8 :: Ax(border_size)
    real*8 :: Ar(border_size)
    real*8 :: Ar_1(border_size)
    real*8 :: Ap(border_size)
    real*8 :: Ap_1(border_size)
    real*8 :: alpha, beta
    integer :: iter, i, j, k, l
    real*8 :: temp0, temp1

    call get_border(macrogrid, num_subgrids, subgrid_size, border_size, u)
    call Cp(macrogrid, u, b, num_subgrids, subgrid_size, border_size, Ax)

    r = - b - Ax
    p = r(:)

    call Cp(macrogrid, r, b, num_subgrids, subgrid_size, border_size, Ar)
    call scalar_product(Ar,r,border_size,temp0)
    call scalar_product(Ar,Ar,border_size,temp1)

    alpha = temp0/temp1
    u = u + alpha*p

    Ap_1 = Ar(:)
    Ar_1 = Ar(:)
    r_1 = r(:)
    p_1 = p(:)

    do iter = 1, max_iter
        r = r_1 - alpha * Ap_1
        call Cp(macrogrid, r, b, num_subgrids, subgrid_size, border_size, Ar)
        call scalar_product(Ar,r,border_size,temp0)
        call scalar_product(Ar_1,r_1,border_size,temp1)
        beta = temp0/temp1
        p = r + beta*p_1
        Ap = Ar + beta*Ap_1
        u_old = u(:)
        u = u + alpha*p_1

        Ap_1 = Ar(:)
        Ar_1 = Ap(:)
        r_1 = r(:)
        p_1 = p(:)

        call scalar_product(Ar,r,border_size,temp0)
        call scalar_product(Ap,Ap,border_size,temp1)

        alpha = temp0/temp1

        norm = 0.0d0
        do k = 1, border_size
            norm = norm + abs(u_old(k) - u(k))
        end do

        if (norm < epsilon) then
            exit
        end if

    end do

    write (*,*) iter

    l = 1
    norm = 0.0d0

    do i = 1, num_subgrids
        do j = 1, num_subgrids-1
            do k = 2, subgrid_size-1
                norm = norm + abs(u(l) - macrogrid(i,j,k,subgrid_size))
                macrogrid(i,j,k,subgrid_size) = u(l)
                macrogrid(i,j+1,k,1) = u(l)
                l = l + 1
            end do
        end do
    end do

    do i = 1, num_subgrids-1
        do j = 1, num_subgrids
            do k = 2, subgrid_size-1
                norm = norm + abs(u(l) - macrogrid(i,j,subgrid_size,k))
                macrogrid(i,j,subgrid_size,k) = u(l)
                macrogrid(i+1,j,1,k) = u(l)
                l = l + 1
            end do
        end do
    end do

    do i = 1, num_subgrids-1
        do j = 1, num_subgrids-1
            temp0 = 0.25d0 * (macrogrid(i,j,subgrid_size,subgrid_size-1) + &
            macrogrid(i,j,subgrid_size-1,subgrid_size) + &
            macrogrid(i+1,j+1,2,1) + &
            macrogrid(i+1,j+1,1,2))

            macrogrid(i,j,subgrid_size,subgrid_size) = temp0
            macrogrid(i+1,j,1,subgrid_size) = temp0
            macrogrid(i,j+1,subgrid_size,1) = temp0
            macrogrid(i+1,j+1,1,1) = temp0
            l = l + 1
        end do
    end do

end subroutine update_boundaries_gmres

subroutine solve_macrogrid_with_gmres(one_sor, macrogrid, num_subgrids, subgrid_size, omega, max_iter, epsilon, iterations)
    implicit none
    logical, intent(in) :: one_sor
    integer, intent(in) :: num_subgrids, subgrid_size, max_iter
    real*8, intent(in) :: omega, epsilon
    real*8, intent(inout) :: macrogrid(num_subgrids, num_subgrids, subgrid_size, subgrid_size)
    integer, intent(out) :: iterations
    integer :: border_size
    real*8 :: b(2*num_subgrids*(num_subgrids-1)*(subgrid_size-2))
    integer :: i, j, iter, max_iter_sor
    real*8 :: norm, norm_sor, norm_boundaries

    border_size = 2*num_subgrids*(num_subgrids-1)*(subgrid_size-2)

    max_iter_sor = max_iter
    if (one_sor) max_iter_sor = 1

    iterations = 0

    norm = 0.0d0

    call print_macrogrid(macrogrid, num_subgrids, subgrid_size)
    print *, " "
    print *, " "

    do i = 1, num_subgrids
        do j = 1, num_subgrids
            call solve_sor(macrogrid(i,j,:,:), subgrid_size, omega, epsilon, max_iter_sor, norm_sor)
            norm = norm + norm_sor
        end do
    end do

    call initialize_b(macrogrid, num_subgrids, subgrid_size, border_size, b)

    call update_boundaries_gmres(macrogrid, b, num_subgrids, subgrid_size, border_size, max_iter, epsilon, norm_boundaries)

    call print_macrogrid(macrogrid, num_subgrids, subgrid_size)
    print *, " "
    print *, " "

    do iter = 1, max_iter
        norm = 0.0d0

        do i = 1, num_subgrids
            do j = 1, num_subgrids
                call solve_sor(macrogrid(i,j,:,:), subgrid_size, omega, epsilon, max_iter_sor, norm_sor)
                norm = norm + norm_sor
            end do
        end do

        if (.not.one_sor) norm = 0.0d0

        call update_boundaries_gmres(macrogrid, b, num_subgrids, subgrid_size, border_size, max_iter, epsilon, norm_boundaries)
        norm = norm + norm_boundaries

        if (norm < epsilon) then
            iterations = iter
            return
        end if
    end do

    call print_macrogrid(macrogrid, num_subgrids, subgrid_size)
    print *, " "
    print *, " "

    iterations = max_iter
end subroutine solve_macrogrid_with_gmres

subroutine compute_error(macrogrid, num_subgrids, subgrid_size)
    implicit none
    integer, intent(in) :: num_subgrids, subgrid_size
    real*8, intent(in) :: macrogrid(num_subgrids, num_subgrids, subgrid_size, subgrid_size)
    real*8 :: max_error, diff
    integer :: i, j, x, y

    max_error = 0.0d0

    do i = 1, num_subgrids
        do j = 1, num_subgrids
            do x = 1, subgrid_size
                do y = 1, subgrid_size
                    diff = abs(macrogrid(i,j,x,y) - 1.0d0)
                    if (diff > max_error) max_error = diff
                end do
            end do
        end do
    end do

    print *, 'Maximum error compared to the exact solution: ', max_error
end subroutine compute_error