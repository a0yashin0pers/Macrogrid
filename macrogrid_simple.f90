program macrogrid_solver
    implicit none

    ! Переключатель — нужно ли накапливать ошибку в подобластях при подсчёте нормы на интерфейсе
    logical :: accumulate_subgrid_error = .true.

    ! Списки оптимальных параметров релаксации и размеров подсеток
    integer, parameter :: NTEST = 8
    integer, dimension(NTEST), parameter :: u_sizes = [10, 18, 34, 66, 130, 258, 514, 1026]

    real*8, dimension(NTEST), parameter :: factors =      &
         [1.52729d0, 1.69504d0, 1.82964d0, 1.90932d0, &
          1.95305d0, 1.97616d0, 1.98798d0, 1.99394d0]

    real*8, dimension(NTEST), parameter :: factors_hard = &
         [1.52129d0, 1.69524d0, 1.82919d0, 1.90899d0, &
          1.95307d0, 1.97618d0, 1.98793d0, 1.99403d0]

    ! Индексы тестов (какие взять размеры из массива u_sizes)
    integer, dimension(NTEST), parameter :: test_size_indices = [1,2,3,4,5,6,7,8]
    ! Число подсеток (по x или y) в каждом тесте — для упрощения делаем одинаковое по x и y
    integer, dimension(NTEST), parameter :: test_subgrids = [4,4,4,4,4,4,4,4]

    ! Параметры итераций и точности
    integer, parameter :: max_iter_subgrid    = 100000   ! Максимум итераций SOR в подобласти
    integer, parameter :: max_iter_interface  = 100000   ! Максимум итераций обновления интерфейсов
    real*8,  parameter :: eps_subgrid         = 1.0d-10  ! Точность для SOR внутри подобластей
    real*8,  parameter :: eps_interface       = 1.0d-5   ! Точность для итераций на интерфейсе

    ! Массив для макросетки. Зададим максимальные границы 10x10 подобластей и до 1026x1026 внутри
    real*8 :: macrogrid(10, 10, 1026, 1026)

    ! Локальные рабочие переменные в главной программе
    integer :: i, idx, subgrid_size, Nx, Ny
    real*8  :: omega
    integer :: one_sor_ind

    print *, "ones task"

    ! Два варианта: accumulate_subgrid_error = .true. и .false.
    do one_sor_ind = 1, 2

        do i = 1, size(test_size_indices)
            ! Текущий размер подзадачи
            subgrid_size = u_sizes(test_size_indices(i))
            ! Параметр релаксации
            omega        = factors(test_size_indices(i))
            ! Число подсеток (оба направления одинаково)
            Nx = test_subgrids(i)
            Ny = test_subgrids(i)

            print *, "accumulate_subgrid_error = ", accumulate_subgrid_error
            print *, "Running test with subgrid size=", subgrid_size, " Nx=", Nx, " Ny=", Ny

            ! Инициализируем "единичную" задачу
            call initialize_macrogrid_ones(macrogrid, Nx, Ny, subgrid_size)

            ! Запускаем измерение времени + решение
            call measure_time(accumulate_subgrid_error,                   &
                             macrogrid, Nx, Ny, subgrid_size,             &
                             omega, max_iter_subgrid, max_iter_interface, &
                             eps_subgrid, eps_interface)

            ! Считаем ошибку (где граничные значения = 1)
            call compute_error(macrogrid, Nx, Ny, subgrid_size)
            print *, " "
        end do

        ! Меняем флаг
        accumulate_subgrid_error = .not. accumulate_subgrid_error
    end do


    print *, "hard task"

    ! Снова два варианта accumulate_subgrid_error
    do one_sor_ind = 1, 2

        do i = 1, size(test_size_indices)
            subgrid_size = u_sizes(test_size_indices(i))
            omega        = factors_hard(test_size_indices(i))
            Nx = test_subgrids(i)
            Ny = test_subgrids(i)

            print *, "accumulate_subgrid_error = ", accumulate_subgrid_error
            print *, "Running test with subgrid size=", subgrid_size, " Nx=", Nx, " Ny=", Ny

            call initialize_hard_macrogrid(macrogrid, Nx, Ny, subgrid_size)

            call measure_time(accumulate_subgrid_error,                   &
                             macrogrid, Nx, Ny, subgrid_size,             &
                             omega, max_iter_subgrid, max_iter_interface, &
                             eps_subgrid, eps_interface)

            call compute_hard_error(macrogrid, Nx, Ny, subgrid_size)
            print *, " "
        end do

        accumulate_subgrid_error = .not. accumulate_subgrid_error
    end do

end program macrogrid_solver


!-----------------------------------------------------------------------
! ИНИЦИАЛИЗАЦИЯ: "ЕДИНИЧКИ" НА ГРАНИЦАХ
!-----------------------------------------------------------------------
subroutine initialize_macrogrid_ones(macrogrid, num_subgrids_x, num_subgrids_y, subgrid_size)
    implicit none
    integer, intent(in) :: num_subgrids_x, num_subgrids_y, subgrid_size
    real*8,  intent(inout) :: macrogrid(num_subgrids_x,num_subgrids_y,subgrid_size,subgrid_size)

    integer :: global_size_x, global_size_y
    integer :: iX, iY, i1, j1, lX, lY

    ! Глобальное число узлов по x и y
    global_size_x = num_subgrids_x * subgrid_size - (num_subgrids_x - 1)
    global_size_y = num_subgrids_y * subgrid_size - (num_subgrids_y - 1)

    do iX = 1, num_subgrids_x
        do iY = 1, num_subgrids_y
            do i1 = 1, subgrid_size
                do j1 = 1, subgrid_size
                    ! Глобальные индексы
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
end subroutine initialize_macrogrid_ones


!-----------------------------------------------------------------------
! ИНИЦИАЛИЗАЦИЯ: "ЛОГАРИФМИЧЕСКИЕ" ГРАНИЧНЫЕ УСЛОВИЯ
!-----------------------------------------------------------------------
subroutine initialize_hard_macrogrid(macrogrid, num_subgrids_x, num_subgrids_y, subgrid_size)
    implicit none
    real*8, parameter :: R1 = 0.1d0, R2 = 1.0d0
    real*8, parameter :: x_min = 0.3d0, y_min = 0.0d0
    real*8, parameter :: len = 0.4d0

    integer, intent(in) :: num_subgrids_x, num_subgrids_y, subgrid_size
    real*8,  intent(inout) :: macrogrid(num_subgrids_x,num_subgrids_y,subgrid_size,subgrid_size)

    integer :: global_size_x, global_size_y
    integer :: iX, iY, i1, j1, lX, lY
    real*8  :: boundary_val

    global_size_x = num_subgrids_x * subgrid_size - (num_subgrids_x - 1)
    global_size_y = num_subgrids_y * subgrid_size - (num_subgrids_y - 1)

    do iX = 1, num_subgrids_x
        do iY = 1, num_subgrids_y
            do i1 = 1, subgrid_size
                do j1 = 1, subgrid_size
                    lX = i1 + (iX-1)*subgrid_size - (iX-1)
                    lY = j1 + (iY-1)*subgrid_size - (iY-1)

                    ! Граница -> логарифмическое значение
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
end subroutine initialize_hard_macrogrid


!-----------------------------------------------------------------------
! ВЫЧИСЛЕНИЕ ОШИБКИ ДЛЯ "ЕДИНИЧКИ" НА ГРАНИЦЕ
!-----------------------------------------------------------------------
subroutine compute_error(macrogrid, num_subgrids_x, num_subgrids_y, subgrid_size)
    implicit none
    integer, intent(in) :: num_subgrids_x, num_subgrids_y, subgrid_size
    real*8,  intent(in) :: macrogrid(num_subgrids_x,num_subgrids_y,subgrid_size,subgrid_size)

    integer :: iX, iY, i1, j1
    real*8  :: max_error, diff

    max_error = 0.0d0

    do iX = 1, num_subgrids_x
        do iY = 1, num_subgrids_y
            do i1 = 1, subgrid_size
                do j1 = 1, subgrid_size
                    diff = abs(macrogrid(iX,iY,i1,j1) - 1.0d0)
                    if (diff > max_error) max_error = diff
                end do
            end do
        end do
    end do

    print *, 'Maximum error compared to the exact solution (ones): ', max_error
end subroutine compute_error


!-----------------------------------------------------------------------
! ВЫЧИСЛЕНИЕ ОШИБКИ ДЛЯ "ЛОГАРИФМИЧЕСКИХ" ГРАНИЧНЫХ УСЛОВИЙ
!-----------------------------------------------------------------------
subroutine compute_hard_error(macrogrid, num_subgrids_x, num_subgrids_y, subgrid_size)
    implicit none
    real*8, parameter :: R1 = 0.1d0, R2 = 1.0d0
    real*8, parameter :: x_min = 0.3d0, y_min = 0.0d0
    real*8, parameter :: len = 0.4d0

    integer, intent(in) :: num_subgrids_x, num_subgrids_y, subgrid_size
    real*8,  intent(in) :: macrogrid(num_subgrids_x,num_subgrids_y,subgrid_size,subgrid_size)

    integer :: global_size_x, global_size_y
    integer :: iX, iY, i1, j1, lX, lY
    real*8  :: max_error, val_exact, diff

    global_size_x = num_subgrids_x * subgrid_size - (num_subgrids_x - 1)
    global_size_y = num_subgrids_y * subgrid_size - (num_subgrids_y - 1)

    max_error = 0.0d0

    do iX = 1, num_subgrids_x
        do iY = 1, num_subgrids_y
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

    print *, 'Maximum error compared to the exact solution (log): ', max_error
end subroutine compute_hard_error


!-----------------------------------------------------------------------
! ПРОЦЕДУРА МЕТОДА SOR В ОДНОЙ ПОДОБЛАСТИ (n x n)
!-----------------------------------------------------------------------
subroutine solve_sor_subdomain(u, n, omega, eps_sor, max_iter_sor, out_error)
    implicit none
    integer, intent(in) :: n, max_iter_sor
    real*8,  intent(in) :: omega, eps_sor
    real*8,  intent(inout) :: u(n*n)
    real*8,  intent(out)   :: out_error

    integer :: iter, ix, iy, idx
    real*8  :: old_val, sum_diff

    do iter = 1, max_iter_sor
        sum_diff = 0.0d0

        ! Обход внутренних узлов (2..n-1)
        do iy = 2, n-1
            do ix = 2, n-1
                idx     = ix + (iy-1)*n
                old_val = u(idx)

                ! Классическая формула SOR (5-точечный лапласиан)
                ! u_new = u_old + omega*(0.25*(левый + правый + верх + низ) - u_old)
                u(idx)  = old_val + omega*0.25d0*( u(idx-1) + u(idx+1) + &
                                                  u(idx-n) + u(idx+n) - 4.0d0*old_val )

                sum_diff = sum_diff + dabs(u(idx) - old_val)
            end do
        end do

        if (sum_diff < eps_sor) then
            out_error = sum_diff
            return
        end if
    end do

    out_error = sum_diff
end subroutine solve_sor_subdomain


!-----------------------------------------------------------------------
! КОПИРОВАНИЕ ДВУМЕРНОГО subgrid(n,n) В 1D u(n*n)
!-----------------------------------------------------------------------
subroutine copy_in_subgrid(block_2d, local_u, n)
    implicit none
    integer, intent(in) :: n
    real*8,  intent(in) :: block_2d(n,n)
    real*8,  intent(out):: local_u(n*n)

    integer :: ix, iy, idx

    do iy = 1, n
        do ix = 1, n
            idx          = ix + (iy-1)*n
            local_u(idx) = block_2d(ix, iy)
        end do
    end do
end subroutine copy_in_subgrid


!-----------------------------------------------------------------------
! КОПИРОВАНИЕ ОБРАТНО: ИЗ 1D u(n*n) В block_2d(n,n)
!-----------------------------------------------------------------------
subroutine copy_out_subgrid(block_2d, local_u, n)
    implicit none
    integer, intent(in) :: n
    real*8,  intent(in) :: local_u(n*n)
    real*8,  intent(inout):: block_2d(n,n)

    integer :: ix, iy, idx

    do iy = 1, n
        do ix = 1, n
            idx                = ix + (iy-1)*n
            block_2d(ix, iy)   = local_u(idx)
        end do
    end do
end subroutine copy_out_subgrid


!-----------------------------------------------------------------------
! ОБНОВЛЕНИЕ УЗЛОВ НА ГРАНИЦАХ (ПРОСТАЯ ИТЕРАЦИЯ)
!-----------------------------------------------------------------------
subroutine update_boundaries(macrogrid, num_subgrids_x, num_subgrids_y, subgrid_size, out_norm)
    implicit none
    integer, intent(in) :: num_subgrids_x, num_subgrids_y, subgrid_size
    real*8,  intent(inout) :: macrogrid(num_subgrids_x,num_subgrids_y,subgrid_size,subgrid_size)
    real*8,  intent(out)   :: out_norm

    integer :: iX, iY, k
    real*8  :: old_value, avg

    out_norm = 0.0d0

    ! Горизонтальные "стыки" между блоками (iY / iY+1)
    do iX = 1, num_subgrids_x
        do iY = 1, num_subgrids_y - 1
            do k = 2, subgrid_size - 1
                old_value = macrogrid(iX, iY, k, subgrid_size)
                avg = ( 4.0d0*macrogrid(iX, iY,   k, subgrid_size-1) + &
                        4.0d0*macrogrid(iX, iY+1, k, 2)               - &
                        macrogrid(iX, iY,   k, subgrid_size-2)        - &
                        macrogrid(iX, iY+1, k, 3 ) ) / 6.0d0

                macrogrid(iX, iY,   k, subgrid_size) = avg
                macrogrid(iX, iY+1, k, 1)            = avg
                out_norm = out_norm + dabs(avg - old_value)
            end do
        end do
    end do

    ! Вертикальные "стыки" между блоками (iX / iX+1)
    do iX = 1, num_subgrids_x - 1
        do iY = 1, num_subgrids_y
            do k = 2, subgrid_size - 1
                old_value = macrogrid(iX, iY, subgrid_size, k)
                avg = ( 4.0d0*macrogrid(iX,   iY, subgrid_size-1, k) + &
                        4.0d0*macrogrid(iX+1, iY, 2,             k)   - &
                        macrogrid(iX,   iY, subgrid_size-2, k)         - &
                        macrogrid(iX+1, iY, 3,             k) ) / 6.0d0

                macrogrid(iX,   iY, subgrid_size, k) = avg
                macrogrid(iX+1, iY, 1,            k) = avg
                out_norm = out_norm + dabs(avg - old_value)
            end do
        end do
    end do

    ! "Углы" между четырьмя блоками
    do iX = 1, num_subgrids_x - 1
        do iY = 1, num_subgrids_y - 1
            old_value = macrogrid(iX, iY, subgrid_size, subgrid_size)
            avg = ( 4.0d0 * macrogrid(iX,   iY,   subgrid_size,   subgrid_size-1) + &
                    4.0d0 * macrogrid(iX,   iY,   subgrid_size-1, subgrid_size)   + &
                    4.0d0 * macrogrid(iX+1, iY+1, 2,              1)              + &
                    4.0d0 * macrogrid(iX+1, iY+1, 1,              2)              - &
                    macrogrid(iX,   iY,   subgrid_size,   subgrid_size-2)         - &
                    macrogrid(iX,   iY,   subgrid_size-2, subgrid_size)           - &
                    macrogrid(iX+1, iY+1, 3,              1)                      - &
                    macrogrid(iX+1, iY+1, 1,              3) ) / 12.0d0

            macrogrid(iX,   iY,   subgrid_size,   subgrid_size)   = avg
            macrogrid(iX+1, iY,   1,              subgrid_size)   = avg
            macrogrid(iX,   iY+1, subgrid_size,   1)              = avg
            macrogrid(iX+1, iY+1, 1,              1)              = avg
            out_norm = out_norm + dabs(avg - old_value)
        end do
    end do

end subroutine update_boundaries


!-----------------------------------------------------------------------
! РЕШЕНИЕ НА ВСЕЙ МАКРОСЕТКЕ (ДЕКОМПОЗИЦИЯ) — ОСНОВНОЙ ЦИКЛ
!-----------------------------------------------------------------------
subroutine solve_macrogrid(accumulate_subgrid_error,                          &
                           macrogrid, num_subgrids_x, num_subgrids_y,         &
                           subgrid_size, omega,                               &
                           max_iter_sor, max_iter_interface,                  &
                           eps_sor, eps_interface,                            &
                           out_iterations)
    implicit none
    logical, intent(in) :: accumulate_subgrid_error
    integer, intent(in) :: num_subgrids_x, num_subgrids_y
    integer, intent(in) :: subgrid_size, max_iter_sor, max_iter_interface
    real*8,  intent(in) :: omega, eps_sor, eps_interface
    real*8,  intent(inout) :: macrogrid(num_subgrids_x,num_subgrids_y,subgrid_size,subgrid_size)
    integer, intent(out) :: out_iterations

    integer :: iter, iX, iY
    real*8  :: norm_total, norm_sub_sor, norm_bound
    real*8, allocatable :: local_u(:)

    out_iterations = 0

    do iter = 1, max_iter_interface

        norm_total = 0.0d0

        ! (1) Решаем SOR в каждой подобласти
        do iX = 1, num_subgrids_x
            do iY = 1, num_subgrids_y

                allocate(local_u(subgrid_size*subgrid_size))
                call copy_in_subgrid(macrogrid(iX,iY,:,:), local_u, subgrid_size)

                call solve_sor_subdomain(local_u, subgrid_size, omega, eps_sor, max_iter_sor, norm_sub_sor)

                if (accumulate_subgrid_error) then
                    norm_total = norm_total + norm_sub_sor
                end if

                call copy_out_subgrid(macrogrid(iX,iY,:,:), local_u, subgrid_size)
                deallocate(local_u)
            end do
        end do

        ! (2) Обновляем границы
        call update_boundaries(macrogrid, num_subgrids_x, num_subgrids_y, subgrid_size, norm_bound)
        norm_total = norm_total + norm_bound

        ! Проверяем сходимость на уровне макросетки
        if (norm_total < eps_interface) then
            out_iterations = iter
            return
        end if
    end do

    out_iterations = max_iter_interface
end subroutine solve_macrogrid


!-----------------------------------------------------------------------
! ПРОСТАЯ РАСЧЁТ ВРЕМЕНИ И ВЫЗОВ solve_macrogrid
!-----------------------------------------------------------------------
subroutine measure_time(accumulate_subgrid_error,                            &
                        macrogrid, num_subgrids_x, num_subgrids_y,           &
                        subgrid_size, omega,                                 &
                        max_iter_sor, max_iter_interface,                    &
                        eps_sor, eps_interface)
    implicit none
    logical, intent(in) :: accumulate_subgrid_error
    integer, intent(in) :: num_subgrids_x, num_subgrids_y
    integer, intent(in) :: subgrid_size, max_iter_sor, max_iter_interface
    real*8,  intent(in) :: omega, eps_sor, eps_interface
    real*8,  intent(inout) :: macrogrid(num_subgrids_x,num_subgrids_y,subgrid_size,subgrid_size)

    real*8 :: start_time, end_time, elapsed_time
    integer :: iterations

    call cpu_time(start_time)

    call solve_macrogrid(accumulate_subgrid_error,                &
                         macrogrid, num_subgrids_x, num_subgrids_y,&
                         subgrid_size, omega,                      &
                         max_iter_sor, max_iter_interface,         &
                         eps_sor, eps_interface,                   &
                         iterations)

    call cpu_time(end_time)

    elapsed_time = end_time - start_time
    print *, 'Time to solve macrogrid = ', elapsed_time, ' seconds'
    print *, 'Number of interface iterations = ', iterations

end subroutine measure_time
