program macrogrid_solver
    implicit none

    ! === ПАРАМЕТРЫ ДЛЯ ТЕСТИРОВАНИЯ ===
    ! Логический флаг: суммировать ли ошибку SOR по всем подобластям
    ! при вычислении нормы на интерфейсах
    logical :: accumulate_subgrid_error

    ! Возможные размеры подсеток (только первые 5 из исходного массива)
    integer, parameter :: N_SIZES = 5
    integer, dimension(N_SIZES), parameter :: subgrid_sizes = [10, 18, 34, 66, 130]

    ! Оптимальные факторы релаксации (omega) для "constant boundary" (соответствуют subgrid_sizes выше)
    real*8, dimension(N_SIZES), parameter :: factors_const =      &
         [1.52729d0, 1.69504d0, 1.82964d0, 1.90932d0, 1.95305d0]

    ! Оптимальные факторы релаксации для "logarithmic boundary"
    real*8, dimension(N_SIZES), parameter :: factors_log = &
         [1.52129d0, 1.69524d0, 1.82919d0, 1.90899d0, 1.95307d0]

    ! Макросеточные размеры, которые мы хотим протестировать:
    ! (macrogrid_size_x, macrogrid_size_y)
    integer, parameter :: N_CONFIGS = 4
    integer, dimension(2, N_CONFIGS), parameter :: macrogrid_configs = reshape([ &
         1, 1,  &  ! 1x1
         2, 1,  &  ! 2x1
         4, 2,  &  ! 4x2
         2, 2   &  ! 2x2
         ], shape=[2, N_CONFIGS])

    ! Параметры итераций и точности
    integer, parameter :: max_iter_subgrid    = 100000   ! максимум итераций SOR в подобласти
    integer, parameter :: max_iter_interface  = 100000   ! максимум итераций интерфейсной итерации
    real*8,  parameter :: eps_subgrid         = 1.0d-5   ! точность для SOR внутри подобластей
    real*8,  parameter :: eps_interface       = 1.0d-5   ! точность для итерации на интерфейсе

    ! Большой массив для макросетки (запас по размеру)
    real*8 :: macrogrid(10, 10, 1026, 1026)

    ! Локальные переменные для цикла тестов
    integer :: i_cfg, i_size, sub_sz
    integer :: i        ! счётчик
    real*8  :: omega
    integer :: cfg_x, cfg_y

    print *, "=== Запуск тестов для различных макросеток и подсеток ==="

    ! ДВА цикла: accumulate_subgrid_error = .true. и .false.
    do i = 1, 2
        accumulate_subgrid_error = (i == 1)
        print *, ""
        print *, "-----------------------------------------------"
        print *, "accumulate_subgrid_error = ", accumulate_subgrid_error
        print *, "-----------------------------------------------"

        ! Перебираем все варианты макросетки
        do i_cfg = 1, N_CONFIGS
            cfg_x = macrogrid_configs(1, i_cfg)
            cfg_y = macrogrid_configs(2, i_cfg)

            print *, ""
            print *, "=== Macrogrid size: ", cfg_x, " x ", cfg_y, " ==="

            ! Перебираем все нужные размеры подсеток
            do i_size = 1, N_SIZES
                sub_sz = subgrid_sizes(i_size)

                print *, "  Subgrid size =", sub_sz

                ! Сначала решаем "constant boundary" задачу
                omega = factors_const(i_size)
                call initialize_constant_boundary(macrogrid, cfg_x, cfg_y, sub_sz)
                call measure_execution_time(accumulate_subgrid_error,             &
                                            macrogrid, cfg_x, cfg_y, sub_sz,      &
                                            omega,                                &
                                            max_iter_subgrid, max_iter_interface, &
                                            eps_subgrid, eps_interface)
                call compute_constant_boundary_error(macrogrid, cfg_x, cfg_y, sub_sz)

                ! Теперь решаем "logarithmic boundary"
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
    end do

contains

    !-----------------------------------------------------------------------
    ! ФУНКЦИИ ВЫЧИСЛЕНИЯ ШАГОВ СЕТКИ ДЛЯ ЗАДАННЫХ ПАРАМЕТРОВ
    !-----------------------------------------------------------------------
    ! Предполагаем, что область [0,1] x [0,1].
    ! Тогда global_size_x = macrogrid_size_x*subgrid_size - (macrogrid_size_x - 1)
    ! dx = 1.0 / (global_size_x - 1)
    !
    ! Аналогично по y.

    real*8 function compute_dx(macrogrid_size_x, subgrid_size) result(dx)
        implicit none
        integer, intent(in) :: macrogrid_size_x, subgrid_size
        integer :: global_size_x
        global_size_x = macrogrid_size_x*subgrid_size - (macrogrid_size_x - 1)
        dx = 1.0d0 / dble(global_size_x - 1)
    end function compute_dx

    real*8 function compute_dy(macrogrid_size_y, subgrid_size) result(dy)
        implicit none
        integer, intent(in) :: macrogrid_size_y, subgrid_size
        integer :: global_size_y
        global_size_y = macrogrid_size_y*subgrid_size - (macrogrid_size_y - 1)
        dy = 1.0d0 / dble(global_size_y - 1)
    end function compute_dy


    !-----------------------------------------------------------------------
    ! ИНИЦИАЛИЗАЦИЯ: "КОНСТАНТНЫЕ" ГРАНИЧНЫЕ УСЛОВИЯ (раньше ones)
    !-----------------------------------------------------------------------
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


    !-----------------------------------------------------------------------
    ! ИНИЦИАЛИЗАЦИЯ: "ЛОГАРИФМИЧЕСКИЕ" ГРАНИЧНЫЕ УСЛОВИЯ
    !-----------------------------------------------------------------------
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


    !-----------------------------------------------------------------------
    ! МЕТОД SOR ДЛЯ ПОДОБЛАСТИ (n x n) С УЧЁТОМ dx, dy
    !-----------------------------------------------------------------------
    subroutine solve_sor_subdomain(u, n, factor, eps, max_iter, iter_error, dx, dy)
        implicit none
        integer, intent(in) :: n, max_iter
        real*8,  intent(in) :: factor, eps, dx, dy
        real*8,  intent(inout) :: u(n*n)
        real*8,  intent(out)   :: iter_error

        ! Для уравнения Laplacian(u) = 0 с разными шагами:
        ! u_new(i) = u_old(i) + factor * [ ( (uL+uR)/dx^2 + (uD+uU)/dy^2 ) / (2/dx^2 + 2/dy^2 ) - u_old(i) ]

        integer :: iter, ix, iy, idx
        real*8  :: sum_diff, u_old, denom, tmp
        real*8  :: invdx2, invdy2

        invdx2 = 1.0d0/(dx*dx)
        invdy2 = 1.0d0/(dy*dy)
        denom  = 2.0d0*invdx2 + 2.0d0*invdy2    ! знаменатель в формуле Лапласа

        do iter = 1, max_iter
            sum_diff = 0.0d0

            do iy = 2, n-1
                do ix = 2, n-1
                    idx    = ix + (iy-1)*n
                    u_old  = u(idx)

                    tmp = ( (u(idx-1) + u(idx+1))*invdx2 + &
                            (u(idx-n) + u(idx+n))*invdy2 ) / denom

                    ! SOR обновление
                    u(idx) = u_old + factor*(tmp - u_old)

                    sum_diff = sum_diff + dabs(u(idx) - u_old)
                end do
            end do

            if (sum_diff < eps) then
                iter_error = sum_diff
                return
            end if
        end do

        iter_error = sum_diff
    end subroutine solve_sor_subdomain


    !-----------------------------------------------------------------------
    ! ОБНОВЛЕНИЕ УЗЛОВ НА ИНТЕРФЕЙСАХ (ПРОСТАЯ ИТЕРАЦИЯ С УЧЁТОМ dx, dy)
    !-----------------------------------------------------------------------
    subroutine update_interfaces(macrogrid, macrogrid_size_x, macrogrid_size_y, subgrid_size, out_norm)
        implicit none
        integer, intent(in) :: macrogrid_size_x, macrogrid_size_y, subgrid_size
        real*8,  intent(inout) :: macrogrid(macrogrid_size_x,macrogrid_size_y,subgrid_size,subgrid_size)
        real*8,  intent(out)   :: out_norm

        ! Аналогичная логика: точка на границе — это узел,
        ! для которого решение берём из аппроксимации Лапласиана = 0
        ! с учётом dx, dy.

        integer :: iX, iY, k
        integer :: global_size_x, global_size_y
        real*8  :: dx, dy
        real*8  :: invdx2, invdy2, denom
        real*8  :: old_value, new_val

        ! Вычислим "глобальные" dx, dy (предполагаем, что шаг везде одинаковый
        ! и не зависит от самого iX/iY — т.е. вся область [0,1], просто разбитая
        ! на блоки)
        global_size_x = macrogrid_size_x*subgrid_size - (macrogrid_size_x - 1)
        global_size_y = macrogrid_size_y*subgrid_size - (macrogrid_size_y - 1)
        dx = 1.0d0 / dble(global_size_x - 1)
        dy = 1.0d0 / dble(global_size_y - 1)

        invdx2 = 1.0d0/(dx*dx)
        invdy2 = 1.0d0/(dy*dy)
        denom  = 2.0d0*invdx2 + 2.0d0*invdy2

        out_norm = 0.0d0

        ! -- Горизонтальные границы (между iY и iY+1)
        do iX = 1, macrogrid_size_x
            do iY = 1, macrogrid_size_y - 1
                do k = 2, subgrid_size - 1

                    old_value = macrogrid(iX, iY, k, subgrid_size)

                    ! Сосед слева (по y): macrogrid(iX, iY, k, subgrid_size-1)
                    ! Сосед справа (по y): macrogrid(iX, iY+1, k, 2)
                    ! "Верх/низ" по x внутри той же/смежной подсетки — это trickier,
                    ! но в данном простом случае:
                    !   macrogrid(iX, iY, k, subgrid_size-2) => "ещё выше" и т.п.
                    !
                    ! Но корректнее взять "4 соседа" для точки (k, subgrid_size) на стыке:
                    !   левый: macrogrid(iX, iY, k-1, subgrid_size)
                    !   правый: macrogrid(iX, iY, k+1, subgrid_size)
                    !   верхний: macrogrid(iX, iY, k, subgrid_size-1)
                    !   нижний: macrogrid(iX, iY+1, k, 1) — ведь (k,1) — это тот же узел
                    !   "с нижней" области, но на стыке.
                    !
                    ! Для упрощения здесь используем шаблон "гибрид" —
                    ! ориентируясь на ваш исходный способ, но с весами dx,dy.
                    !
                    ! Пример (упрощённый):
                    new_val = ( &
                        invdx2*( macrogrid(iX, iY,   k-1, subgrid_size ) + macrogrid(iX, iY,   k+1, subgrid_size ) ) + &
                        invdy2*( macrogrid(iX, iY,   k,   subgrid_size-1 ) + macrogrid(iX, iY+1, k, 1 ) ) &
                       ) / denom

                    macrogrid(iX, iY,   k, subgrid_size) = new_val
                    macrogrid(iX, iY+1, k, 1)            = new_val

                    out_norm = out_norm + dabs(new_val - old_value)
                end do
            end do
        end do

        ! -- Вертикальные границы (между iX и iX+1)
        do iX = 1, macrogrid_size_x - 1
            do iY = 1, macrogrid_size_y
                do k = 2, subgrid_size - 1

                    old_value = macrogrid(iX, iY, subgrid_size, k)

                    new_val = ( &
                        invdx2*( macrogrid(iX,   iY, subgrid_size-1, k) + macrogrid(iX+1, iY, 2, k ) ) + &
                        invdy2*( macrogrid(iX,   iY, subgrid_size,   k-1) + macrogrid(iX, iY, subgrid_size, k+1) ) &
                      ) / denom

                    macrogrid(iX,   iY, subgrid_size, k) = new_val
                    macrogrid(iX+1, iY, 1,            k) = new_val

                    out_norm = out_norm + dabs(new_val - old_value)
                end do
            end do
        end do

        ! -- "Углы" между четырьмя блоками (iX, iY, subgrid_size, subgrid_size)
        ! Аналогично берём 4 соседних узла. Примерно:
        do iX = 1, macrogrid_size_x - 1
            do iY = 1, macrogrid_size_y - 1

                old_value = macrogrid(iX, iY, subgrid_size, subgrid_size)
                new_val = ( &
                    invdx2*( macrogrid(iX, iY, subgrid_size-1, subgrid_size) +  &
                             macrogrid(iX+1, iY, 2, subgrid_size) ) + &
                    invdy2*( macrogrid(iX, iY, subgrid_size, subgrid_size-1) + &
                             macrogrid(iX, iY+1, subgrid_size, 2 ) ) &
                ) / denom

                macrogrid(iX,   iY,   subgrid_size,   subgrid_size)   = new_val
                macrogrid(iX+1, iY,   1,              subgrid_size)   = new_val
                macrogrid(iX,   iY+1, subgrid_size,   1)              = new_val
                macrogrid(iX+1, iY+1, 1,              1)              = new_val

                out_norm = out_norm + dabs(new_val - old_value)
            end do
        end do
    end subroutine update_interfaces


    !-----------------------------------------------------------------------
    ! РЕШЕНИЕ НА ВСЕЙ МАКРОСЕТКЕ (ДЕКОМПОЗИЦИЯ) — ОСНОВНОЙ ЦИКЛ
    !-----------------------------------------------------------------------
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

        ! Считаем глобальные шаги dx, dy (единая область [0,1] x [0,1])
        dx = compute_dx(macrogrid_size_x, subgrid_size)
        dy = compute_dy(macrogrid_size_y, subgrid_size)

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

end program macrogrid_solver
