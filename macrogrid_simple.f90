program macrogrid_solver
    implicit none
    
    ! Логический параметр, нужно ли накапливать ошибку от всех подобластей
    ! при вычислении нормы на уровне макроитераций
    logical :: accumulate_subgrid_error = .true.
    
    ! В данном примере зададим некий верхний предел на размеры
    integer, parameter :: MAX_SUBGRID_SIZE = 2048
    integer, parameter :: MAX_NUM_SUBGRIDS_X = 20
    integer, parameter :: MAX_NUM_SUBGRIDS_Y = 20
    
    ! Разные точности и максимальные итерации
    integer, parameter :: max_iter_subgrid    = 10000   ! Макс. итераций SOR в подобласти
    integer, parameter :: max_iter_interface  = 5000    ! Макс. итераций итерации по интерфейсу
    real*8,  parameter :: eps_subgrid         = 1.0d-8  ! Точность для SOR в подобласти
    real*8,  parameter :: eps_interface       = 1.0d-5  ! Точность для итераций по интерфейсу
    
    ! Для примера зададим количество подсеток по x и y
    ! и размер каждой подсетки (в узлах)
    integer :: Nx, Ny, subgrid_size
    
    ! Параметр релаксации
    real*8  :: omega
    
    ! Трёхмерный массив (Nx,Ny) и каждая ячейка содержит subgrid_size x subgrid_size
    ! Чтобы хранить максимум, возьмём верхние границы
    real*8 :: macrogrid(MAX_NUM_SUBGRIDS_X, MAX_NUM_SUBGRIDS_Y, &
                        MAX_SUBGRID_SIZE,     MAX_SUBGRID_SIZE)
                        
    integer :: iterations
    
    ! Ниже — просто пример инициализации:
    Nx           = 5      ! Число подсеток по оси X
    Ny           = 3      ! Число подсеток по оси Y
    subgrid_size = 16     ! Размер каждой квадратной подобласти (кол-во узлов)
    omega        = 1.9d0  ! Параметр релаксации для SOR
    
    print *, "Пример: прямоугольная макросетка ", Nx, "x", Ny, &
             " c сабгридами ", subgrid_size, "x", subgrid_size
    
    ! Инициализируем макросетку для простейшей граничной задачи (единицы по периметру)
    call initialize_macrogrid_ones(macrogrid, Nx, Ny, subgrid_size)
    
    ! Запустим решение методом декомпозиции с заданными параметрами:
    call measure_time( accumulate_subgrid_error,               &
                       macrogrid, Nx, Ny, subgrid_size,        &
                       omega,                                  &
                       max_iter_subgrid, max_iter_interface,   &
                       eps_subgrid, eps_interface )
    
    ! Проверим погрешность (если у нас задача с граничным условием "1" по краям):
    call compute_error_ones(macrogrid, Nx, Ny, subgrid_size)
    
    ! Аналогичный пример для "логарифмических" граничных условий:
    call initialize_macrogrid_logs(macrogrid, Nx, Ny, subgrid_size)
    
    ! Снова решаем
    call measure_time( accumulate_subgrid_error,               &
                       macrogrid, Nx, Ny, subgrid_size,        &
                       omega,                                  &
                       max_iter_subgrid, max_iter_interface,   &
                       eps_subgrid, eps_interface )
    
    ! И считаем ошибку с логарифмическим решением
    call compute_error_logs(macrogrid, Nx, Ny, subgrid_size)
    
end program macrogrid_solver


!-----------------------------------------------------------------------
! ИНИЦИАЛИЗАЦИЯ: "ЕДИНИЧКИ" НА ГРАНИЦАХ
!-----------------------------------------------------------------------
subroutine initialize_macrogrid_ones(macrogrid, num_subgrids_x, num_subgrids_y, subgrid_size)
    implicit none
    integer, intent(in)            :: num_subgrids_x, num_subgrids_y, subgrid_size
    real*8,  intent(inout)         :: macrogrid(:,:,:,:)
    
    integer :: global_size_x, global_size_y
    integer :: iX, iY, lX, lY, i1, j1
    
    ! Число узлов в глобальной сетке по x и y
    global_size_x = num_subgrids_x * subgrid_size - (num_subgrids_x - 1)
    global_size_y = num_subgrids_y * subgrid_size - (num_subgrids_y - 1)
    
    do iX = 1, num_subgrids_x
        do iY = 1, num_subgrids_y
            do i1 = 1, subgrid_size
                do j1 = 1, subgrid_size
                    ! Посчитаем глобальные индексы
                    lX = i1 + (iX - 1)*subgrid_size - (iX - 1)
                    lY = j1 + (iY - 1)*subgrid_size - (iY - 1)
                    
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
subroutine initialize_macrogrid_logs(macrogrid, num_subgrids_x, num_subgrids_y, subgrid_size)
    implicit none
    real*8, parameter :: R1 = 0.1d0, R2 = 1.0d0
    real*8, parameter :: x_min = 0.3d0, y_min = 0.0d0
    real*8, parameter :: length_domain = 0.4d0
    
    integer, intent(in)            :: num_subgrids_x, num_subgrids_y, subgrid_size
    real*8,  intent(inout)         :: macrogrid(:,:,:,:)
    
    integer :: global_size_x, global_size_y
    integer :: iX, iY, lX, lY, i1, j1
    
    real*8 :: boundary_value
    
    global_size_x = num_subgrids_x * subgrid_size - (num_subgrids_x - 1)
    global_size_y = num_subgrids_y * subgrid_size - (num_subgrids_y - 1)
    
    do iX = 1, num_subgrids_x
        do iY = 1, num_subgrids_y
            do i1 = 1, subgrid_size
                do j1 = 1, subgrid_size
                    lX = i1 + (iX - 1)*subgrid_size - (iX - 1)
                    lY = j1 + (iY - 1)*subgrid_size - (iY - 1)
                    
                    if (lX == 1 .or. lX == global_size_x .or. &
                        lY == 1 .or. lY == global_size_y) then
                        
                        ! Значение точного решения на границе
                        boundary_value = log( sqrt( (x_min + length_domain*(lX-1)/(global_size_x-1))**2 + &
                                                    (y_min + length_domain*(lY-1)/(global_size_y-1))**2 ) * R2 / (R1*R1) ) / &
                                         log(R2/R1)
                        
                        macrogrid(iX,iY,i1,j1) = boundary_value
                    else
                        macrogrid(iX,iY,i1,j1) = 0.0d0
                    end if
                end do
            end do
        end do
    end do
end subroutine initialize_macrogrid_logs


!-----------------------------------------------------------------------
! ВЫЧИСЛЕНИЕ ОШИБКИ ДЛЯ "ЕДИНИЧКИ" НА ГРАНИЦЕ
!-----------------------------------------------------------------------
subroutine compute_error_ones(macrogrid, num_subgrids_x, num_subgrids_y, subgrid_size)
    implicit none
    integer, intent(in)    :: num_subgrids_x, num_subgrids_y, subgrid_size
    real*8, intent(in)     :: macrogrid(:,:,:,:)
    
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
    
    print *, 'Maximum error (ones-boundary) = ', max_error
end subroutine compute_error_ones


!-----------------------------------------------------------------------
! ВЫЧИСЛЕНИЕ ОШИБКИ ДЛЯ "ЛОГАРИФМИЧЕСКИХ" ГРАНИЧНЫХ УСЛОВИЙ
!-----------------------------------------------------------------------
subroutine compute_error_logs(macrogrid, num_subgrids_x, num_subgrids_y, subgrid_size)
    implicit none
    real*8, parameter :: R1 = 0.1d0, R2 = 1.0d0
    real*8, parameter :: x_min = 0.3d0, y_min = 0.0d0
    real*8, parameter :: length_domain = 0.4d0
    
    integer, intent(in)    :: num_subgrids_x, num_subgrids_y, subgrid_size
    real*8, intent(in)     :: macrogrid(:,:,:,:)
    
    integer :: global_size_x, global_size_y
    integer :: iX, iY, i1, j1, lX, lY
    real*8  :: exact_val, err, max_error
    
    global_size_x = num_subgrids_x * subgrid_size - (num_subgrids_x - 1)
    global_size_y = num_subgrids_y * subgrid_size - (num_subgrids_y - 1)
    
    max_error = 0.0d0
    
    do iX = 1, num_subgrids_x
        do iY = 1, num_subgrids_y
            do i1 = 1, subgrid_size
                do j1 = 1, subgrid_size
                    lX = i1 + (iX - 1)*subgrid_size - (iX - 1)
                    lY = j1 + (iY - 1)*subgrid_size - (iY - 1)
                    
                    exact_val = log( sqrt( (x_min + length_domain*(lX-1)/(global_size_x-1))**2 + &
                                           (y_min + length_domain*(lY-1)/(global_size_y-1))**2 ) * R2 / (R1*R1) ) / &
                                log(R2/R1)
                    
                    err = abs( macrogrid(iX,iY,i1,j1) - exact_val )
                    if (err > max_error) max_error = err
                end do
            end do
        end do
    end do
    
    print *, 'Maximum error (log-boundary) = ', max_error
end subroutine compute_error_logs


!-----------------------------------------------------------------------
! ПРОЦЕДУРА МЕТОДА SOR В ОДНОЙ ПОДОБЛАСТИ
!-----------------------------------------------------------------------
subroutine solve_sor_subdomain(u, n, omega, eps_sor, max_iter_sor, out_error)
    ! Решает простейшую 2D-задачу (лапласиан) на квадрате n x n
    ! u — массив размером n*n (развёрнут в 1D), границы фиксированы
    ! omega — параметр релаксации
    ! eps_sor — требуемая точность для SOR
    ! max_iter_sor — максимум итераций
    ! out_error — накопленная разница за последнюю итерацию
    implicit none
    
    integer, intent(in) :: n, max_iter_sor
    real*8,  intent(in) :: omega, eps_sor
    real*8,  intent(inout) :: u(n*n)
    real*8,  intent(out) :: out_error
    
    integer :: iter, i, ix, iy, idx
    real*8  :: old_val, sum_diff
    
    do iter = 1, max_iter_sor
        sum_diff = 0.0d0
        
        ! Простой шаблон: (ix,iy) – внутренние узлы
        do iy = 2, n-1
            do ix = 2, n-1
                idx     = ix + (iy-1)*n
                old_val = u(idx)
                ! классическая формула SOR для 5-точечного лапласиана
                u(idx)  = old_val + omega*0.25d0*(u(idx-1) + u(idx+1) + &
                                                 u(idx-n) + u(idx+n) - 4*old_val)
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
! ОБНОВЛЕНИЕ УЗЛОВ НА ГРАНИЦАХ ПОДОБЛАСТЕЙ (ПРОСТАЯ ИТЕРАЦИЯ)
!-----------------------------------------------------------------------
subroutine update_boundaries( macrogrid, num_subgrids_x, num_subgrids_y, subgrid_size, &
                              out_norm )
    implicit none
    
    integer, intent(in) :: num_subgrids_x, num_subgrids_y, subgrid_size
    real*8,  intent(inout) :: macrogrid(:,:,:,:)
    real*8,  intent(out) :: out_norm
    
    integer :: iX, iY, k
    real*8  :: old_value, avg
    
    out_norm = 0.0d0
    
    ! Горизонтальные границы между subgrid(iX,iY) и subgrid(iX,iY+1)
    do iX = 1, num_subgrids_x
        do iY = 1, num_subgrids_y-1
            do k = 2, subgrid_size-1
                old_value = macrogrid(iX, iY, k, subgrid_size)
                ! Простейшее "усреднение" (как и было в старом коде)
                avg = ( macrogrid(iX,iY,   k, subgrid_size-1)*4.d0 + &
                        macrogrid(iX,iY+1, k, 2)*4.d0            - &
                        macrogrid(iX,iY,   k, subgrid_size-2)   - &
                        macrogrid(iX,iY+1, k, 3 ) )/6.d0
                
                macrogrid(iX,iY,   k, subgrid_size) = avg
                macrogrid(iX,iY+1, k, 1)            = avg
                
                out_norm = out_norm + dabs(avg - old_value)
            end do
        end do
    end do
    
    ! Вертикальные границы между subgrid(iX,iY) и subgrid(iX+1,iY)
    do iX = 1, num_subgrids_x-1
        do iY = 1, num_subgrids_y
            do k = 2, subgrid_size-1
                old_value = macrogrid(iX,   iY, subgrid_size, k)
                avg       = ( 4.d0*macrogrid(iX,   iY, subgrid_size-1,k) + &
                              4.d0*macrogrid(iX+1, iY, 2,k)               - &
                              macrogrid(iX,   iY, subgrid_size-2,k)       - &
                              macrogrid(iX+1, iY, 3,k) ) / 6.d0
                
                macrogrid(iX,   iY, subgrid_size, k) = avg
                macrogrid(iX+1, iY, 1,            k) = avg
                out_norm = out_norm + dabs(avg - old_value)
            end do
        end do
    end do
    
    ! Узлы в "углах" между четырьмя блоками
    do iX = 1, num_subgrids_x-1
        do iY = 1, num_subgrids_y-1
            old_value = macrogrid(iX,iY, subgrid_size, subgrid_size)
            avg       = ( 4.d0 * macrogrid(iX,iY,   subgrid_size,   subgrid_size-1) + &
                          4.d0 * macrogrid(iX,iY,   subgrid_size-1, subgrid_size)   + &
                          4.d0 * macrogrid(iX+1,iY+1, 2, 1)                          + &
                          4.d0 * macrogrid(iX+1,iY+1, 1, 2 )                          - &
                          macrogrid(iX,iY,   subgrid_size,   subgrid_size-2)         - &
                          macrogrid(iX,iY,   subgrid_size-2, subgrid_size)           - &
                          macrogrid(iX+1,iY+1, 3, 1 )                                 - &
                          macrogrid(iX+1,iY+1, 1, 3 ) ) / 12.d0
            
            macrogrid(iX,   iY,   subgrid_size,   subgrid_size)   = avg
            macrogrid(iX+1, iY,   1,              subgrid_size)   = avg
            macrogrid(iX,   iY+1, subgrid_size,   1)              = avg
            macrogrid(iX+1, iY+1, 1,              1)              = avg
            
            out_norm = out_norm + dabs(avg - old_value)
        end do
    end do
    
end subroutine update_boundaries


!-----------------------------------------------------------------------
! РЕШЕНИЕ НА ВСЕЙ МАКРОСЕТКЕ (ДЕКОМПОЗИЦИЯ)
!-----------------------------------------------------------------------
subroutine solve_macrogrid( accumulate_subgrid_error,                        &
                            macrogrid, num_subgrids_x, num_subgrids_y,       &
                            subgrid_size, omega,                              &
                            max_iter_sor, max_iter_interface,                &
                            eps_sor, eps_interface,                           &
                            out_iterations )
    implicit none
    
    logical, intent(in) :: accumulate_subgrid_error
    integer, intent(in) :: num_subgrids_x, num_subgrids_y
    integer, intent(in) :: subgrid_size, max_iter_sor, max_iter_interface
    real*8,  intent(in) :: omega, eps_sor, eps_interface
    real*8,  intent(inout) :: macrogrid(:,:,:,:)
    integer, intent(out) :: out_iterations
    
    integer :: iter, iX, iY
    real*8  :: norm_total, norm_sub_sor, norm_bound
    
    ! Развёртка поддиапазона размером subgrid_size * subgrid_size
    real*8, allocatable :: local_u(:)
    
    out_iterations = 0
    
    do iter = 1, max_iter_interface
        
        norm_total = 0.0d0
        
        ! 1) Запускаем (возможно) SOR в каждой подобласти
        do iX = 1, num_subgrids_x
            do iY = 1, num_subgrids_y
                ! Собираем подблок в local_u (или работаем напрямую)
                ! Здесь для примера работаем напрямую с pointer/transfer, упрощённо
                ! Но ниже показан вариант просто через alloc:
                allocate(local_u(subgrid_size*subgrid_size))
                call copy_in_subgrid( macrogrid(iX,iY,:,:), local_u, subgrid_size )
                
                ! Решаем методом SOR в локальной области
                call solve_sor_subdomain( local_u, subgrid_size, omega, &
                                          eps_sor, max_iter_sor, norm_sub_sor )
                
                ! Накопим ошибку
                norm_total = norm_total + norm_sub_sor
                
                ! Запишем решение обратно
                call copy_out_subgrid( macrogrid(iX,iY,:,:), local_u, subgrid_size )
                
                deallocate(local_u)
            end do
        end do
        
        ! Если НЕ накапливаем ошибку по подобластям — обнуляем norm_total
        if (.not.accumulate_subgrid_error) norm_total = 0.0d0
        
        ! 2) Обновляем значения на границах подобластей (простая итерация)
        call update_boundaries(macrogrid, num_subgrids_x, num_subgrids_y, subgrid_size, norm_bound)
        
        norm_total = norm_total + norm_bound
        
        ! Проверяем, достигли ли требуемой точности
        if (norm_total < eps_interface) then
            out_iterations = iter
            return
        end if
        
    end do
    
    out_iterations = max_iter_interface
end subroutine solve_macrogrid


!-----------------------------------------------------------------------
! КОПИРОВАНИЕ ИЗ 4D-СЕГМЕНТА В 1D-МАССИВ
!-----------------------------------------------------------------------
subroutine copy_in_subgrid(block_2d, local_u, n)
    implicit none
    integer, intent(in)            :: n
    real*8,  intent(in)            :: block_2d(n,n)
    real*8,  intent(out)           :: local_u(n*n)
    
    integer :: ix, iy, idx
    
    do iy = 1, n
        do ix = 1, n
            idx           = ix + (iy-1)*n
            local_u(idx)  = block_2d(ix, iy)
        end do
    end do
end subroutine copy_in_subgrid


!-----------------------------------------------------------------------
! КОПИРОВАНИЕ ОБРАТНО ИЗ 1D МАССИВА В 4D-БЛОК
!-----------------------------------------------------------------------
subroutine copy_out_subgrid(block_2d, local_u, n)
    implicit none
    integer, intent(in)           :: n
    real*8,  intent(in)           :: local_u(n*n)
    real*8,  intent(inout)        :: block_2d(n,n)
    
    integer :: ix, iy, idx
    
    do iy = 1, n
        do ix = 1, n
            idx             = ix + (iy-1)*n
            block_2d(ix,iy) = local_u(idx)
        end do
    end do
end subroutine copy_out_subgrid


!-----------------------------------------------------------------------
! ПРОСТАЯ РАСЧЁТ ВРЕМЕНИ
!-----------------------------------------------------------------------
subroutine measure_time( accumulate_subgrid_error,                        &
                         macrogrid, num_subgrids_x, num_subgrids_y,       &
                         subgrid_size, omega,                             &
                         max_iter_sor, max_iter_interface,                &
                         eps_sor, eps_interface )
    implicit none
    logical, intent(in) :: accumulate_subgrid_error
    integer, intent(in) :: num_subgrids_x, num_subgrids_y
    integer, intent(in) :: subgrid_size, max_iter_sor, max_iter_interface
    real*8,  intent(in) :: omega, eps_sor, eps_interface
    real*8,  intent(inout) :: macrogrid(:,:,:,:)
    
    real*8 :: t1, t2
    integer :: iterations
    
    call cpu_time(t1)
    
    call solve_macrogrid( accumulate_subgrid_error,               &
                          macrogrid, num_subgrids_x, num_subgrids_y, &
                          subgrid_size, omega,                      &
                          max_iter_sor, max_iter_interface,         &
                          eps_sor, eps_interface,                   &
                          iterations )
    
    call cpu_time(t2)
    
    print *, 'Time for decomposition = ', (t2 - t1), ' seconds.'
    print *, 'Number of interface iterations = ', iterations
end subroutine measure_time
