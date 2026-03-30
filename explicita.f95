!==============================================================
! Author: Mylena Carvalho Silva
! Date: 2023
!
! Problem:
! 1D Burgers' equation (nonlinear advection-diffusion):
!     du/dt + u du/dx = nu d²u/dx²
!
! Methods:
! - Analytical solution (validation)
! - Explicit finite difference scheme
!
! Dependencies:
! - analitica.f95
! - explicita.f95
! - main.f95
!
! Notes:
! Developed for academic and research purposes in fluid mechanics
! and numerical methods.
!==============================================================

subroutine explicita(ni, nx, nt, dt, dx)
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)
    real(dp) :: ni, dt, dx
    real(dp) :: A, B
    real(dp), allocatable, dimension(:) :: x, u, ut
    integer :: i, nx, nt, t, ierr
    character(len=50) :: nomearq, arqler

    !Determinar tamanho do vetor
    allocate(x(nx), u(nx), ut(nx))

    !Lendo dados de entrada para o t = 0
    if (nx < 100) then
        write (arqler, '("src/N0", I2,"/analitica/analitica_000.txt")') nx
        write (nomearq, '("src/N0", I2,"/explicita/exp_000.txt")') nx
    else
        write (arqler, '("src/N", I3,"/analitica/analitica_000.txt")') nx
        write (nomearq, '("src/N", I3,"/explicita/exp_000.txt")') nx
    end if
    call ler(arqler, x, u, nx)
    call escrever(nomearq, x, u, nx)

    A = dt/(8*dx)
    B = ni*dt/(dx*dx)

    !Calcular os proximos passos
    do t = 1, nt
        ut(1) = u(1) - A*(((u(2)+u(1))*(u(2)+u(1))) &
                          -((u(1)+u(nx))*(u(1)+u(nx)))) &
                          + B*(u(2)-2*u(1)+u(nx))
        ut(nx) = u(nx) - A*(((u(1)+u(nx))*(u(1)+u(nx))) &
                          -((u(nx)+u(nx-1))*(u(nx)+u(nx-1)))) &
                          + B*(u(1)-2*u(nx)+u(nx-1))
        do i = 2, (nx-1)
            ut(i) = u(i) - A*(((u(i+1)+u(i))*(u(i+1)+u(i))) &
                            -((u(i)+u(i-1))*(u(i)+u(i-1)))) &
                            + B*(u(i+1)-2*u(i)+u(i-1))
        end do

    if (nx < 100) then
        if (t < 10) then
            write (nomearq, '("src/N0", I2,"/explicita/exp_00", I1, ".txt")') nx, t
        elseif (t < 100) then
            write (nomearq, '("src/N0", I2,"/explicita/exp_0", I2, ".txt")')  nx, t
        else
            write (nomearq, '("src/N0", I2,"/explicita/exp_", I3, ".txt")')  nx, t
        end if
    else
        if (t < 10) then
            write (nomearq, '("src/N", I3,"/explicita/exp_00", I1, ".txt")') nx, t
        elseif (t < 100) then
            write (nomearq, '("src/N", I3,"/explicita/exp_0", I2, ".txt")')  nx, t
        else
            write (nomearq, '("src/N", I3,"/explicita/exp_", I3, ".txt")')  nx, t
        end if
    end if
    call escrever(nomearq, x, ut, nx)
    u = ut
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Funções!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    contains

    ! Ler arquivo txt e transformar em vetor
    subroutine ler(nome1, x1, x2, n)

        character(len=37), intent(in) :: nome1
        integer, intent(in) :: n
        real(dp), dimension(:), intent(out) :: x1(n), x2(n)
        integer :: i

        open(1, file=nome1, status="old", action="read")
        do i = 1, n
            read(1, *) x1(i), x2(i)
        end do
        close(1)

    end subroutine ler

    ! Criar arquivo com dados
    subroutine escrever(nome2, x1, x2, n)

        character(len=30), intent(in) :: nome2
        integer, intent(in) :: n
        real(dp), dimension(:), intent(in) :: x1(n), x2(n)
        integer :: i

        open(1, file=nome2, status="replace", action="write")
        ! Escrever os dados
        do i = 1, n
            write(1, *) x1(i), x2(i)
        end do
        close(1)
        return
    end subroutine escrever

end subroutine explicita
