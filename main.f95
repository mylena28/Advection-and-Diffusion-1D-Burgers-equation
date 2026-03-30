
! Título: Equação de Burgers por diferentes métodos
! Descrição: Coaparação entre solução analítica e numérica com Euler Explícito, Runge-Kutta
! Autor: Mylena Carvalho Silva

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

program burgers

    !Verificar variaveis
    implicit none

    integer, parameter :: dp = selected_real_kind(15, 307)

    !Declarar tipo de  variaveis 
    real(dp) :: pi, c, ni, dx, dt, CFL, t, tf
    integer :: n, q

    !Valor das constantes  
    c = 4.0   !m/s -- velocidade
    ni = 0.3  !m^2/s -- viscosidade
    n = 32    !Pontos de divisão = 32, 64, 128, 256
    CFL = 0.1_dp     
    pi = 4.0*atan(1.0)  !Criando pi 3.1415927410125732
    tf = pi/4    !Tempo = 0.78
    dx = (2*pi)/(n-1)  !Variação de x para n pontos
    dt = CFL*dx/c      
    t = 0.0

    do while (t <= tf)
        call analitica(c, ni, n, t, pi, dx, q)
        t = t + dt
        q = q + 1
    end do

    call explicita(ni, n, int(tf/dt), dt, dx)
end program burgers
