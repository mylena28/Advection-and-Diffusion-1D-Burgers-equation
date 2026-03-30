!==============================================================
! File: main.f95
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
!
! Notes:
! Developed for academic and research purposes in fluid mechanics
! and numerical methods.
!==============================================================


subroutine analitica(c, ni, n, t, pi, dx, q)
    implicit none

    integer, parameter :: dp = selected_real_kind(15, 307)
    real(dp) :: pi, c, ni, dx, cont, velocidade, phi_derivada, soma_phi_derivada, phi, soma_phi, t
    real(dp), allocatable, dimension(:) :: x,u
    integer :: m, i, b, n, q
    character(len=50) :: nomearq

    !Determinar tamanho do vetor
    allocate(x(n), u(n))

    !Contadores
    cont = -pi
    b = 3     !Pontos da somatoria

   !valor das variaveis
   soma_phi_derivada = 0.0
   soma_phi = 0.0

    if (n<100) then
        if (q < 10) then
            write (nomearq, '("src/N0", I2,"/analitica/analitica_00", I1, ".txt")' ) n, q
        elseif (q < 100) then
            write (nomearq, '("src/N0", I2,"/analitica/analitica_0", I2, ".txt")' ) n, q
        else
            write (nomearq, '("src/N0", I2,"/analitica/analitica_", I3, ".txt")' ) n, q
        end if
    else
        if (q < 10) then
            write (nomearq, '("src/N", I3,"/analitica/analitica_00", I1, ".txt")' ) n, q
        elseif (q < 100) then
            write (nomearq, '("src/N", I3,"/analitica/analitica_0", I2, ".txt")' ) n, q
        else
            write (nomearq, '("src/N", I3,"/analitica/analitica_", I3, ".txt")' ) n, q
        end if
    end if

    ! Abra o arquivo
    open(unit=10, file=trim(nomearq), status="replace")

    !Calculando o valor da velocidade 
    do i = 1,n

         !Criando vetor x
         !Cada valor de i, dentro de l, em um intervalo de dx é acrescentado ao vetor x
         x(i) = cont
         cont = cont+dx

        !Loop com velocidade constnte para somatoria de n pontos
        do m = -b,b

            !Somatorio da derivada
            phi_derivada = -(((x(i)-c*t)-(2*m+1)*pi)/(2*ni*(t+1))) &
                            *exp(-(((x(i)-c*t)-(2*m+1)*pi)**2)/(4*ni*(t+1)))
            soma_phi_derivada = soma_phi_derivada + phi_derivada

            !Somatorio de phi
            phi = exp(-(((x(i)-c*t)-(2*m+1)*pi)**2)/(4*ni*(t+1)))
            soma_phi = soma_phi + phi

            !Limpando a variavel
            phi_derivada = 0.0
            phi = 0.0

        end do

        !Velocidade para o tempo 0
        velocidade = c - 2 * ni * soma_phi_derivada / soma_phi
        u(i) = velocidade

        !Arquivo txt
        write(10,'(5E12.4)') x(i),u(i)

        !Limpando contador
        soma_phi_derivada = 0.0
        soma_phi= 0.0

    end do
    close(10)

end subroutine analitica
