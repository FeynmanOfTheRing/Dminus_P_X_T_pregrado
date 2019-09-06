! ============================================================================
! Name        : D_menos_con_B.f90
! Author      : 
! Version     :
! Copyright   : Your copyright notice
! Description : Hello World in Fortran
! ============================================================================

module subrutinas_funciones_parametros

    implicit none
    real(8),parameter ::PI=4*atan(1.),d=0.00001d0, R1=1, R2=1,dseta=1d5, GAM = 0,h=0.06 !0.06 radianes !parametros del problema
    real(8):: mu,Er,En,En_ant,Uee,Ue1d,Ue2d,phi,tol,Fdospi,Fcero,F,Fant,yini!masa reducida, energ�as, potenciales,to
    real(8),dimension(3):: vector_entra,vector_sale!vectores que utilizo para hacer RK4

    integer :: i,j,k,bflag

contains
    subroutine potential(phi,Uee)
        real(8), intent(in)  :: phi
        real(8), intent(out) :: Uee
        Uee=(2.)/sqrt(R1**2 + R2**2 -2*R1*R2*dcos(phi) +d**2) !Potencial Uee electr�n electr�n

    end subroutine
    subroutine RUNGE_KUTTA(xt,xth) !subrutina RK4
        implicit none
        real(8), dimension (3), intent(IN) :: xt
        real(8), dimension (3), intent(OUT) :: xth
        real(8), dimension (3) :: F1,F2_in,F2,F3_in,F3,F4,F4_in, suma

        call EDO(xt,F1)
        F2_in= xt + F1/2d0
        call EDO(F2_in,F2)
        F3_in= xt + F2/2d0
        call EDO(F3_in,F3)
        F4_in=xt+F3
        call EDO(F4_in,F4)
        xth= xt + (F1+ 2d0*F2 +2d0*F3+ F4)/6d0
    end subroutine
    subroutine EDO(vector_in,vector_out)!campo vectorial del problema para realizar RK4

        implicit none
        real(8), dimension(3), intent(IN) :: vector_in
        real(8), dimension(3), intent(OUT) :: vector_out


        call potential(vector_in(1),Uee)


        vector_out(1)=1;
        vector_out(2)=vector_in(3)
        vector_out(3)= (Uee+Ue1d+Ue2d-En)*(1/mu)*vector_in(2)
        vector_out=h*vector_out
    end subroutine
    !***************Biseccion***************
    subroutine biseccion(xi,xd,fxi,fxd,xr)

        real(8), INTENT(INOUT) :: xi,xd
        REAL(8), INTENT(IN)    ::fxi,fxd
        REAL(8):: fxr
        REAL(8), INTENT(OUT)   ::xr


        xr=(xi+xd)/2d0

        En=xr

        call Rungekutta_anillo

        fxr=Fcero-Fdospi

        if(fxi*fxr.lt.0)then
            xd=xr
        else if(fxi*fxr.gt.0) then
            xi=xr
        end if

        xr=(xi+xd)/2d0
    end subroutine
    subroutine Rungekutta_anillo
        vector_entra(1)=0.0         !condiciones iniciales
        vector_entra(2)=yini
        vector_entra(3)=0.000001
        DO WHILE(vector_entra(1).LE.2*PI) !ciclo para hacer RK4, para cada energ�a de la iteraci�n correspondiente


            CALL RUNGE_KUTTA(vector_entra,vector_sale)
            vector_entra=vector_sale

            if(vector_sale(1) .le. (2*pi)) then !valor de la funci�n evaluada en 2pi
                Fdospi=vector_sale(2)
            end if
            if(k.eq.1) then
                write(3,*) Uee+Ue1d+Ue2d , vector_entra(1) !escribo los valores del potencial vs phi
            end if
        END DO
    end subroutine
end module

program D_menos_con_B

    use subrutinas_funciones_parametros

    implicit none
    real(8)::Estep,Efin,Eini
    integer(8)::ndata

    mu = (1./(R1)**2)+(1./(R2)**2)  !masa reducida
    Eini=1                         !energ�as inicial y final
    Efin=11
    ndata=int(Efin-Eini)*1d4
    write(*,*) ndata          !total de energ�as
    Estep=(Efin-Eini)/real(ndata)  !paso
    Ue1d = -2./sqrt(R1**2+dseta**2)!potenciales con la impureza
    Ue2d = -2./sqrt(R2**2+(d-dseta)**2)
    tol=1d-9              !tolerancia que funciona bisecci�n /18/01/19 = 1d-4
    vector_sale=0                  !invicio el vector que sale
    yini=0.1 !0.00001 funciona
    Fcero=yini
    open(unit=1,file ="F")              !aqu� guardo la funci�n para graficar
    Open(unit=2, file="energias.txt")   !aqu� guardo las energ�as
    open(unit=3, file="potencial.txt")  !aqu� guardo el potencial una vez obtenida la funci�n
    En=Eini
    k=0
    j=0
    Do while(En.le.Efin) !Ciclo que varia las energ�as

        call Rungekutta_anillo


        F=(Fcero-Fdospi)
        !


        !                   write(1,*) F
        !                   write(*,*) F
        if((F*Fant.lt.0d0)) then
            j=j+1

            call biseccion(En_ant,En,Fant,F,En)
            call Rungekutta_anillo
            if(abs(Fcero-Fdospi) .lt. tol) then

!                 k=k+1
!                write(*,*) "Er= ", ((GAM**2)/4d0)*(R1**2+R2**2) -real(k)*GAM+ (real(k)**2)/(R1**2+R2**2) + En
!                write(2,*) ((GAM**2)/4d0)*(R1**2+R2**2) -real(k)*GAM+ (real(k)**2)/(R1**2+R2**2) + En


                write(*,*) "Er= ", En
                write(2,*) En

                Fant=F
                En_ant=En
                j=0
                En=En+Estep
            end if

            if(j.gt.1d7)then
                write(*,*) "la biseccion no converge satisfactoriamente"
                write(*,*) "un aproximado del autovalor es: ", En

                En=En+Estep
                j=0
            end if




        else
            Fant=F
            En_ant=En
            En=En+Estep
            write(1,*) En
            bflag = 0
        end if







    End do
    Close(1)
    Close(2)
    close(3)


end program
