! ============================================================================
! Name        : D_minus_P_T_X.f90
! Author      : Richard Feynman
! Version     :
! Copyright   : Your copyright notice
! Description : Hello World in Fortran
! ============================================================================
module subrutinas_funciones_parametros

    implicit none
    real(8),parameter ::PI=4*atan(1.), GAM = 0,h=0.06,m_o=9.1d-31!0.06 radianes !parametros del problema
    real(8),parameter :: S11=1.16d-3,S12=-3.7d-4
    real(8):: d=0.000001d0, R1=1, R2=1,dseta=1d10
    real(8):: mu,Er,En,En_ant,Uee,Ue1d,Ue2d,phi,tol,Fdospi,Fcero,F,Fant,yini!masa reducida, energías, potenciales,to
    real(8)::X,P,Pfin,Pini,Pstep,T
    real(8),dimension(3):: vector_entra,vector_sale!vectores que utilizo para hacer RK4

    integer :: i,j,k,bandera_imprimir

contains
    real function rho(P)
        real(8)::rho_o
        real(8),intent(in)::P
        rho=sqrt((1d0-2d0*(S11+2d0*S12)*P))
        return
    end function
    real function H_p(P)
        real(8)::Ho
        real(8),intent(in)::P
        H_p=(1d0-(S11+2d0*S12)*P)
        return
    end function
    subroutine potential(phi,Uee)
        real(8), intent(in)  :: phi
        real(8), intent(out) :: Uee
        Uee=(2.)/sqrt(R1**2 + R2**2 -2*R1*R2*dcos(phi) +d**2) !Potencial Uee electrón electrón

    end subroutine

    real function m(X,P,T)
        implicit none
        real(8), intent(in) :: X,P,T
        m=m_o/(1d0+(PI_2(X)/3d0)*( (2d0/E_g_gm(X,P,T))+( 1d0 / (E_g_gm(X,P,T)+delta_o(X)) ) ) +delta(X) )
        RETURN
    end function
    real function PI_2(X)
        implicit none
        real(8),intent(in) ::X
        PI_2=28900d0-6290d0*X
        return
    end function
    real function delta_o(X)
        implicit none
        real(8),intent(in) ::X
        delta_o=341d0-66d0*X
        return
    end function

    real function delta(X)
        implicit none
        real(8),intent(in) ::X
        delta=-3.935d0+0.488d0*X+4.938d0*(X**2)
        return
    end function
    real function E_g_gm(X,P,T)
        implicit none
        real(8),intent(in) ::X,P,T
        E_g_gm=1519.4d0 + 1360d0*X+220d0*(X**2)+10.7d0*P - (0.5405d0*(T**2)/(204d0+T))
        return
    end function
    real function epsil(P,T)
        implicit none
        real(8), intent(in)::P,T

        if (T .le. 200) then
            epsil=12.74d0*exp((9.4d-5)*T-(1.67d-3)*P)
            return
        else if(T.gt.200) then
            epsil=12.29d0*exp((20.4d-5)*T-(1.73d-3)*P)
            return
        end if
    end function

    real function f_(X,P,T)
        implicit none
        real(8), intent(in)::X,P,T

        f_=m(X,P,T)/m(0d0,0d0,4d0)
        return
    end function

    real function g_(P,T)
        implicit none
        real(8), intent(in)::P,T

        g_ =epsil(P,T)/epsil(0d0,4d0)
        return
    end function
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


        vector_out(1)=1
        vector_out(2)=vector_in(3)
        vector_out(3)= ( f_(X,P,T)/((-1d0)*mu) )*( En-(Uee+Ue1d+Ue2d)/g_(P,T) )*vector_in(2)
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
        bandera_imprimir=0
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
        vector_entra(1)=-2d0*PI      !condiciones iniciales
        vector_entra(2)=yini
        vector_entra(3)=0.01d0
        DO WHILE(vector_entra(1).LE.(2d0*PI)) !ciclo para hacer RK4, para cada energía de la iteración correspondiente


            CALL RUNGE_KUTTA(vector_entra,vector_sale)
            vector_entra=vector_sale

            !valor de la función evaluada en 2pi
            Fdospi=vector_sale(2)

            if(bandera_imprimir.eq.1) then

                write(3,*) Uee+Ue1d+Ue2d ,";", vector_entra(1) !escribo los valores del potencial vs phi
                write(1,*) vector_entra(2),";", vector_entra(1)
                write(2,*) (vector_entra(2))**2,";",vector_entra(1)
            end if

        END DO
    end subroutine
end module


program D_minus_P_T_X

    use subrutinas_funciones_parametros

    implicit none
    real(8)::Estep,Efin,Eini
    integer(8)::ndata

    X=0d0
    T=4d0
    P=0d0
    R1=R1*rho(P)
    R2=R2*rho(P)
    dseta=dseta*H_p(P)
    d=d*H_p(P)
    write(*,*) Pstep
    mu = (1d0/(R1)**2d0)+(1d0/(R2)**2d0)  !masa reducida
    Eini=0d0                         !energías inicial y final
    Efin=11d0
    ndata=int(Efin-Eini)*1d4
    write(*,*) ndata          !total de energías
    Estep=(Efin-Eini)/real(ndata)  !paso
    Ue1d = -2./sqrt(R1**2+dseta**2)!potenciales con la impureza
    Ue2d = -2./sqrt(R2**2+(d-dseta)**2)
    tol=1d-9
    vector_sale=0                  !inicio el vector que sale
    yini=0.1d0
    Fcero=yini
    En=Eini
    k=0
    j=0
    bandera_imprimir=0
    write(*,*) En,Efin,"energía inicial y final"
    WRITE(*,*) m(0d0,0d0,4d0)
    
    open(unit=1,file ="Funcion.txt")              !aquí guardo la función para graficar
    Open(unit=2, file="density.txt")   !aquí guardo la densidad de probabilidad
    open(unit=3, file="potencial.txt")  !aquí guardo el potencial una vez obtenida la función
    open(unit=4, file="energias.txt")


    Do while(En.le.Efin) !Ciclo que varia las energías
        bandera_imprimir=0
        call Rungekutta_anillo
        F=(Fcero-Fdospi)
        !                   write(1,*) F
        !                   write(*,*) F
        if((F*Fant.lt.0d0)) then
            j=j+1
            bandera_imprimir=0
            call biseccion(En_ant,En,Fant,F,En)
            call Rungekutta_anillo
            if(abs(Fcero-Fdospi) .lt. tol) then

                k=k+1
                !                write(*,*) "Er= ", ((GAM**2)/4d0)*(R1**2+R2**2) -real(k)*GAM+ (real(k)**2)/(R1**2+R2**2) + En


                if(k.eq.2) then
                    bandera_imprimir =1
                end if
                write(*,*) "Er= ", En
                call Rungekutta_anillo

                write(4,*) En

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


        end if

    End do

    Close(1)
    Close(2)
    close(3)
    close(4)

    open(unit=5, file="f_fo x 0 vs P.txt")
    open(unit=6, file="f_fo x 0_2 vs P.txt")
    open(unit=7, file="f_fo x 0_4 vs P.txt")
    open(unit=8, file="g_go vs P.txt")
    open(unit=9,file="rho vs P.txt")
    open(unit=10,file="H vs P.txt")
    Pini=0d0
    P=Pini
    Pfin=30d0
    Pstep=0.4d0

    do while(P.le.Pfin)
        P=P+Pstep
        write(5,*) f_(0d0,P,4d0),";",f_(0d0,P,200d0),";",f_(0d0,P,400d0),";",P
        write(6,*) f_(0.2d0,P,4d0),";",f_(0.2d0,P,200d0),";",f_(0.2d0,P,400d0),";",P
        write(7,*) f_(0.4d0,P,4d0),";",f_(0.4d0,P,200d0),";",f_(0.4d0,P,400d0),";",P
        write(8,*) g_(P,4d0),";",g_(P,200d0),";",g_(P,400d0),";",P
    end do
    P=Pini
    do while(P .le. Pfin)
        P=P+Pstep
        write(9,*) rho(P),";",P
        write(10,*)  H_p(P),";",P
    end do
    close(5)
    close(6)
    close(7)
    close(8)
    close(9)
    close(10)

end program
