module subrutinas_funciones_parametros
    
    implicit none
    real(8),parameter ::PI=4*atan(1.), GAM = 0,h=0.06,m_o=9.1d-31!0.06 radianes !parametros del problema
    real(8),parameter :: S11=1.16d-3,S12=-3.7d-4
    real(8):: d=2d0, R1=1, R2=1,dseta=1d0,R,ds_ini=-1d0,ds_fin=3d0,ds_step=0.1d0
    real(8):: mu,Ecm,Er,En,En_ant,Uee,Ue1d,Ue2d,phi,tol,Fdospi,Fcero,F,Fant,yini!masa reducida, energías, potenciales,to
    real(8)::X=0.2d0,P=0,Pfin,Pini,Pstep,T=4 
    real(16)::cte_norm
    real(8),dimension(3):: vector_entra,vector_sale!vectores que utilizo para hacer RK4
    real(8),DIMENSION(420)::psi_value
    INTEGER :: i,j,k,bandera_imprimir,l
    INTEGER :: M_cm !!número cuántico de centro de masa,
    
    contains
    real(8) function rho(P)
    real(8)::rho_o
    real(8),intent(in)::P
    rho=sqrt((1d0-2d0*(S11+2d0*S12)*P))
    return
end function
real(8) function H_p(P)
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

real(8) function m(X,P,T)
implicit none
real(8), intent(in) :: X,P,T
m=m_o/(1d0+(PI_2(X)/3d0)*( (2d0/E_g_gm(X,P,T))+( 1d0 / (E_g_gm(X,P,T)+delta_o(X)) ) ) +delta(X) )
RETURN
end function
real(8) function PI_2(X)
implicit none
real(8),intent(in) ::X
PI_2=28900d0-6290d0*X
return
end function
real(8) function delta_o(X)
implicit none
real(8),intent(in) ::X
delta_o=341d0-66d0*X
return
end function

real(8) function delta(X)
implicit none
real(8),intent(in) ::X
delta=-3.935d0+0.488d0*X+4.938d0*(X**2)
return
end function
real(8) function E_g_gm(X,P,T)
implicit none
real(8),intent(in) ::X,P,T
E_g_gm=1519.4d0 + 1360d0*X+220d0*(X**2)+10.7d0*P - (0.5405d0*(T**2)/(204d0+T))
return
end function
real(8) function epsil(P,T)
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
!***************Biseccion***************!
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
    !vector_entra(1)=-2*PI   !condiciones iniciales
    vector_entra(1)=0.000001
    vector_entra(2)=yini
    vector_entra(3)=0.0d0
    i=0
    if(bandera_imprimir.eq.1)then
        !!!normalización
        DO WHILE(vector_entra(1).LE.(2d0*PI)) !ciclo para hacer RK4, para cada energía de la iteración correspondiente
            CALL RUNGE_KUTTA(vector_entra,vector_sale)
            vector_entra=vector_sale
            !valor de la función evaluada en 2pi
            Fdospi=vector_sale(2)
            i=i+1
            psi_value(i)=vector_sale(2)
            
            call integration_simpson(cte_norm)
        END DO
        WRITE(*,*) cte_norm
    end if
    
    !!impersión
    !vector_entra(1)=-2*PI   !condiciones iniciales
    vector_entra(1)=0.000001 
    vector_entra(2)=yini
    vector_entra(3)=0.01d0
    DO WHILE(vector_entra(1).LE.(2d0*PI)) !ciclo para hacer RK4, para cada energía de la iteración correspondiente
        
        CALL RUNGE_KUTTA(vector_entra,vector_sale)
        vector_entra=vector_sale
        !valor de la función evaluada en 2pi
        Fdospi=vector_sale(2)
        
        if(bandera_imprimir.eq.1) then  
            WRITE(1,*) vector_entra(2)/sqrt(cte_norm), vector_entra(1)
        else if(bandera_imprimir.eq.2)then
            WRITE(2,*) vector_entra(2)/sqrt(cte_norm), vector_entra(1)
        else if(bandera_imprimir.eq.3)then
            WRITE(7,*) vector_entra(2)/sqrt(cte_norm), vector_entra(1)
        end if
        
    END DO
    bandera_imprimir=0
end subroutine
!*************Integración Simpson**************
subroutine integration_simpson(valor)!!!!ya integra la (función)^2
    implicit none
    real(kind=16)::h,a,b,suma_1,suma_2
    real(kind=16),intent(out)::valor
    integer(4)::i,N
    suma_1=0
    suma_2=0
    !a= -2*PI
    a=0
    b= 2*PI
    N= 104!!! REMEMBER THAT N HAS TO BE EVEN
    h = abs(b-a)/N
    
    do i=1,((N/2)-1)
        suma_1 = suma_1 + (psi_value(2*i))**2
    end do
    do i=1,N/2
        suma_2 = suma_2 + (psi_value((2*i-1)))**2
    end do
    valor= (h/3.)*( psi_value(1)**2 + psi_value(105)**2 +2*suma_1 +4*suma_2 )
    
end subroutine
end module


program D_minus_P_T_X
    
    use subrutinas_funciones_parametros
    
    implicit none
    real(8)::Estep,Efin,Eini
    integer(8)::ndata
    
    R1=R1*rho(P)
    R2=R2*rho(P)
    R=R1
    d=d*H_p(P)
    mu = (1d0/(R1)**2d0)+(1d0/(R2)**2d0)  !masa reducida
    Ue1d = -2./sqrt(R1**2+dseta**2)!potenciales con la impureza
    Ue2d = -2./sqrt(R2**2+(d-dseta)**2)
    tol=1e-15
    vector_sale=0                  !inicio el vector que sale
    yini=0.001d0
    Fcero=yini
    M_cm=0
    k=0
    j=0
    l=0
    bandera_imprimir=0

   
    Eini=-20d0                !energías inicial y final
    Efin=50d0
    write(*,*) Eini,Efin,"energía inicial y final"
                 !aquí guardo la función para graficar
    open(unit=1,file ="funcion_base.txt")
    open(unit=2,file="psi.txt")
    open(unit=4, file="energias.txt")
    open(unit=3, file="energias_2.txt")  
    OPEN(unit=5, file="energias_3.txt")
    OPEN(unit=7, file="segundo_exi.txt")
    dseta=ds_ini
    Do while(dseta.le.ds_fin)
        
        bandera_imprimir=0
        dseta=dseta+ds_step
       
        dseta=dseta*H_p(P)
        Ue1d = -2./sqrt(R1**2+dseta**2)!potenciales con la impureza
        Ue2d = -2./sqrt(R2**2+(d-dseta)**2)
        
        ndata=int(Efin-Eini)*1d4 !# de datos de energía
        Estep=(Efin-Eini)/real(ndata)  !paso
        
        En=Eini
        k=0
        j=0
        Do while(En.le.Efin) !Ciclo que varia las energías
            
            call Rungekutta_anillo
            F=(Fcero-Fdospi)
            if((F*Fant.lt.0d0)) then
                j=j+1
                bandera_imprimir=0
                call biseccion(En_ant,En,Fant,F,En)
                call Rungekutta_anillo
                if(abs(Fcero-Fdospi) .lt. tol) then
                       
                    k=k+1
                    Ecm=((real(M_cm)**2)/(f_(X,P,T)*(R1**2+R2**2)) &
                    +  real(M_cm)*GAM/f_(X,P,T)  &
                    +  ((GAM**2)/(f_(X,P,T)*4d0))*(R1**2+R2**2)&
                    -(Ue1d+Ue2d)/(g_(P,T)) + En)*R**2
                    WRITE(*,*) "E+Ecm", Ecm
                    if(k.eq.1) then
                        write(4,*) Ecm,dseta !escribo los valores E vs R
                        write(*,*) Ecm,dseta ,"k=1"
                        bandera_imprimir=1
                        call Rungekutta_anillo
                    end if
                    if(k.eq.2) then
                        write(3,*) Ecm,dseta!escribo los valores E vs R
                        write(*,*) Ecm,dseta,"k=2"
                        bandera_imprimir=2
                        call Rungekutta_anillo
                    end if              
                    if(k.eq.3) then
                        write(5,*) Ecm,dseta !escribo los valores E vs R
                        write(*,*) Ecm,dseta ,"k=3"
                        bandera_imprimir=3
                        call Rungekutta_anillo
                        exit
                    end if
                    
                    call Rungekutta_anillo
                    
                    Fant=F
                    En_ant=En
                    j=0
                    En=En+Estep
                else if(j.gt.5d4)then
                    k=k+1
                    write(*,*) "la biseccion no converge satisfactoriamente"
                    
                    Ecm=((real(M_cm)**2)/(f_(X,P,T)*(R1**2+R2**2)) &
                    +  real(M_cm)*GAM/f_(X,P,T)  &
                    +  ((GAM**2)/(f_(X,P,T)*4d0))*(R1**2+R2**2)&
                    -(Ue1d+Ue2d)/(g_(P,T)) + En)*R**2
                    write(*,*) "un aproximado del autovalor es: ", En+Ecm
                    if(k.eq.1) then
                        write(4,*) Ecm,dseta !escribo los valores E vs R
                        write(*,*) Ecm,dseta ,"k=1"
                        bandera_imprimir=1
                        call Rungekutta_anillo
                    end if
                    if(k.eq.2) then
                        write(3,*) Ecm,dseta !escribo los valores E vs R
                        write(*,*) Ecm,dseta,"k=2"
                        bandera_imprimir=2
                        call Rungekutta_anillo
                        
                    end if              
                    if(k.eq.3) then
                        write(5,*) Ecm,dseta !escribo los valores E vs R
                        write(*,*) Ecm,dseta ,"k=3"
                        bandera_imprimir=3
                        call Rungekutta_anillo
                        exit
                    end if
                    Fant=F
                    En=En+Estep
                    
                    j=0
                end if
                
            else
                Fant=F
                En_ant=En
                En=En+Estep
            end if
        End Do
        
    End do
    close(1)
    close(2)
    close(3)
    close(4)
    close(5)
    close(7)
end program
