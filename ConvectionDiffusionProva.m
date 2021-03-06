% convection and diffusion of a scalar field T
function [Tnew, error]=ConvectionDiffusionProva(T,u)
global IMAX dx dt TR TL
%Tnew = T;
E=0;
cT=100000000;
rho=1000;
cW=1;
for i=1:IMAX
    if(i==1)

        fp =   rho*cW*( 0.5*u(i+1)*(T(i+1)+T(i)) - 0.5*abs(u(i+1))*(T(i+1)-T(i)) );
        fm =   rho*cW*( 0.5*u(i  )*(T(i  )+TL  ) - 0.5*abs(u(i  ))*(T(i  )-TL  ) );
%         fp = 0;
%         fm = 0;
        lambdap = 1;
        lambdam = 1;
        
        rhs(i) = cT*T(i) - dt/dx*( fp - fm ) + 2*dt/dx^2*lambdam*TL;
        a(i) = 0;
        b(i) = cT + dt/dx^2*(lambdap+2*lambdam);
        c(i) = -dt/dx^2*lambdap;
        heatFluxBottom = fm + 2/dx*lambdam*TL;
        E = E + cT*T(i);
    elseif(i==IMAX)
        
        fp =  rho*cW*( 0.5*u(i+1)*(TR  +T(i  )) - 0.5*abs(u(i+1))*(TR  -T(i  )) );
        fm =  rho*cW*( 0.5*u(i  )*(T(i)+T(i-1)) - 0.5*abs(u(i  ))*(T(i)-T(i-1)) );
%         fp = 0;
%         fm = 0;
        lambdap = 1;
        lambdam = 1;
       
        rhs(i) = cT*T(i) - dt/dx*( fp - fm ) + 2*dt/dx^2*lambdap*TR;
        a(i) = -dt/dx^2*lambdam;
        b(i) = cT + dt/dx^2*(2*lambdap+lambdam);
        c(i) = 0;
        heatFluxTop = fp - 2/dx*lambdap*TR;
        E = E + cT*T(i);
    else

        fp = rho*cW*( 0.5*u(i+1)*(T(i+1)+T(i  )) - 0.5*abs(u(i+1))*(T(i+1)-T(i  )) );
        fm = rho*cW*( 0.5*u(i  )*(T(i  )+T(i-1)) - 0.5*abs(u(i  ))*(T(i  )-T(i-1)) );
%         fp = 0;
%         fm = 0;
        lambdap = 1;
        lambdam = 1;
        
        rhs(i) = cT*T(i) - dt/dx*( fp - fm );
        a(i) = -dt/dx^2*lambdam;
        b(i) = cT + dt/dx^2*(lambdap+lambdam);
        c(i) = -dt/dx^2*lambdap;
        E = E + cT*T(i);
    end
    
end
Tnew=Thomas(a,b,c,rhs);
Enew=0;
for i=1:IMAX
    Enew = Enew + cT*Tnew(i);
end
heatFluxTop    = heatFluxTop    + 2*(0.5*(1+1))*Tnew(IMAX)/dx;
heatFluxBottom = heatFluxBottom - 2*(0.5*(1+1))*Tnew(1)/dx;
Echeck = E -dt/dx*(heatFluxTop-heatFluxBottom);
error=(Enew-Echeck)/Enew
-(rho*cW)/cT*u(1)



