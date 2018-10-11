%% Convection and diffusion with coefficients depending on water content
% convection explicit in time, Godunov flux (upwind)
% diffusion implicit in time (time restriction due to diffusion is too severe)

function [Tnew, error]=ConvectionDiffusionImplicit(T,u,theta,thetaNew,volume,volumeNew)
global IMAX dx dt rho cW TR TL Rain
% Total internal energy at time level n
E=0;
for i=1:IMAX+1
    if(i==1)       
        lambdap = 0.5*(Lambda(theta(i))+Lambda(theta(i+1)));
        lambdam = 0.5*(Lambda(theta(i))+Lambda(theta(i)));
        % convection flux fp at interface i+1/2 and fm at i-1/2
        %fp =   rho*cW*( 0.5*u(i+1)*(T(i+1)+T(i)) - 0.5*abs(u(i+1))*(T(i+1)-T(i)) );
        %fm =   rho*cW*( 0.5*u(i  )*(T(i  )+TL  ) - 0.5*abs(u(i  ))*(T(i  )-TL  ) );

        rhs(i) = CT(theta(i))*T(i)*dx - dt*rho*cW*( -0.5*u(i  )*TL - 0.5*abs(u(i  ))*TL ) + 2*dt/dx*lambdam*TL;
        a(i) = 0;
        b(i) = CT(thetaNew(i))*dx + dt/dx*(lambdap+2*lambdam) + dt*rho*cW*(  0.5*u(i+1) + 0.5*abs(u(i+1)) - ( 0.5*u(i) - 0.5*abs(u(i)) ) );
        c(i) = -dt/dx*lambdap + dt*rho*cW*( 0.5*u(i+1) - 0.5*abs(u(i+1)) );
        
        heatFluxBottom = 2/dx*lambdam*TL + rho*cW*(0.5*u(i  )*TL + 0.5*abs(u(i  ))*TL );
        E(i) = + CT(theta(i))*T(i)*dx;
    elseif(i==IMAX+1)
        %lambdap = 0.5*(Lambda(theta(i-1))+Lambda(theta(i-1)));
        lambdam = 0.5*(Lambda(theta(i-1))+Lambda(theta(i-1)));
        % convection flux fp at interface i+1/2 and fm at i-1/2
        %fp =  rho*cW*( 0.5*u(i+1)*(TR  +T(i  )) - 0.5*abs(u(i+1))*(TR  -T(i  )) );
        %fm =  rho*cW*( 0.5*u(i  )*(T(i)+T(i-1)) - 0.5*abs(u(i  ))*(T(i)-T(i-1)) );

        rhs(i) = volume(i)*rho*cW*T(i) + dt*0.00 + dt*rho*cW*( Rain*TR );
        a(i) = -2*dt/dx*lambdam + dt*rho*cW*(  -0.5*u(i  ) - 0.5*abs(u(i  )) );
        b(i) = volumeNew(i)*rho*cW + 2*dt/dx*(lambdam) + dt*rho*cW*( - 0.5*u(i) + 0.5*abs(u(i) ) );
        c(i) = 0;
        
        heatFluxTop = 0 + rho*cW*Rain*TR;
        E(i) = + volumeNew(i)*rho*cW*T(i);
    elseif(i==IMAX)
        lambdap = 0.5*(Lambda(theta(i))+Lambda(theta(i)));
        lambdam = 0.5*(Lambda(theta(i-1))+Lambda(theta(i)));
        % convection flux fp at interface i+1/2 and fm at i-1/2
        %fp =  rho*cW*( 0.5*u(i+1)*(TR  +T(i  )) - 0.5*abs(u(i+1))*(TR  -T(i  )) );
        %fm =  rho*cW*( 0.5*u(i  )*(T(i)+T(i-1)) - 0.5*abs(u(i  ))*(T(i)-T(i-1)) );

        rhs(i) = CT(theta(i))*T(i)*dx;
        a(i) = -dt/dx*lambdam + dt*rho*cW*(  -0.5*u(i  ) - 0.5*abs(u(i  )) );
        b(i) = CT(thetaNew(i))*dx + dt/dx*(2*lambdap+lambdam) + dt*rho*cW*(  0.5*u(i+1) + 0.5*abs(u(i+1)) - 0.5*u(i) + 0.5*abs(u(i) ) );
        c(i) = -2*dt/dx*lambdap + dt*rho*cW*( 0.5*u(i+1) - 0.5*abs(u(i+1)) );
        
        %heatFluxTop = - 2/dx*lambdap*TR + rho*cW*( 0.5*u(i+1)*TR  - 0.5*abs(u(i+1))*TR );
        E(i) = + CT(theta(i))*T(i)*dx;
    else
        lambdap = 0.5*(Lambda(theta(i))+Lambda(theta(i+1)));
        lambdam = 0.5*(Lambda(theta(i-1))+Lambda(theta(i)));
        % convection flux fp at interface i+1/2 and fm at i-1/2
%         fp = rho*cW*( 0.5*u(i+1)*(T(i+1)+T(i  )) - 0.5*abs(u(i+1))*(T(i+1)-T(i  )) );
%         fm = rho*cW*( 0.5*u(i  )*(T(i  )+T(i-1)) - 0.5*abs(u(i  ))*(T(i  )-T(i-1)) );

        rhs(i) = CT(theta(i))*T(i)*dx;
        a(i) = -dt/dx*lambdam + dt*rho*cW*( -0.5*u(i) - 0.5*abs(u(i)) );
        b(i) = CT(thetaNew(i))*dx + dt/dx*(lambdap+lambdam)  + dt*rho*cW*(  0.5*u(i+1) + 0.5*abs(u(i+1)) - 0.5*u(i) + 0.5*abs(u(i) ) );
        c(i) = -dt/dx*lambdap + dt*rho*cW*( 0.5*u(i+1) - 0.5*abs(u(i+1)) );
        
        E(i) = + CT(theta(i))*T(i)*dx;
    end
    
end
Tnew=Thomas(a,b,c,rhs);

% compute the total internal energy at time level n+1
Enew=0;
for i=1:IMAX
    Enew(i) = + CT(thetaNew(i))*Tnew(i)*dx;
end
Enew(IMAX+1) =+volumeNew(IMAX+1)*rho*cW*Tnew(IMAX+1);
heatFluxTop    = heatFluxTop;%    + 2*(0.5*( Lambda(theta(IMAX)) + Lambda(theta(IMAX)) ) ) *Tnew(IMAX)/dx;
heatFluxBottom = heatFluxBottom - 2*(0.5*(Lambda(theta(1))+Lambda(theta(1))))*Tnew(1)/dx;

% Absolute error of internal energy E^n+1 = E^n - dt/dx*(fluxTop - fluxBottom)
%error = Enew - E + dt/dx*(heatFluxTop-heatFluxBottom);
Enew-E;
error = sum(Enew) - sum(E) + dt*( rho*cW*( Rain*TR )  ...
                            -rho*cW*0.5*( u(1)*(Tnew(1)+TL) - abs(u(1))*(Tnew(1)-TL) ) + Lambda(theta(1))*2*(Tnew(1)-TL)/dx );
error = error;
end