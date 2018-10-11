%% Convection and diffusion with coefficients depending on water content
% convection explicit in time, Godunov flux (upwind)
% diffusion implicit in time (time restriction due to diffusion is too severe)

function [Tnew, error]=ConvectionDiffusion(T,u,theta,thetaNew,volume,volumeNew)
global IMAX dx dt rho cW TR TL Rain
% Total internal energy at time level n
%E=0;
for i=1:IMAX+1
    if(i==1)       
        lambdap = 0.5*(Lambda(theta(i)+Lambda(theta(i+1))));
        lambdam = 0.5*(Lambda(theta(i)+Lambda(theta(i))));
        % convection flux fp at interface i+1/2 and fm at i-1/2
        fp =   rho*cW*( 0.5*u(i+1)*(T(i+1)+T(i)) - 0.5*abs(u(i+1))*(T(i+1)-T(i)) );
        fm =   rho*cW*( 0.5*u(i  )*(T(i  )+TL  ) - 0.5*abs(u(i  ))*(T(i  )-TL  ) );

        rhs(i) = CT(theta(i))*T(i)*dx - dt*( fp - fm ) + 2*dt/dx*lambdam*TL;
        a(i) = 0;
        b(i) = CT(thetaNew(i))*dx + dt/dx*(lambdap+2*lambdam);
        c(i) = -dt/dx*lambdap;
        
        heatFluxBottom = fm + 2/dx*lambdam*TL;
        E(i) = CT(theta(i))*T(i)*dx;
    elseif(i==IMAX+1)
        %lambdap = 0.5*(Lambda(theta(i)+Lambda(theta(i))));
        lambdam = 0.5*(Lambda(theta(i-1)+Lambda(theta(i-1))));
        % convection flux fp at interface i+1/2 and fm at i-1/2
        fp =  rho*cW*( 0.5*(-Rain)*(TR  +T(i  )) - 0.5*abs(-Rain)*(TR  -T(i  )) );
        fm =  rho*cW*( 0.5*u(i  )*(T(i)+T(i-1)) - 0.5*abs(u(i  ))*(T(i)-T(i-1)) );
        
        rhs(i) =  volume(i)*rho*cW*T(i) - dt*( fp - fm ) + 0.0;%2*dt/dx^2*lambdap*TR;
        a(i) = -2*dt/dx*lambdam;
        b(i) =  volumeNew(i)*rho*cW + 2*dt/dx*(lambdam);
        c(i) = 0;
        
        heatFluxTop = fp - 0.0;%2/dx*lambdap*TR;
        E(i) = volume(i)*rho*cW*T(i);
    elseif(i==IMAX)
        lambdap = 0.5*(Lambda(theta(i)+Lambda(theta(i))));
        lambdam = 0.5*(Lambda(theta(i-1)+Lambda(theta(i))));
        % convection flux fp at interface i+1/2 and fm at i-1/2
        fp =  rho*cW*( 0.5*u(i+1)*(T(i+1)  +T(i  )) - 0.5*abs(u(i+1))*(T(i+1)  -T(i  )) );
        fm =  rho*cW*( 0.5*u(i  )*(T(i)+T(i-1)) - 0.5*abs(u(i  ))*(T(i)-T(i-1)) );

        rhs(i) = CT(theta(i))*T(i)*dx - dt*( fp - fm );
        a(i) = -dt/dx*lambdam;
        b(i) = CT(thetaNew(i))*dx + dt/dx*(2*lambdap+lambdam);
        c(i) = -2*dt/dx*lambdap;
        
        %heatFluxTop = fp - 2/dx*lambdap*TR;
        E(i) = CT(theta(i))*T(i)*dx;
    else
        lambdap = 0.5*(Lambda(theta(i)+Lambda(theta(i+1))));
        lambdam = 0.5*(Lambda(theta(i-1)+Lambda(theta(i))));
        % convection flux fp at interface i+1/2 and fm at i-1/2
        fp = rho*cW*( 0.5*u(i+1)*(T(i+1)+T(i  )) - 0.5*abs(u(i+1))*(T(i+1)-T(i  )) );
        fm = rho*cW*( 0.5*u(i  )*(T(i  )+T(i-1)) - 0.5*abs(u(i  ))*(T(i  )-T(i-1)) );

        rhs(i) = CT(theta(i))*T(i)*dx - dt*( fp - fm );
        a(i) = -dt/dx*lambdam;
        b(i) = CT(thetaNew(i))*dx + dt/dx*(lambdap+lambdam);
        c(i) = -dt/dx*lambdap;
        
        E(i) = CT(theta(i))*T(i)*dx;
    end
    
end
Tnew=Thomas(a,b,c,rhs);

% compute the total internal energy at time level n+1
%Enew=0;
for i=1:IMAX
    Enew(i) = CT(thetaNew(i))*Tnew(i)*dx;
end
Enew(IMAX+1) = volumeNew(IMAX+1)*rho*cW*Tnew(i);

heatFluxTop    = heatFluxTop    + 0.0;%2*(0.5*(Lambda(theta(IMAX)+Lambda(theta(IMAX)))))*Tnew(IMAX)/dx;
heatFluxBottom = heatFluxBottom - 2*(0.5*(Lambda(theta(1)+Lambda(theta(1)))))*Tnew(1)/dx;

% Absolute error of internal energy E^n+1 = E^n - dt/dx*(fluxTop - fluxBottom)
error = sum(Enew) - sum(E) + dt*(heatFluxTop-heatFluxBottom);

end