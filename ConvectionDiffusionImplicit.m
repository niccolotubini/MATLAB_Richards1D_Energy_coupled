%% Convection and diffusion with coefficients depending on water content
% convection explicit in time, Godunov flux (upwind)
% diffusion implicit in time (time restriction due to diffusion is too severe)

function [Tnew, error]=ConvectionDiffusionImplicit(T,u,theta,thetaNew)
global IMAX dx dt rho cW TR TL
% Total internal energy at time level n
E=0;
for i=1:IMAX
    if(i==1)       
        lambdap = 0.5*(Lambda(theta(i))+Lambda(theta(i+1)));
        lambdam = 0.5*(Lambda(theta(i))+Lambda(theta(i)));
        % convection flux fp at interface i+1/2 and fm at i-1/2
        %fp =   rho*cW*( 0.5*u(i+1)*(T(i+1)+T(i)) - 0.5*abs(u(i+1))*(T(i+1)-T(i)) );
        %fm =   rho*cW*( 0.5*u(i  )*(T(i  )+TL  ) - 0.5*abs(u(i  ))*(T(i  )-TL  ) );

        rhs(i) = CT(theta(i))*T(i) - dt/dx*rho*cW*( -0.5*u(i  )*TL - 0.5*abs(u(i  ))*TL ) + 2*dt/dx^2*lambdam*TL;
        a(i) = 0;
        b(i) = CT(thetaNew(i)) + dt/dx^2*(lambdap+2*lambdam) + dt/dx*rho*cW*(  0.5*u(i+1) + 0.5*abs(u(i+1)) - ( 0.5*u(i) - 0.5*abs(u(i)) ) );
        c(i) = -dt/dx^2*lambdap + dt/dx*rho*cW*( 0.5*u(i+1) - 0.5*abs(u(i+1)) );
        
        heatFluxBottom = 2/dx*lambdam*TL + rho*cW*(0.5*u(i  )*TL + 0.5*abs(u(i  ))*TL );
        E = E + CT(theta(i))*T(i);
    elseif(i==IMAX)
        lambdap = 0.5*(Lambda(theta(i))+Lambda(theta(i)));
        lambdam = 0.5*(Lambda(theta(i-1))+Lambda(theta(i)));
        % convection flux fp at interface i+1/2 and fm at i-1/2
        %fp =  rho*cW*( 0.5*u(i+1)*(TR  +T(i  )) - 0.5*abs(u(i+1))*(TR  -T(i  )) );
        %fm =  rho*cW*( 0.5*u(i  )*(T(i)+T(i-1)) - 0.5*abs(u(i  ))*(T(i)-T(i-1)) );

        rhs(i) = CT(theta(i))*T(i) + 2*dt/dx^2*lambdap*TR - dt/dx*rho*cW*( 0.5*u(i+1)*TR  - 0.5*abs(u(i+1))*TR );
        a(i) = -dt/dx^2*lambdam + dt/dx*rho*cW*(  -0.5*u(i  ) - 0.5*abs(u(i  )) );
        b(i) = CT(thetaNew(i)) + dt/dx^2*(2*lambdap+lambdam) + dt/dx*rho*cW*(  0.5*u(i+1) + 0.5*abs(u(i+1)) - ( 0.5*u(i) - 0.5*abs(u(i)) ) );
        c(i) = 0;
        
        heatFluxTop = - 2/dx*lambdap*TR + rho*cW*( 0.5*u(i+1)*TR  - 0.5*abs(u(i+1))*TR );
        E = E + CT(theta(i))*T(i);
    else
        lambdap = 0.5*(Lambda(theta(i))+Lambda(theta(i+1)));
        lambdam = 0.5*(Lambda(theta(i-1))+Lambda(theta(i)));
        % convection flux fp at interface i+1/2 and fm at i-1/2
%         fp = rho*cW*( 0.5*u(i+1)*(T(i+1)+T(i  )) - 0.5*abs(u(i+1))*(T(i+1)-T(i  )) );
%         fm = rho*cW*( 0.5*u(i  )*(T(i  )+T(i-1)) - 0.5*abs(u(i  ))*(T(i  )-T(i-1)) );

        rhs(i) = CT(theta(i))*T(i);
        a(i) = -dt/dx^2*lambdam + dt/dx*rho*cW*( -0.5*u(i) - 0.5*abs(u(i)) );
        b(i) = CT(thetaNew(i)) + dt/dx^2*(lambdap+lambdam)  + dt/dx*rho*cW*(  0.5*u(i+1) + 0.5*abs(u(i+1)) - ( 0.5*u(i) - 0.5*abs(u(i)) ) );
        c(i) = -dt/dx^2*lambdap + dt/dx*rho*cW*( 0.5*u(i+1) - 0.5*abs(u(i+1)) );
        
        E = E + CT(theta(i))*T(i);
    end
    
end
Tnew=Thomas(a,b,c,rhs);

% compute the total internal energy at time level n+1
Enew=0;
for i=1:IMAX
    Enew = Enew + CT(thetaNew(i))*Tnew(i);
end
heatFluxTop    = heatFluxTop    + 2*(0.5*( Lambda(theta(IMAX)) + Lambda(theta(IMAX)) ) ) *Tnew(IMAX)/dx;
heatFluxBottom = heatFluxBottom - 2*(0.5*(Lambda(theta(1))+Lambda(theta(1))))*Tnew(1)/dx;

% Absolute error of internal energy E^n+1 = E^n - dt/dx*(fluxTop - fluxBottom)
%error = Enew - E + dt/dx*(heatFluxTop-heatFluxBottom);

error = Enew - E + dt/dx*( rho*cW*0.5*( u(IMAX+1)*(TR+Tnew(IMAX)) - abs(u(IMAX+1))*(TR-Tnew(IMAX)) ) -Lambda(theta(IMAX))*2*(TR-Tnew(IMAX))/dx ...
                            -rho*cW*0.5*( u(1)*(Tnew(1)+TL) - abs(u(1))*(Tnew(1)-TL) ) + Lambda(theta(1))*2*(Tnew(1)-TL)/dx );
error = error;
end