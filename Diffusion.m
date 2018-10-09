%% Pure diffusion with coefficients depending on water content
% diffusion implicit in time (time restriction due to diffusion is too severe)

%% CORREGGERE ERROR HEAT
%%
function [Tnew, error]=Diffusion(T,theta)
global IMAX dx dt TR TL
% Total internal energy at time level n
E=0;
for i=1:IMAX
    if(i==1)
        lambdap = 0.5*( Lambda(theta(i))+Lambda(theta(i+1)) );
        lambdam = 0.5*( Lambda(theta(i))+Lambda(theta(i  )) );
        
        rhs(i) = CT(theta(i))*T(i) + 2*dt/dx^2*TL*lambdam;
        a(i) = 0;
        b(i) = CT(theta(i))+dt/dx^2*(lambdap+2*lambdam);
        c(i) = -dt/dx^2*lambdap;
        
        heatFluxBottom =+ 2/dx*lambdam*TL;
        E = E + CT(theta(i))*T(i);
    elseif(i==IMAX)
        lambdap = 0.5*(Lambda(theta(i))+Lambda(theta(i)));
        lambdam = 0.5*(Lambda(theta(i-1))+Lambda(theta(i)));
        
        rhs(i) = CT(theta(i))*T(i) + 2*dt/dx^2*TR*lambdap;
        a(i) = -dt/dx^2*lambdam;
        b(i) = CT(theta(i))+dt/dx^2*(2*lambdap+lambdam);
        c(i) = 0;
        heatFluxTop = - 2/dx*lambdap*TR;
        E = E + CT(theta(i))*T(i);
    else
        lambdap = 0.5*(Lambda(theta(i))+Lambda(theta(i+1)));
        lambdam = 0.5*(Lambda(theta(i-1))+Lambda(theta(i)));
        
        rhs(i) =CT(theta(i))*T(i);
        a(i) = -dt/dx^2*lambdam;
        b(i) = CT(theta(i))+dt/dx^2*(lambdap+lambdam);
        c(i) = -dt/dx^2*lambdap;
        E = E + CT(theta(i))*T(i);
    end
end
Tnew = Thomas(a,b,c,rhs);

% compute the total internal energy at time level n+1
Enew=0;
for i=1:IMAX
    Enew = Enew + CT(theta(i))*Tnew(i);
end

heatFluxTop    = heatFluxTop    + 2*Lambda(theta(IMAX))*Tnew(IMAX)/dx;
heatFluxBottom = heatFluxBottom - 2*Lambda(theta(1))*Tnew(1)/dx;

% Absolute error of internal energy E^n+1 = E^n - dt/dx*(fluxTop - fluxBottom)
error = Enew - E + dt/dx*(heatFluxTop-heatFluxBottom);

end