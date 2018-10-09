%% Pure diffusion with coefficients depending on water content
% diffusion implicit in time (time restriction due to diffusion is too severe)
% it is possible to perform iPicard iteration

function [Tnew, error]=DiffusionNonLinear(T)
global dt dx TR TL IMAX

for iPicard=1:1
    % Total internal energy at time level n
    E=0;
    
    for i=1:IMAX
        if(i==1)
            lambdap = 0.5*( LambdaT(T(i))+LambdaT(T(i+1)) );
            lambdam = 0.5*( LambdaT(TL)+LambdaT(T(i)) ) ;
            
            rhs(i) = T(i) + 2*dt/dx^2*TL*lambdam;
            a(i) = 0;
            b(i) = 1  +dt/dx^2*(lambdap+2*lambdam);
            c(i) = -dt/dx^2*lambdap;
            
            heatFluxBottom = + 2/dx*lambdam*TL;
            E = E + 1*T(i);
        elseif(i==IMAX)
            lambdap = 0.5*( LambdaT(TR)+LambdaT(T(i)) );
            lambdam = 0.5*( LambdaT(T(i-1))+LambdaT(T(i)) );
            
            rhs(i) = T(i) + 2*dt/dx^2*TR*lambdap;
            a(i) = -dt/dx^2*lambdam;
            b(i) = 1 + dt/dx^2*(2*lambdap+lambdam);
            c(i) = 0;
            
            heatFluxTop = - 2/dx*lambdap*TR;
            E = E + 1*T(i);
        else
            lambdap = 0.5*( LambdaT(T(i))+LambdaT(T(i+1)) );
            lambdam = 0.5*( LambdaT(T(i-1))+LambdaT(T(i)) );
            
            rhs(i) = T(i);
            a(i) = -dt/dx^2*lambdam;
            b(i) = 1 + dt/dx^2*(lambdap+lambdam);
            c(i) = -dt/dx^2*lambdap;
            
            E = E + 1*T(i);
        end
    end
    Tnew = Thomas(a,b,c,rhs);
    T = Tnew;
end

% compute the total internal energy at time level n+1
Enew=0;
for i=1:IMAX
    Enew = Enew + 1*Tnew(i);
end

heatFluxTop    = heatFluxTop    + 2*0.5*( LambdaT(TR)+LambdaT(Tnew(IMAX)) )*Tnew(IMAX)/dx;
heatFluxBottom = heatFluxBottom - 2*0.5*( LambdaT(TL)+LambdaT(Tnew(1))    )*Tnew(1)/dx;

% Absolute error of internal energy E^n+1 = E^n - dt/dx*(fluxTop - fluxBottom)
error = Enew - E + dt/dx*(heatFluxTop-heatFluxBottom);

end