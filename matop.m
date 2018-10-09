%% matop computes the product between the coefficient matrix and the unknown vector for Richards equation

function Apsi = matop(psi)
global di dx  dt K IMAX  KL KR
% diagonal part including the derivatives of the nonlinear functions
% linear part
% la matrice T è quadrata con solo elementi sulla diagonale
Apsi = di.*psi;
for i=1:IMAX
        if(i==1)
            Kp = 0.5*(K(i)+K(i+1));
            Km = 0.5*(K(i)+KL);
            Apsi(i) =  Apsi(i) - dt/dx^2*Kp*psi(i+1) + dt/dx^2*(2*Km+Kp)*psi(i);
            % no flux boundary
%             Apsi(i) =  Apsi(i) - dt/dx^2*Kp*psi(i+1) + dt/dx^2*(Kp)*psi(i) ;
        elseif(i==IMAX)
            Kp = 0.5*(K(i) + KR);
            Km = 0.5*(K(i) + K(i-1));
            Apsi(i) =  Apsi(i) + dt/dx^2*(Km+2*Kp)*psi(i) - dt/dx^2*Km*psi(i-1);
            % no flux boundary
%             Apsi(i) =  Apsi(i) + dt/dx^2*(Km)*psi(i) - dt/dx^2*Km*psi(i-1);
        else
            Kp = 0.5*(K(i) + K(i+1));
            Km = 0.5*(K(i) + K(i-1));
            Apsi(i) =  Apsi(i) - dt/dx^2*Kp*psi(i+1) + dt/dx^2*(Km+Kp)*psi(i) - dt/dx^2*Km*psi(i-1);
        end
end
