function Apsi = matopCoupled(psi)
global di dx  dt K IMAX  KL KR 
% diagonal part including the derivatives of the nonlinear functions
 %questo è quello che nel metodo di thomas vado a sommare a b
% linear part
% la matrice T è quadrata con solo elementi sulla diagonale
Apsi = di.*psi;
for i=1:IMAX+1
        if(i==1)
            Kp = 0.5*(K(i)+K(i+1));
            Km = 0.5*(K(i)+KL);
            Apsi(i) =  Apsi(i)+ dt/dx*(2*Km+Kp)*psi(i) - dt/dx*Kp*psi(i+1);
        elseif(i==IMAX+1)
            %Kp = 0.5*(K(i) + KR);
            Km = 0.5*(K(i) + K(i-1));
            Apsi(i) =  Apsi(i) + dt/dx*(Km)*psi(i) -dt/dx*Km*psi(i-1);
        else
            Kp = 0.5*(K(i) + K(i+1));
            Km = 0.5*(K(i) + K(i-1));
            Apsi(i) =  Apsi(i) -dt/dx*Kp*psi(i+1) + dt/dx*(Km+Kp)*psi(i) - dt/dx*Km*psi(i-1);
        end
end
