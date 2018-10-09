%% DUMBSER'S LESSON

% BTCS scheme for the 1D Richards equation using the nested Newton method of
% Casulli & Zanolli.

clear all
close all
clc
global alpha thetas thetar n m Ks psic K KL KR di IMAX dx dt
%epsilon= 0.1; % regularization parameter: size of the regularization interval
%Phisical model parameters in SI units
day     = 24*3600;
Ks      = 0.062/day;    %[meter/second]
thetas  = 0.41;         %[-] saturated water content
thetar  = 0.095;        %[-] residuel water content
n       = 1.31;         %[-] parameter n
m       = 1-1/n;        %[-] parameter m
alpha   = 1.9;          %[m^(-1)]
psic    = -1/alpha*((n-1)/n)^(1/n);  %critical value of psi where the maximum of the derivative is located
%Domain
xL = 0;                 %bottom
xR = 2;                 %surface
IMAX = 320;             %number of control volumes
dx = (xR-xL)/IMAX;      %mesh spacing
x = linspace(xL+dx/2,xR-dx/2,IMAX);
tend = 3e5;             %set the final simulation time
time = 0;               %initial time
%set the initial condition (the right temperature-than it start freezing from the left)
for i=1:IMAX
    psi(i) = -x(i);      % hydrostatic pressure
    %psi(i) = 0;
end
NMAX=1000;
for nit=1:NMAX
    dt = 100;   %BTCS is unconditionalli stable (choose what you want)
    if(time+dt>tend)
        dt=tend-time;
    end
    %right boundary condition for pressure
    %if(time<=1e5)
    %    psiR = -0.05+0.03*sin(2*pi*time/1e5);
    %elseif(time>1e5 && time<=1.8e5)
    %    psiR = +0.1;
    %else
    %    psiR = -0.05+2952.45*exp(-time/18204.8);
    %end
    psiR=-2;
    %left
    psiL = 0;
    KR = kappa(psiR);
    KL = kappa(psiL);
    %Compute lambda (heat conduction coefficient) and the internal energy at the old time level
    for i=1:IMAX
        theta(i)=Thetaf(psi(i));
        K(i) = kappa(psi(i));
    end
    subplot(2,1,1)
    plot(x,psi,'o')
    subplot(2,1,2)
    plot(x,theta,'o')
    title(sprintf('Current time = %f',time))
    drawnow
    if(time>=tend)
        break
    end
    %compute right hand side and the linear part of the system
    
    for i=1:IMAX
        if(i==1)
            Kp = 0.5*(K(i)+K(i+1));
            Km = 0.5*(K(i)+KL);
            %dirichlet (pressure) boundary condition
            rhs(i) = theta(i) + dt/dx*(Kp-Km) + 2*Km*dt/dx^2*psiL;
        elseif(i==IMAX)
            Kp = 0.5*(K(i)+KR);
            Km = 0.5*(K(i)+K(i-1));
            %dirichlet (pressure boundary condition
            rhs(i) = theta(i) + dt/dx*(Kp-Km) + 2*Kp*dt/dx^2*psiR;
        else
            Kp = 0.5*(K(i)+K(i+1));
            Km = 0.5*(K(i)+K(i-1));
            rhs(i) = theta(i) + dt/dx*(Kp-Km);
        end
    end
    %now we start with the nested Newton method
    tol = 5e-15;
    % initial guess
    psi = min(psi,psic);
    disp(sprintf('Inizio loop'));
    for iNewton=1:100 %outer Newton iterations
        %The task of the outer iterations is ONLY to linearize one of the
        %two nonlinear functions q1 or q2  (è il ciclo in k sul quaderno)
        di = zeros(1,IMAX); %set the derivative of the nonlinear function to 0
        Mpsi = matop(psi);
        for i=1:IMAX
                f(i)= Theta1(psi(i))+ Mpsi(i)-rhs(i);
        end
        outres = sqrt(sum(f.*f)); %outer residual
%         disp(sprintf('Outer iteration %d, outres=%e',iNewton, outres));
        if(outres<tol)
            break %tolerance has been reached
        end
        psik = psi;             % save the value at the current outer iteration
        psi = max(psi,psic); % initial guess for the inner iteration
        for inner=1:100
            di = zeros(1,IMAX);
            Mpsi = matop(psi);
            for i=1:IMAX
                fk(i)=Theta1(psi(i))-(Theta2(psik(i))+dTheta2(psik(i))*(psi(i)-psik(i)))+ Mpsi(i)-rhs(i);  %Tk is frozen
                di(i)=dTheta1(psi(i))-dTheta2(psik(i));
            end
            inres=sqrt(sum(fk.*fk));
%             disp(sprintf(' -Inner iteration %d, inres= %e', inner,inres));
            if(inres<tol)
                break
            end
            dpsi = CGop(fk); %inner Newton step
            psi = psi -dpsi; %update psi at the inner iteration
        end
        
    end
    for i=1:IMAX
         if i==1
            k=0.5*(KL+K(i));
            vL = -k*(psi(i)-psiL)/(dx/2)-k;
         elseif i==IMAX
             k=0.5*(KR+K(i-1));
            vR = -k*(psiR-psi(i-1))/(dx/2)-k;
         end
         thetaNew(i) = Thetaf(psi(i));
    end
    error = (sum(thetaNew)-sum(theta)+dt/dx*(vR-vL))/sum(thetaNew)
             
    time = time+dt;     %advance time
end










