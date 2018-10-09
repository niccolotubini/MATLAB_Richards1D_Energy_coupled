% BTCS scheme for the 1D Richards equation using the nested Newton method of
% Casulli & Zanolli.
% Ricahrds equation is coupled with the thermal problem
% THOMAS Algorithm

clear all
close all
clc
global alpha thetas thetar n m Ks psic cW rho K KL KR TR TL di IMAX dx dt

%Phisical model parameters in SI units
day     = 24*3600;
Ks      = 0.062/day;    %[meter/second]
thetas  = 0.41;         %[-] saturated water content
thetar  = 0.095;        %[-] residuel water content
n       = 1.31;         %[-] parameter n
m       = 1-1/n;        %[-] parameter m
alpha   = 1.9;          %[m^(-1)]
psic    = -1/alpha*((n-1)/n)^(1/n);  %critical value of psi where the maximum of the derivative is located
rho = 1000;
cW = 4188;

%Domain
xL = 0;                 %bottom
xR = 2;                 %surface
IMAX = 320;             %number of control volumes
dx = (xR-xL)/IMAX;      %mesh spacing
x = linspace(xL+dx/2,xR-dx/2,IMAX);
xv = linspace(xL,xR,IMAX+1);
tend = 207000;          %set the final simulation time
time = 0;               %initial time

% initialize variables
theta = zeros(1,IMAX);
thetaNew = zeros(1,IMAX);
psi = zeros(1,IMAX);
T = zeros(1,IMAX);
K = zeros(1,IMAX);
v = zeros(1,IMAX+1);
rhs = zeros(1,IMAX);
a = zeros(1,IMAX);
b = zeros(1,IMAX);
c = zeros(1,IMAX);
f = zeros(1,IMAX);
fk = zeros(1,IMAX);

%set the initial condition
for i=1:IMAX
    psi(i) = -x(i);      % hydrostatic pressure
    T(i) = 19;
end

% time cycle
NMAX=100000;
for nit=1:NMAX
    disp(sprintf('time iteration:%d', nit ));
    % Right boundary conditions
    psiR=1.0;
    TR=220;
    % Left boundary conditions
    psiL =0.0;
    TL=19;
    
    KR = kappa(psiR,TR);
    KL = kappa(psiL,TL);
    for i=1:IMAX
        theta(i)=Thetaf(psi(i));
        K(i) = kappa(psi(i),T(i));
    end
    
    if(time>=tend)
        break
    end
    
    % compute velocities at time level n to determine the maximum time
    % step (this is an approximation)
    for i=1:IMAX+1
        if i==1
            k=0.5*(KL+K(i));
            v(i) = -k*(psi(i)-psiL)/(dx/2)-k;
        elseif i==IMAX+1
            k=0.5*(KR+K(i-1));
            v(i) = -k*(psiR-psi(i-1))/(dx/2)-k;
        else
            k=0.5*(K(i)+K(i-1));
            v(i) = -k*(psi(i)-psi(i-1))/dx-k;
        end
    end
    %
    vmax=max(abs(v));
    if vmax==0
           dt=1000; % no advection therefore no stability requirement on dt
    else
        dt = min(1000, 0.9*CT(thetas)/( 1/dx*0.5*rho*cW*(2*vmax) ));
    end
    dt = 300;
    % store all dt
    %dtn(nit)=dt;
    if(time+dt>tend)
        dt=tend-time;
    end
    
    % plot at time level n
    subplot(3,1,1)
    plot(x,psi,'o')
    ylabel('Psi [m]')
    title(sprintf('Current time = %f',time))
    subplot(3,1,2)
    plot(x,theta,'o')
    ylabel('Theta [-]')
%     subplot(4,1,3)
%     plot(xv,v,'o')
    subplot(3,1,3)
    plot(x,T,'r o')
    ylabel('Temperature [C]')
    xlabel('z [m]')
    drawnow
    
    % compute right hand side and the linear part of the system
    for i=1:IMAX
        if(i==1)
            Kp = 0.5*(K(i)+K(i+1));
            Km = 0.5*(K(i)+KL);
            % dirichlet (pressure) boundary condition
            rhs(i) = theta(i) + dt/dx*(Kp-Km) + 2*Km*dt/dx^2*psiL;
            a(i) = 0;
            b(i) = dt/dx^2*(2*Km+Kp);
            c(i) = -Kp*dt/dx^2;
            
            % no flux boundary
            % rhs(i) = theta(i) + dt/dx*(Kp);
        elseif(i==IMAX)
            Kp = 0.5*(K(i)+KR);
            Km = 0.5*(K(i)+K(i-1));
            % dirichlet (pressure boundary condition
            rhs(i) = theta(i) + dt/dx*(Kp-Km) + 2*Kp*dt/dx^2*psiR;
            a(i) = -Km*dt/dx^2;
            b(i) = dt/dx^2*(Km+2*Kp);
            c(i) = 0;
           
            % no flux boundary
            % rhs(i) = theta(i) + dt/dx*(-Km);
        else
            Kp = 0.5*(K(i)+K(i+1));
            Km = 0.5*(K(i)+K(i-1));
            rhs(i) = theta(i) + dt/dx*(Kp-Km);
            a(i) = -Km*dt/dx^2;
            b(i) = dt/dx^2*(Km+Kp);
            c(i) = -Kp*dt/dx^2;
        end
    end
    
    tol = 5e-15; %now we start with the nested Newton method
    
    psi = min(psi,psic); % initial guess
    for iNewton=1:100 %outer Newton iterations
        
        for i=1:IMAX
            %compute the true non linear function f(psi)
            f(i) = Thetaf(psi(i))-rhs(i);
            if(i==1)
                f(i) = f(i) + b(i)*psi(i) + c(i)*psi(i+1);
            elseif(i==IMAX)
                f(i) = f(i)+a(i)*psi(i-1)+b(i)*psi(i);
            else
                f(i) = f(i)+a(i)*psi(i-1)+b(i)*psi(i)+c(i)*psi(i+1);
            end
        end
        outres = sqrt(sum(f.*f)); %outer residual
        %disp(sprintf('Outer iteration %d, outres=%e',iNewton, outres));
        if(outres<tol)
            break %tolerance has been reached
        end
        
        psik = psi;             % save the value at the current outer iteration
        psi = max(psi,psic); % initial guess for the inner iteration
        for inner=1:100
            
            for i=1:IMAX
                fk(i)=Theta1(psi(i)) - ( Theta2(psik(i)) + dTheta2(psik(i))*(psi(i)-psik(i)) ) - rhs(i);  %Tk is frozen
                if(i==1)
                    fk(i)=fk(i)+b(i)*psi(i)+c(i)*psi(i+1);
                elseif(i==IMAX)
                    fk(i)=fk(i)+a(i)*psi(i-1)+b(i)*psi(i);
                else
                    fk(i)=fk(i)+a(i)*psi(i-1)+b(i)*psi(i)+c(i)*psi(i+1);
                end
                di(i)=dTheta1(psi(i))-dTheta2(psik(i));
            end
            inres=sqrt(sum(fk.*fk));
            %disp(sprintf(' -Inner iteration %d, inres= %e', inner,inres));
            if(inres<tol)
                break
            end
            dpsi = Thomas(a,b+di,c,fk); %inner Newton step
            psi = psi(:)-dpsi(:); %update temperature at the inner iteration
        end
        
    end
    
    
    % compute new theta and velocities to be used as convection coefficients and error
    % mass computation
    for i=1:IMAX+1
        if i==1
            k=0.5*(KL+K(i));
            v(i) = -k*(psi(i)-psiL)/(dx/2)-k;
            thetaNew(i)=Thetaf(psi(i));
            % no flux boundary
            % v(i) = 0;
        elseif i==IMAX+1
            k=0.5*(KR+K(i-1));
            v(i) = -k*(psiR-psi(i-1))/(dx/2)-k;
            % no flux boundary
            % v(i) = 0;
        else
            k=0.5*(K(i)+K(i-1));
            v(i) = -k*(psi(i)-psi(i-1))/dx-k;
            thetaNew(i)=Thetaf(psi(i));
        end
    end
    
    errorMass(nit)=(sum(thetaNew) - sum(theta) + dt/dx*(v(IMAX+1)-v(1)));
    disp(sprintf('    Error mass:%e', errorMass(nit) ));
    
    [T, errorHeat(nit)]= ConvectionDiffusionImplicit(T,v,theta,thetaNew);
    %[T, errorHeat(nit)]= Diffusion(T,theta);
    disp(sprintf('    Error heat:%e', errorHeat(nit) ));
    disp(sprintf('    Tmax:%f', max(T) ));
    if(max(T)>TR)
        break
    end
    time = time+dt;     %advance time
end


