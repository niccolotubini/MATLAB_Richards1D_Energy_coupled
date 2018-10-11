%% BTCS scheme for the 1D Richards equation using the nested Newton method of
% Casulli & Zanolli 2005, 2010, 2017.
% This code integrate the mass conservation equation (Richards' equation
% coupled with surface flow) and the energy conservation equation written
% in conservative form.
% Resulting system are solved with Thomas algorithm

clear all
close all
clc
global alpha thetas thetar n m Ks psic cW rho K KL KR TR TL di IMAX dx dt Rain

%Phisical model parameters in SI units
thetaMethod = 1;
day     = 24*3600;
Ks      = 3.3333e-03;    %[meter/second]
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
x = [x,xR];
xv = linspace(xL,xR,IMAX+1);
tend = 207000;          %set the final simulation time
time = 0;               %initial time

% initialize variables
theta = zeros(1,IMAX);
thetaNew = zeros(1,IMAX);
volume = zeros(1,IMAX+1);
volumeNew = zeros(1,IMAX+1);
psi = zeros(1,IMAX+1);
T = zeros(1,IMAX+1);
a = zeros(1,IMAX+1);
b = zeros(1,IMAX+1);
c = zeros(1,IMAX+1);
rhs = zeros(1,IMAX+1);
K = zeros(1,IMAX+1);
v = zeros(1,IMAX+1);
f  = zeros(1,IMAX+1);
fk = zeros(1,IMAX+1);

%set the initial condition
for i=1:IMAX+1
    if(i==IMAX+1)
        psi(i) = -2.0;%psi(IMAX)-dx/2;   % psi at surface: if >0 ponding, otherwise not
        T(i) = 19;
    else
        psi(i) = -x(i);      % hydrostatic pressure
        T(i) = 19;
    end
    
end

% time cycle
NMAX=100000;
for nit=1:NMAX
    disp(sprintf('time iteration:%d', nit ));
    dt = 300;
    % Right boundary conditions
%      if (nit<=30)
        Rain=0/1000/300;
%      elseif( nit>30 && nit<90)
%          Rain = 3.3333e-03*2;
%      else 
%          Rain = 0;
%      end
    TR=19;
    % Left boundary conditions
    psiL =0.0;
    TL=19;
 
    KL = kappa(psiL,TL);
    for i=1:IMAX+1
        if(i==IMAX+1)
            volume(i) = Hf(psi(i));
            K(i) = kappa(psi(i)); % hydaraulic conductivity at soil surface
        else
            volume(i) = Thetaf(psi(i))*dx;
            theta(i) = Thetaf(psi(i));
            K(i) = kappa(psi(i));
        end
    end
    
    if(time>=tend)
        break
    end
    
   
    if(time+dt>tend)
        dt=tend-time;
    end
    
    % plot at time level n
    subplot(3,1,1)
    plot(x,psi,'o')
    ylabel('Psi [m]')
    title(sprintf('Current time = %f',time))
    subplot(3,1,2)
    plot(x,volume,'o')
    ylabel('Theta [-]')
     %subplot(4,1,3)
     %plot(xv,v,'o')
    subplot(3,1,3)
    plot(x,T,'r o')
    ylabel('Temperature [C]')
    xlabel('z [m]')
    drawnow
    
    %compute right hand side and the linear part of the system
    for i=1:IMAX+1
        if(i==1)
            Kp = 0.5*(K(i)+K(i+1));
            Km = 0.5*(K(i)+KL);
            %Kp = max(K(i),K(i+1));
            %Km = max(K(i),KL);
            a(i) = 0;                   % lower diagonal -> moved to the rhs multiplied by (-2)*psiL
            b(i) = +(2*Km+Kp)*dt/dx*thetaMethod;  % diagonal(not rhs!)
            c(i) = -Kp*dt/dx*thetaMethod;         % upper diagonal
            rhs(i) = volume(i) + dt*(Kp-Km) + dt*Kp*(1-thetaMethod)*(psi(i+1)-psi(i))/dx...
                    - 2*Km*dt/dx*(1-thetaMethod)*psi(i) + 2*Km*dt/dx*psiL;  % right hand side <- the "b" of lecture notes
        elseif(i==IMAX+1)
            Kp = 0;
            Km = 0.5*(K(i)+K(i-1));
            %Kp = 0;
            %Km = max(K(i),K(i-1));
            %rainfall rate (Rain) boundary condition
            a(i) = - dt/dx*Km*2*thetaMethod;
            b(i) = + dt/dx*(2*Km+0)*thetaMethod;
            c(i) = 0;
            rhs(i) = volume(i) + dt*(Rain-Km) - 2*dt*Km*(1-thetaMethod)*(psi(i)-psi(i-1))/dx;
        % questo if serve perche` a differenza del codice java in cui ho le
        % distanze dei centroidi qui ho solo un dx costante
        elseif(i==IMAX)
            Kp = 0.5*(K(i)+K(i+1));
            Km = 0.5*(K(i)+K(i-1));
            a(i) = - dt/dx*Km*thetaMethod;
            b(i) = + dt/dx*(Km+2*Kp)*thetaMethod;
            c(i) = - 2*dt/dx*Kp*thetaMethod;
            rhs(i) = volume(i) + dt*(Kp-Km) + 2*dt*Kp*(1-thetaMethod)*(psi(i+1)-psi(i))/dx...
                    - dt*Km*(1-thetaMethod)*(psi(i)-psi(i-1))/dx;    
        else
            Kp = 0.5*(K(i)+K(i+1));
            Km = 0.5*(K(i)+K(i-1));
            a(i) = -Km*dt/dx*thetaMethod;         % lower diagonal
            b(i) = +(Km+Kp)*dt/dx*thetaMethod;    % diagonal(not rhs!)
            c(i) = -Kp*dt/dx*thetaMethod;         % upper diagonal
            rhs(i) = volume(i) + dt*(Kp-Km) + dt*Kp*(1-thetaMethod)*(psi(i+1)-psi(i))/dx...
                    - dt*Km*(1-thetaMethod)*(psi(i)-psi(i-1))/dx;  % right hand side <- the "b" of lecture notes
        end
    end
    
    
    %% Newton method
    tol = 5e-10; %now we start with the nested Newton method
    
    % initial guess
    psi(IMAX+1) = max(psi(IMAX+1),0.1);
    psi(1:IMAX) = min(psi(1:IMAX),psic);%psic1
    
    for iNewton=1:100 %outer Newton iterations
        
        di = zeros(1,IMAX+1);
       
        for i=1:IMAX+1
            if(i==IMAX+1)
                f(i) = Hf(psi(i)) + a(i)*psi(i-1) + b(i)*psi(i) - rhs(i);
            elseif(i==1)
                f(i) = Thetaf(psi(i))*dx + a(i)*psiL + b(i)*psi(i) + c(i)*psi(i+1) - rhs(i);
            else
                f(i) = Thetaf(psi(i))*dx + a(i)*psi(i-1) + b(i)*psi(i) + c(i)*psi(i+1) - rhs(i);
            end
        end

        outres = sqrt(sum(f.*f)); %outer residual
        %disp(sprintf('    Outer iteration %d, outres=%e',iNewton, outres));
        if(outres<tol)
            break %tolerance has been reached
        end
        
        psik = psi;  % save the value at the current outer iteration
        psi(IMAX+1) = max(psi(IMAX+1),0.1);
        psi(1:IMAX) = max(psi(1:IMAX),psic);%psic1
        for inner=1:100 % inner Newton step
            
            di = zeros(1,IMAX+1);
            
            for i=1:IMAX+1
                if(i==IMAX+1)
                    fk(i) = H1(psi(i)) + a(i)*psi(i-1) + b(i)*psi(i) - rhs(i);
                    di(i) = dH1(psi(i));
                elseif(i==1)
                    fk(i) = Theta1(psi(i))*dx - ( Theta2(psik(i))*dx + dTheta2(psik(i))*dx*(psi(i)-psik(i)) ) + a(i)*psiL + b(i)*psi(i) + c(i)*psi(i+1) - rhs(i);
                    di(i) = dTheta1(psi(i))*dx - dTheta2(psik(i))*dx;
                else
                    fk(i) = Theta1(psi(i))*dx - ( Theta2(psik(i))*dx + dTheta2(psik(i))*dx*(psi(i)-psik(i)) ) + a(i)*psi(i-1) + b(i)*psi(i) + c(i)*psi(i+1) - rhs(i);
                    di(i) = dTheta1(psi(i))*dx - dTheta2(psik(i))*dx;
                end
            end

            inres=sqrt(sum(fk.*fk));
            %disp(sprintf('        Inner iteration %d, inres= %e', inner,inres));
            if(inres<tol)
                break
            end

            dpsi = Thomas(a,b+di,c,fk);
            psi(:) = psi(:) - dpsi(:);  %update psi at the inner iteration
            
        end
        
        
    end
    
    
    % compute new theta and velocities to be used as convection coefficients and error
    % mass computation
    for i=1:IMAX+1
        if i==1
            k=0.5*(KL+K(i));
            v(i) = -k*(psi(i)-psiL)/(dx/2)-k;
            thetaNew(i)=Thetaf(psi(i));
            volumeNew(i) = thetaNew(i)*dx;
            % no flux boundary
            % v(i) = 0;
        elseif i==IMAX+1
            k=0.5*(K(i)+K(i-1));
            v(i) = -k*(psi(i)-psi(i-1))/(dx/2)-k;
            volumeNew(i) = Hf(psi(i));
            % no flux boundary
            % v(i) = 0;
        else
            k=0.5*(K(i)+K(i-1));
            v(i) = -k*(psi(i)-psi(i-1))/dx-k;
            thetaNew(i)=Thetaf(psi(i));
            volumeNew(i) = thetaNew(i)*dx;
        end
    end
    disp(sprintf('    Ponding:%e', volumeNew(IMAX+1) ));
    errorMass(nit)=(sum(volumeNew) - sum(volume) - dt*(Rain - 0.5*( KL+K(1) )*( (psi(1)-psiL)/(dx/2) +1 ) ) );
    disp(sprintf('    Error mass:%e', errorMass(nit) ) );
    
    [T, errorHeat(nit)]= ConvectionDiffusionImplicit(T,v,theta,thetaNew,volume,volumeNew);
    %[T, errorHeat(nit)]= Diffusion(T,theta);
    disp(sprintf('    Error heat:%e', errorHeat(nit) ));
    %disp(sprintf('    Tmax:%f', max(T) ));
    minT(nit) = min(T);
    maxT(nit) = max(T);
    if(max(T)>TR+1*10^(-7) || min(T)<TR+1*10^(-7) )
        disp(sprintf('    differenze:%e    %e', max(T)-19, 19-min(T)  ))
        %break
    end
    time = time+dt;     %advance time
end


