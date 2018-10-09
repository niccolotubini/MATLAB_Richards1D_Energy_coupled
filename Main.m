% BTCS scheme for the 1D Richards equation using the nested Newton method of
% Casulli & Zanolli.
% Ricahrds equation is coupled with the thermal problem

clear all
close all
clc
global alpha thetas thetar n m Ks psic cW rho K KL KR TR TL di IMAX dx dt
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
rho = 1000;
cW = 4188;

%Domain
xL = 0;                 %bottom
xR = 1;                 %surface
IMAX = 160;             %number of control volumes
dx = (xR-xL)/IMAX;      %mesh spacing
x = linspace(xL+dx/2,xR-dx/2,IMAX);
xv = linspace(xL,xR,IMAX+1);
tend = 207000;             %set the final simulation time
time = 0;               %initial time
%set the initial condition (the right temperature-than it start freezing from the left)
for i=1:IMAX
    psi(i) = -x(i);      % hydrostatic pressure
    %psi(i) = 0;
    T(i) = 10;
end
NMAX=100000;
for nit=1:100
    nit
    psiR=-0.5;
    TR=15;
    %left
    psiL =0.0;
    TL=10;
    KR = kappa(psiR);
    KL = kappa(psiL);
    for i=1:IMAX
        theta(i)=Thetaf(psi(i));
        K(i) = kappa(psi(i));
    end
    psi1(nit)=psi(1);
    theta1(nit)=theta(1);
    if(time>=tend)
        break
    end
    
%     for i=1:IMAX+1
%         if i==1
%             k=0.5*(KL+K(i));
%             v(i) = -k*(psi(i)-psiL)/(dx/2)-k;
%         elseif i==IMAX+1
%             k=0.5*(KR+K(i-1));
%             v(i) = -k*(psiR-psi(i-1))/(dx/2)-k;
%         else
%             k=0.5*(K(i)+K(i-1));
%             v(i) = -k*(psi(i)-psi(i-1))/dx-k;
%         end
% %                 if abs(v(i))<10^(-14)
% %                     v(i)=0;
% %                 end
%     end
%     
%     vmax=max(abs(v));
    %     vmin=min(v);
    %dt = 0.9/( vmax/dx + 2*(lambda)*(1/dx^2) );
    %     if vmax==0 && vmin==0
    %         dt=5;
    %     else
    %         %dt =  min( 300, 0.9/( rho*cW/CT(thetas)*(-vmin+vmax)/dx ) );
    %dt = min(100, 0.9*CT(thetas)/( 1/dx*0.5*rho*cW*(2*vmax) ));
    %     end
           dt = 100;   %BTCS is unconditionalli stable (choose what you want)
    % store all dt
    dtn(nit)=dt;
    if(time+dt>tend)
        dt=tend-time;
    end
    
    
    subplot(3,1,1)
    plot(x,psi,'o')
    subplot(3,1,2)
    plot(x,theta,'o')
%     subplot(4,1,3)
%     plot(xv,v,'o')
    subplot(3,1,3)
    plot(x,T,'o')
    title(sprintf('Current time = %f',time))
    drawnow
    %compute right hand side and the linear part of the system
    
    for i=1:IMAX
        if(i==1)
            Kp = 0.5*(K(i)+K(i+1));
            Km = 0.5*(K(i)+KL);
            %% dirichlet (pressure) boundary condition
            rhs(i) = theta(i) + dt/dx*(Kp-Km) + 2*Km*dt/dx^2*psiL;
            fm=-dt/dx*Km+2*Km*dt/dx^2*psiL;
            %% no flux boundary
%             rhs(i) = theta(i) + dt/dx*(Kp);
        elseif(i==IMAX)
            Kp = 0.5*(K(i)+KR);
            Km = 0.5*(K(i)+K(i-1));
            %% dirichlet (pressure boundary condition
            rhs(i) = theta(i) + dt/dx*(Kp-Km) + 2*Kp*dt/dx^2*psiR;
            fp=dt/dx*Kp+ 2*Kp*dt/dx^2*psiR;
            %% no flux boundary
%             rhs(i) = theta(i) + dt/dx*(-Km);
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
        disp(sprintf('Outer iteration %d, outres=%e',iNewton, outres));
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
            %disp(sprintf(' -Inner iteration %d, inres= %e', inner,inres));
            if(inres<tol)
                break
            end
            dpsi = CGop(fk); %inner Newton step
            psi = psi -dpsi; %update psi at the inner iteration
        end
        
    end
    %compute new theta for u^n+1
    for i=1:IMAX
        thetaNew(i)=Thetaf(psi(i));
    end

    % compute new velocities to be used as advection coefficients and error
    % mass computation
    for i=1:IMAX+1
        if i==1
            k=0.5*(KL+K(i));
            v(i) = -k*(psi(i)-psiL)/(dx/2)-k;
            %% no flux boundary
%             v(i) = 0;
        elseif i==IMAX+1
            k=0.5*(KR+K(i-1));
            v(i) = -k*(psiR-psi(i-1))/(dx/2)-k;
            %% no flux boundary
%             v(i) = 0;
        else
            k=0.5*(K(i)+K(i-1));
            v(i) = -k*(psi(i)-psi(i-1))/dx-k;
        end
%                 if abs(v(i))<10^(-14)
%                     v(i)=0;
%                 end
    end
    
    errorMass(nit)=(sum(thetaNew)-sum(theta)+dt/dx*(v(IMAX+1)-v(1)));
    errorMass(nit)
    
    [T, errorHeat(nit)]= ConvectionDiffusion(T,v,theta,thetaNew);
    errorHeat(nit)
    
    TIMAX(nit)=T(IMAX);
    TMAX(nit)=max(T);
    TMIN(nit)=min(T);
        
    time = time+dt;     %advance time
end








%                                 for i=1:IMAX
%                                     if(i==1)
%                                         Kp = 0.5*(K(i)+K(i+1));
%                                         Km = 0.5*(K(i)+KL);
%                                         dirichlet (pressure) boundary condition
%                                         checkBis(i) = theta(i) - dt/dx*( v(i+1) - v(i) );
%                                     elseif(i==IMAX)
%                                         Kp = 0.5*(K(i)+KR);
%                                         Km = 0.5*(K(i)+K(i-1));
%                                         dirichlet (pressure boundary condition
%                                         checkBis(i) = theta(i) - dt/dx*( v(i+1) - v(i) );
%                                     else
%                                         Kp = 0.5*(K(i)+K(i+1));
%                                         Km = 0.5*(K(i)+K(i-1));
%                                         checkBis(i) = theta(i) - dt/dx*( v(i+1) - v(i) );
%                                     end
%                                 end


%                                 for i=1:IMAX
%                                     if(i==1)
%                                         Kp = 0.5*(K(i)+K(i+1));
%                                         Km = 0.5*(K(i)+KL);
%                                         %dirichlet (pressure) boundary condition
%                                         check(i) = theta(i) - dt/dx*( -Kp*(psi(i+1)-psi(i))/dx -Kp + Km*(psi(i)-psiL)/(dx/2) + Km );
%                                     elseif(i==IMAX)
%                                         Kp = 0.5*(K(i)+KR);
%                                         Km = 0.5*(K(i)+K(i-1));
%                                         %dirichlet (pressure boundary condition
%                                         check(i) = theta(i) - dt/dx*( -Kp*(psiR-psi(i))/(dx/2) -Kp + Km*(psi(i)-psi(i-1))/dx + Km );
%                                     else
%                                         Kp = 0.5*(K(i)+K(i+1));
%                                         Km = 0.5*(K(i)+K(i-1));
%                                         check(i) = theta(i) - dt/dx*( -Kp*(psi(i+1)-psi(i))/dx -Kp + Km*(psi(i)-psi(i-1))/dx + Km );
%                                     end
%                                 end




