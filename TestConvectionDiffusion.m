%% TEST Convection and diffusion with constant coefficient and fixed velocities

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
xR = 2;                 %surface
IMAX = 20;             %number of control volumes
dx = (xR-xL)/IMAX;      %mesh spacing
x = linspace(xL+dx/2,xR-dx/2,IMAX);
xv = linspace(xL,xR,IMAX+1);
tend = 20000000;             %set the final simulation time
time = 0;               %initial time
%set the initial condition (the right temperature-than it start freezing from the left)
for i=1:IMAX
    %psi(i) = -x(i);      % hydrostatic pressure
    %psi(i) = 0;
    T(i) = 10+x(i);
end
NMAX=100000;
for i=1:IMAX
    theta(i)=Thetaf(psi(i));
end
for nit=1:NMAX
    nit
    psiR=0.2;
    TR=10+xR;
    %left
    psiL = 0;
    TL=10;
    
    
    if(time>=tend)
        break
    end
    
    for i=1:IMAX+1
        if(i<IMAX/2)
            v(i)=-0.00001;
        else
            v(i)=-0.00001;
        end
    end
    
    vmax=max(v);
    %dt = 0.9/( vmax/dx + 2*(lambda)*(1/dx^2) );
    %if vmax==0
    dt=5;
    %     else
    %         dt = 0.9/( rho*cW/CT(thetas)*abs(vmax)/dx )
    %     end
    dt = 1;   %BTCS is unconditionalli stable (choose what you want)
    if(time+dt>tend)
        dt=tend-time;
    end
    
    
    plot(x,T,'o')
    title(sprintf('Current time = %f',time))
    drawnow
    %compute right hand side and the linear part of the system
    Tnew = ConvectionDiffusionProva(T,v);
    %T = Diffusion(T,theta);
    check=(Tnew-T')/dt
    time = time+dt;     %advance time
end










