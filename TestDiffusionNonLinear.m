%% Test DiffusionNonLinear.m for non-linear conduction equation
% Page 105 Dall'Amico Ph.D. thesis
% analytical solution Dumbser et al. 2008
%       T(x,t) = min{ 1, max[0,0.5*(1- (x-x0)/sqrt(epsilon t) ) ] }
% in DiffusionNonLinear.m it is possible to test Picard iteration in order
% to see how the solution improves. Picard iteration works on the thermal
% coefficient: lambda is always computed at time level n (semi-implicit in time)

% clear all
% close all
% clc
global TR TL IMAX dx dt epsilon

%Parameters
delta = 10^(-6);
epsilon = 10^(-3);
cT = 1; % has to be constant

%Domain
xL = 0;         %bottom
xR = 20;        %surface
x0 = 10;        % discuntinuity position
IMAX = 1000;             %number of control volumes
dx = (xR-xL)/IMAX;      %mesh spacing
x = linspace(xL+dx/2,xR-dx/2,IMAX);
tend = 200;             %set the final simulation time
time = 0;               %initial time

% Initialize variables
T = zeros(1,IMAX);
Texact = zeros(1,IMAX);
lambda = zeros(1,IMAX);

%set the initial condition (the right temperature-than it start freezing from the left)
for i=1:IMAX
    if(x(i)<x0)
        T(i) = 1 - delta;
    else
        T(i) = delta;
    end
end
figure(1)
plot(x,T,'o')
title(sprintf('Initial condition'))



NMAX=200000;
for nit=1:NMAX
    % Boundary conditions
    TR = delta;
    TL = 1-delta;
    
    dt = 2; % the sheme BTCS is unconditionally stable, increasing dt the accuracy decreases
    if(time>=tend)
        break
    end 
    if(time+dt>tend)
        dt=tend-time;
    end

    for i=1:IMAX
        % analytical solution at t=time
        Texact(i) = min( 1,max( 0, 0.5*( 1-(x(i)-x0)/sqrt(epsilon*time) ) ) );
    end
    
    
    
    figure(2)
    plot(x,T,'o')
    hold on
    plot(x,Texact,'r')
    title(sprintf('Current time = %f',time))
    drawnow
    hold off
    
    %compute right hand side and the linear part of the system
    [T, error] = DiffusionNonLinear(T);

    disp(sprintf('   errorHeat:%e\n', error) )
    
    relativeError(nit) = sum(abs(T'-Texact))/sum(abs(Texact));
    norm(nit) = sqrt( sum( (T'-Texact).*(T'-Texact) ) ); % norma del vettore degli scarti T-Texact
    
    disp(sprintf('   relative error:%e\n\n', relativeError(nit) ) )
    disp(sprintf('   norm:%e\n\n', norm(nit) ) )
    time = time+dt;     %advance time
end

max(norm)
max(relativeError)





