%% This script is meant to test how temperature affects hydraulic conductivity
% Hydraulic conductivity is computed by the function kappa.m
% Here you can comment/uncomment
% - isothermal hydraulic conductivity (no dependence on temperature)
% - non-isothermal hydraulic conductivity following  Ronan ar al. 1998
% - non-isothermal hydraulic conductivity following  Constantz 1990, Kestin et al. 1978
% The relation psi-K is described by Van Genucthen parametrization
% accordingly with Mualem's model.
% the reference hydraulic conductivity is assumed to be at T=20 C
global alpha thetas thetar n m Ks

Ks      = 0.062/day;    %[meter/second]
thetas  = 0.41;         %[-] saturated water content
thetar  = 0.095;        %[-] residuel water content
n       = 1.31;         %[-] parameter n
m       = 1-1/n;        %[-] parameter m
alpha   = 1.9;

for i=1:200
    % increase the denominator to see what happens near saturation. Here
    % there are the greatest differences
    psi(i)= -i/1000;
    ks(i)=Ks;
    k10(i) = kappa(psi(i),10); % T= 10 celsius
    k15(i) = kappa(psi(i),15); % T= 15 celsius
    k20(i) = kappa(psi(i),20); % T= 20 celsius
    k25(i) = kappa(psi(i),25); % T= 25 celsius
end
    
plot(psi,k10,'b')
hold on
plot(psi,k15,'g')
plot(psi,k15-k10, 'k')

plot(psi,k20,'r')
% plot(psi,k10-k20, 'k')

plot(psi,ks)