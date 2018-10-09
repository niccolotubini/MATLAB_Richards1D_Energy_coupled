%% Van Genucthen following Mualem's assumpution

function K=kappa(psi,T)

global alpha thetas thetar n m Ks psic
sat = (Thetaf(psi)-thetar)/(thetas-thetar); %saturation is between 0 1 

%% Isothermal unsaturated conductivity
K = Ks*sqrt(sat)*(1-(1-sat^(1/m))^m)^2;

%% K(psi,T) Ronan ar al. 1998
% mu = 0.00002414*10^(247.8/(T+133.16));
% K = 0.00002414*10^(247.8/(20+133.16))/mu*Ks*sqrt(sat)*(1-(1-sat^(1/m))^m)^2;

%% K(psi,T) Constantz 1990, Kestin et al. 1978
% nu2Overnu1 = (20-T)/(T+96)*(1.2364-1.303*10^(-3)*(20-T)...
%             + 5.7*10^(-6)*(20-T)^2);
% nu2 = 100210*10^(-6)* exp(nu2Overnu1);
% K = Ks*sqrt(sat)*(1-(1-sat^(1/m))^m)^2 * 100210*10^(-6)/nu2;
% 

    