%% Test Diffusion.m for non-linear conduction equation
% Page 105 Dall'Amico Ph.D. thesis
% analytical solution Dumbser et al. 2008
% termal conductivity as function of temperature

function lambda = LambdaT(T)

global epsilon
lambda = epsilon*T*(1-T);