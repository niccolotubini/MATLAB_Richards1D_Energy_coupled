%% Compute thermal conductivity as a function of water content

function y=Lambda(theta)
global thetas thetar
lambdaG = 1.5;
lambdaW = 0.6;
lambdaSat = lambdaG^(1-thetas)*lambdaW^(theta);
lambdaDry = (0.135*(2700*(1-thetas))+64.7)/(2700-0.947*(2700*(1-thetas)));

%Ke = log(theta/thetas)+1;
Ke = theta/thetas;
y = Ke*lambdaSat+(1-Ke)*lambdaDry;