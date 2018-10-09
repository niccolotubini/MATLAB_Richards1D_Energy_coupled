function y = CT(theta)
global thetas rho cW
cGS = (2.128*0.6+2.385*0.4)/(0.6+0.4)*10^6;
y = cGS*(1-thetas) + rho*cW*theta;