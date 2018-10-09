%% Van Genucthen SWRC Model

function y=Thetaf(psi)
global alpha thetas thetar n m
if (psi<=0)
    y = thetar + (thetas-thetar)/( (1+abs(alpha*psi)^n)^m );
else
    y = thetas;
end

