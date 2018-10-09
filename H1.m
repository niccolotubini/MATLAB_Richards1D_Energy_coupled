%% Jordan decomposition: Hf = H1-H2

function y = H1(psi)
if (psi<0)
    y = 0;
else
    y = Hf(psi);    
end
