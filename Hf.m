%% Total water depth

function y=Hf(psi)

if (psi<0)
    y = 0;
else
    y = psi;
end