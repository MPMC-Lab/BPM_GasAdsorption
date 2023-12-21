%% Langmuir isotherm
function g =Isotherm_Langmuir(c,b,qm)
g= qm*(b*c)./(1+b*c);
end