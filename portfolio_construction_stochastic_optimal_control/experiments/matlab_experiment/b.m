% Right now b cannot have off diagonal terms which means 
function y = b(x, sig)
    y = {sig(1)*(0.05-x(1)), 0*x(1); 0*x(2), sig(2)*(0.05-x(2))};
end