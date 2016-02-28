function out = heaviside(t)
    if t>0
        out = 1;
    elseif t==0
        out = 0.5;
    else 
        out = 0;
    end
end