%%
% *MIDPOINT METHOD*

function [Tout,Vout] = midpoint(qc0,h,tf,Vin,R,C) 
t = 0;
Nsteps = round((tf-t)/h); % number of steps to take
Tout = t:h:tf;
q = zeros(1,Nsteps+1);
f = @(t,q) (1/R)*(Vin(t) - (q/C));

% store intial condition 
q(1) = qc0;

%%
for i = 1:Nsteps
    k1 = h*f(Tout(i),q(i));                 
    k2 = h*f(Tout(i)+0.5*h,q(i)+0.5*k1);    
    q(i+1) = q(i) + k2;   
end

Vout = q./C;

end
