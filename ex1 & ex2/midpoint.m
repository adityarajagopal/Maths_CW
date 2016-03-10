%%
% *MIDPOINT METHOD*

function [Tout,Vout] = midpoint(qc0,h,tf,Vin,R,C) 
t = 0;
Nsteps = round((tf-t)/h); % number of steps to take
Tout = t:h:tf;
q = zeros(1,Nsteps+1);

%%
% We first have to write the given equation for V_in in terms of q':
%
% $$ V_in(t)=\frac{1}{C}q_{c}(t) + R \dot{q}_{C}(t) $$ 
%
% is equivalent to
%
% $$ \dot{q}_{C}(t) = \frac{1}{R}(V_in(t) - \frac{1}{C}q_c(t)) $$

f = @(t,q) (1/R)*(Vin(t) - (q/C));

% store intial condition 
q(1) = qc0;

%%
% To implement the midpoint method we must find the coefficients $$k_1$ and $$k_2$
% of the increment function $$ \phi=a_1*k_1+a_2*k_2 $$, having the midpoint
% method $$a_1=0$ and $$a_2=1$. We still have to compute $$k_1$ since
% $$k_2$ is defined in terms of $$k_1$. Then the value of the next
% iteration will be $$y_{i+1}=y_i+h*\phi=y_i+h*k_2$$
%%
for i = 1:Nsteps
    k1 = f(Tout(i),q(i));                 
    k2 = f(Tout(i)+0.5*h,q(i)+h*0.5*k1);    
    q(i+1) = q(i) + h*k2;   
end

Vout = q./C;

end
