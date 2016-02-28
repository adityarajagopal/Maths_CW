
%% MIDPOINT METHOD
function [Tout,Vout] = midpoint(qc0,h,tf,Vin,R,C) 
t = 0;
Nsteps = round((tf-t)/h); %% number of steps to take
%Tout = zeros(Nsteps+1:1);
Tout = t:h:tf;
%Vout(1:Nsteps+1) = 0;
q = zeros(1,Nsteps+1);
f = @(t,q) (1/R)*(Vin(t) - (q/C));
%% store intial condition 
%qc = qc0;
%Vout(1) = qc0/C;
q(1) = qc0;
%%
%{
for i =1:Nsteps+1
%grad = gradient(f,time,y,R,C); %% evaluate the initial derivatives 
grad = ( Vin(t) - qc/C ) /R;

qcH = qc + grad*h/2; %% take Euler step to midpoint
timeH = t+h/2;

%dy = gradietn(f,(time+dt/2),(y + dy*dt/2),R,C); %% re-evaluate the derivs
grad = ( Vin(timeH) - qcH/C ) /R;

qc = qc + grad*h; 
t = t+h; 

Tout(i+1) = t; 
Vout(i+1) = qc/C;
end
%}
%%
for i = 1:Nsteps
    k1 = h*f(Tout(i),q(i));
    k2 = h*f(Tout(i)+0.5*h,q(i)+0.5*k1);
    q(i+1) = q(i) + k2;   
end

Vout = q./C;

end

%%f(x,y) = dy/dx
%function f = gradient(Vin,t,y,R,C)
%f = ( Vin(t) - y/C ) /R;
%end