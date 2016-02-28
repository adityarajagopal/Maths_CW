
%% MIDPOINT METHOD
function [Tout,Vout] = midpoint(qc0,h,tf,Vin,R,C) 
t = 0;
Nsteps = round((tf-t)/h); %% number of steps to take
%Tout = zeros(Nsteps+1:1);
Tout = t:h:tf;
Vout(1:Nsteps+1) = 0;

%% store intial condition 
qc = qc0;
Vout(1) = qc0/C;

%%
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
end

%%f(x,y) = dy/dx
%function f = gradient(Vin,t,y,R,C)
%f = ( Vin(t) - y/C ) /R;
%end