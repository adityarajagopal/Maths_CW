%% 2) RLC circuit - Excercise 3
%
% *Introduction:*
% 
% The following script tests the 4th order classic Runge-Kutta for an RLC
% circuit, which is a second order ODE. To simulate the system we use the
% rukasecond.m function, which calculates the next iteration of the numeric
% ODE solution. This function is called in the auxilary N_step_rk(), which
% repeats and saves the output for N number of steps
%
% *Mathematics involved:*
%
% The system is characterized by the following equations:
%
% $$ L \frac{d^2}{dt^2}q_{C}(t) + R\frac{d}{dt}q_{C}(t) + \frac{1}{C} q_{C}(t)  $$
%
% $$ Vout = \frac{d}{dt}q_{C}(t) $$
%
% We first need to rewrite the first equation as two simultaneous first
% order equations to be solved by our Runge-Kutta algorithm. Since we are dealing with 
% variances in charge (the derivative of q in terms of t) this can
% convinently be represented as the current i.
%
% $$ i'=\frac{1}{L}(Vin-Ri-\frac{1}{c}q) $$
%
% $$ q' = i $$
%
% Finally we can also rewrite our voltage output in terms of the current,
% which leads to the more recognizable equation
%   
% $$ Vout = R * i $$


function RLC_script()
%%
% Local function to produce Vout for the current set Vin and conditions for N steps. 
% 
% _Note: since we have set the whole script to be a function we can
% use local functions that hae access to the variables within the
% script. This means it will be evaluated with the value Vin, R, C
% etc have when the function is called_
function [Tout,Vout] = N_step_rk()    
func=@(q, i, t) (1/L)*(v_in(t)-(R*i)-(q/C)); 
t=(0:h:h*(N-1));
for ind = 1:N-1
   [q(ind+1), i(ind+1)]=rukasecond(q(ind), i(ind), t(ind), h, func); 
end
Tout=t;
Vout=R*i;
end

%%
% Set given conditions
R=250;
L=600e-3;
C=3.5e-6;
h=0.0000015;
N=10000;
q=zeros(N,1);
i=zeros(N,1);
q(1)=500e-9;
i(1)=0;

%%
% *Step signal input*
v_in = @(t) 5*heaviside(t);
[Tout,Vout] = N_step_rk();
figure(1);
hold on;
plot(Tout,Vout);

title('Graph of Vout against time for a Step Signal');
xlabel('Time(s)');
ylabel('Vout(V)');
legend('Amplitude Vin = 5V')
%%
% The input is a step function with amplitude 5V. Initially, the Vout is 0V
% as there is no current through the circuit. As the voltage across the
% capacitor cannot change instantaneously, the increase in the voltage is
% gradual initially. As t tends to infinity, the input is a steady 5V,
% hence Vout dies down while oscillating about 0V. The peak voltage it
% reaches is about 2V. 
%%
% *Decay signal input*
v_in = @(t) 5*heaviside(t)*exp(-t^2/(3e-6));
[Tout,Vout] = N_step_rk(); 
figure(2)
plot(Tout,Vout,'-r.');
hold on;

title('Graph of Vout against time for Decay Signal input');
xlabel('Time(s)');
ylabel('Vout(V)');
legend('Amplitude tau = 3(ms)^2');
%%
% The input here is a deccaying exponential, hence also has a transient at
% t=0. The output follows a similar shape to before but doesn't reach as
% high a peak voltage as the exponential input decreases with time as
% opposed to staying constant as in the previous example.
%%
% *Square signal inputs*
f=109;
v_in = @(t) 5*square((2*pi*f*t));
[Tout,Vout] = N_step_rk();
figure(3);
subplot(2, 1, 1)
hold on;
plot(Tout,Vout,'-r.');
subplot(2, 1, 2)
plot(Tout, 5*square((2*pi*f*t)),'-r.');
hold on;


f=5;
v_in = @(t) 5*square((2*pi*f*t));
[Tout,Vout] = N_step_rk();
subplot(2, 1, 1)
plot(Tout,Vout,'-b.');
subplot(2, 1, 2)
plot(Tout, 5*square((2*pi*f*t)),'-b.');

f=500;
v_in = @(t) 5*square((2*pi*f*t));
[Tout,Vout] = N_step_rk();
subplot(2, 1, 1)
plot(Tout,Vout,'-g.');
subplot(2, 1, 2)
plot(Tout, 5*square((2*pi*f*t)),'-g.');

subplot(2, 1, 1)
title('Graph of Vout against time for Square Signal inputs');
xlabel('Time(s)');
ylabel('Vout(V)'); 
legend('f=109hz','f=5hz','f=500hz');

subplot(2, 1, 2)
title('Graph of Vin');
xlabel('Time(s)');
ylabel('Vin(V)'); 
legend('f=109hz','f=5hz','f=500hz');
%%
% DESCRIBE GRAPH

%%
% *SINE SIGNAL*
f=109;
v_in = @(t) 5*sin((2*pi*f*t));
[Tout,Vout] = N_step_rk();
figure(4);
plot(Tout,Vout,'-r.');
hold on;

f=9;
v_in = @(t) 5*sin((2*pi*f*t));
[Tout,Vout] = N_step_rk();
plot(Tout,Vout,'-b.');

f=500;
v_in = @(t) 5*sin((2*pi*f*t));
[Tout,Vout] = N_step_rk();
plot(Tout,Vout,'-g.');

title('Graph of Vout against time of Sine Signal inputs');
xlabel('Time(s)');
ylabel('Vout(V)'); 
legend('f= 109Hz','f = 9Hz', 'f = 500hz');
%%
% DESCRIBE GRAPH

end
