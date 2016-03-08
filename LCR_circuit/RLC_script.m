%% RUNGE-Kutta_SCRIPT PLOT
function RLC_script()
%%CONDITIONS
R=250;
L=600e-3;
C=3.5e-6;
h=0.0000008;
N=10000;
q(1)=500e-9;
i(1)=0;


%%
%STEP SIGNAL
v_in = @(t) 5*heaviside(t);
[Tout,Vout] = N_step_rk();
figure(1);
plot(Tout,Vout);

title('Graph of Vout against time for a Step Signal');
xlabel('Time(s)');
ylabel('Vout(V)');
legend('Amplitude Vin = 5V')

%%
%DECAY SIGNAL
v_in = @(t) 5*heaviside(t)*exp(-t^2/(0.003));
[Tout,Vout] = N_step_rk(); 
figure(2)
plot(Tout,Vout,'-r.');

title('Graph of Vout against time for Decay Signal input');
xlabel('Time(s)');
ylabel('Vout(V)');
legend('Amplitude tau = 3ms');

%%
%SQUARE SIGNAL
f=109;
v_in = @(t) 5*square((2*pi*f*t));
[Tout,Vout] = N_step_rk();
figure(3);
plot(Tout,Vout,'-r.');
hold on

f=5;
v_in = @(t) 5*square((2*pi*f*t));
[Tout,Vout] = N_step_rk();
plot(Tout,Vout,'-b.');
hold on

f=500;
v_in = @(t) 5*square((2*pi*f*t));
[Tout,Vout] = N_step_rk();
plot(Tout,Vout,'-g.');
hold on


title('Graph of Vout against time for Square Signal inputs');
xlabel('Time(s)');
ylabel('Vout(V)'); 
legend('f=109hz','f=5hz', 'f=500hz');

%%
%SINE SIGNAL
f=109;
v_in = @(t) 5*sin((2*pi*f*t));
[Tout,Vout] = N_step_rk();
figure(4);
plot(Tout,Vout,'-r.');
hold on

f=9;
v_in = @(t) 5*sin((2*pi*f*t));
[Tout,Vout] = N_step_rk();
plot(Tout,Vout,'-b.');
hold on

f=500;
v_in = @(t) 5*sin((2*pi*f*t));
[Tout,Vout] = N_step_rk();
plot(Tout,Vout,'-g.');
hold on

title('Graph of Vout against time of Sine Signal inputs');
xlabel('Time(s)');
ylabel('Vout(V)'); 
legend('f=hz','f=Hz', 'f=hz');

%%

%CLASSIC FOURTH-ORDER RUNGE-KUTTA ALGORITHM FOR N STEPS
function [Tout,Vout] = N_step_rk()
    q=zeros(N,1);
    i=zeros(N,1);
    func=@(q, i, t) (1/L)*(v_in(t)-(R*i)-(q/C)); 
    t=(0:h:h*(N-1));
    for ind = 1:N-1
       [q(ind+1), i(ind+1)]=rukasecond(q(ind), i(ind), t(ind), h, func); 
    end
    Tout=t;
    Vout=R*i;
end
end
