%% RUNGE-Kutta_SCRIPT PLOT

%%CONDITIONS
R=250;
L=600e-3;
C=3.5e-6;
f=109;
h=0.0000008;
N=100000;
q(1)=500e-9;
i(1)=0;


%%
%STEP SIGNAL
step=@(t) 5*heaviside(t);
%IMPULSIVE SIGNAL WITH DECAY
impulsive_decay=@(t) step(t)*exp(-t/3e-3);
%SQUARE WAVE
sin_f=@(t) sin(2*pi*f*t);
%SINE WAVE
sqr_f=@(t) heaviside(sin(2*pi*f*t));

v_in=step;          %Assign one of above signals to v_in

q=zeros(N,1);
i=zeros(N,1);
in=zeros(N,1);

%t=zeros(500,1);
t=(0:h:h*(N-1));    %Vector from 0 to h*(N-1) in in increments of h
t=t';               %Transpose vector


%CLASSIC FOURTH-ORDER RUNGE-KUTTA ALGORITHM 
func=@(q, i, t) (1/L)*(v_in(t)-(R*i)-(q/C)); 
for ind = 1:N-1
   [q(ind+1), i(ind+1)]=rukasecond(q(ind), i(ind), t(ind), h, func); 
   in(ind)=v_in(t(ind));
end

%plot(sin(2*pi*109*t));
%VOLTAGE PLOT
subplot(3, 1, 1);
plot(t, q);
title('Graph of Vout against time of Signal');
xlabel('Time (s)');
ylabel('Vout (V)');
%CURRENT PLOT
subplot(3, 1, 2);
plot(t, i);
title('Graph of Iout against time of Signal');
xlabel('Time (s)');
ylabel('Iout (A)');
%INPUT SIGNAL PLOT
subplot(3, 1, 3);
plot(t, in);
% plot(t, q, t, i);