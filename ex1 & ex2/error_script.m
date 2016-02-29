function error_script()

%%CONDITIONS
R = 1000;
C = 100e-9;
qc0 = 500e-9;
t0 = 0;
tf = 0.001;
h  = 0.000005;

T = 100e-6    ;
Vin = @(t)5*cos(2*pi*t/T);

%%MIDPOINT METHOD
[Tout,Vout] = midpoint(qc0,h,tf,Vin,R,C); 
figure(1)
plot(Tout,Vout,'-r');
grid on
hold on 


%%EXACT SOLUTION
Vamp = 5;
[Texact,Vexact] = exactcosine(qc0,h,tf,Vamp,R,C,T); 
figure(1);
plot(Texact,Vexact,'-b');
grid on

title('Comparison between Midpoint and exact solution for h=5us') 
xlabel('Time(s)');
ylabel('Vout(V)');
legend('Midpoint','Exact');

%% Error Analysis

error_iter = 5;
error = zeros(1, error_iter);
h=zeros(1, error_iter);
CM=jet(error_iter);

figure(2);
hold on;
axis([3.222e-4, 3.225e-4,0.9796 0.9804])
for l = 1:1:error_iter
    h(l)=0.00001*2^-(l+2);
    [Tout,Vout] = midpoint(qc0,h(l),tf,Vin,R,C); 
    [Texact,Vexact] = exactcosine(qc0,h(l),tf,Vamp,R,C,T);
    error(l)=max(abs(Vout-Vexact));
    plot(Tout, Vout, 'color', CM(l,:));
end
plot(Texact, Vexact, 'color', 'k');
title('Supperposed approaches of solution');
xlabel('Time(s)');
ylabel('Vout(V)');
legend('h=12.5us','h=6.25us','h=3.13us','h=1.56us','h=0.78us','exact solution');

figure(3);
plot(log(h), log(error), 'b*');
title('Error Analysis for given Cosine') 
xlabel('h value(s)');
ylabel('error(V)')

function [Texact,Vexact] = exactcosine(qc0,h,tf,Vin,R,C,T)
    t = 0;
    Nsteps = round((tf-t)/h); %% number of steps to take
    Texact = t:h:tf;
    Vexact(1:Nsteps) = 0;

    %%store intial condition 
    Vexact(1) = qc0/C;

    w = 2*pi/T;
    k = qc0 - (Vin*C)/(1+(w*R*C)^2);
    %%
    for i =1:Nsteps
        qc = (Vin*C)/(1+(w*R*C)^2)*cos(w*t) + (w*Vin*R*C^2)/(1+(w*R*C)^2)*sin(w*t) + k*exp(-t/(R*C)) ; 
        t = t+h; 
        Vexact(i+1) = qc/C;
    end
end
end
