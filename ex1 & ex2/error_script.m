


function error_script()

%%CONDITIONS
R = 1000;
C = 100e-9;
qc0 = 500e-9;
%qc0 = 0;
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
plot(Texact,Vexact,'-b');
grid on

title('Comparison between Midpoint and exact solution') 
xlabel('Time(s)');
ylabel('Vout(V)');
legend('Midpoint','Exact');

%%ERROR%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



end


function [Texact,Vexact] = exactcosine(qc0,h,tf,Vin,R,C,T)
t = 0;
Nsteps = round((tf-t)/h); %% number of steps to take
%Texact = zeros(Nsteps+1,1);
%Vexact = zeros(Nsteps+1,length(qc0));
Texact = t:h:tf;
Vexact(1:Nsteps) = 0;

%% store intial condition 
Vexact(1) = qc0/C;

tau = 1/(R*C);
w = 2*pi/T;
k = qc0 - (Vin*C)/(1+(w*R*C)^2);
%%
for i =1:Nsteps
qc = (Vin*C)/(1+(w*R*C)^2)*cos(w*t) + (w*Vin*R*C^2)/(1+(w*R*C)^2)*sin(w*t) + k*exp(-t/(R*C)) ; 
t = t+h; 
Vexact(i+1) = qc/C;
end
end
