


function error_script()

%%CONDITIONS
R = 1000;
C = 100 * 10^-3;
qc0 = 500 * 10^-9;
%qc0 = 0;
t0 = 0;
tf = 0.001;
h  = 0.00001;

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

%%
for i =1:Nsteps+1
    
qc = (Vin/R)*(1/(tau^2 + w^2))*(tau*cos(w*t) + w*sin(w*t) + (qc0*R/Vin*(tau^2 + w^2)-tau)*exp(-t*tau));

t = t+h; 

Texact(i+1) = t; 
Vexact(i+1) = qc/C;

end
end
