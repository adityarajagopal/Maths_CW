%initialising parameters
h = 0.02; 
A = 1;
B = -1; 
C = -2; 

t0 = 0; 
tf = 1; 
x0 = 1; 
xf = 2; 

func = @(t) t; %this is the forcing function f(t) on the rhs of the differential equation

%solving using method of finite differences 
a = (A/h^2)-(B/(2*h)); %coefficient of x(ti-1)
b = C - ((2*A)/h^2); %coefficient of x(ti)
c = (A/h^2) + (B/(2*h)); %coefficient of x(ti+1)

N = (tf-t0)/h; %number of steps int he interval

F = zeros(N-1,1); %matrix on rhs of eqn

%setting up F
for i = 1:(N-1)
    F(i) = func(i);
end
F(1) = F(1) - (a*x0);
F(N-1) = F(N-1) - (a*xf);

%calling solvetridiag
U = solvetridiag(N-1,a,b,c,F);

%plotting 
U = [x0;U;xf]; 
t = [t0:h:tf];
figure;
plot(t,U,'b*');


