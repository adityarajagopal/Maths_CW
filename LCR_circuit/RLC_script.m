R=250;
L=600e-3;
C=3.5e-6;

f=109;

step=@(t) 5*heaviside(t);
impulsive_decay=@(t) step(t)*exp(-t/3e-3);
sin_f=@(t) sin(2*pi*f*t);
sqr_f=@(t) heaviside(sin(2*pi*f*t));

v_in=step;

h=0.0000008;
N=100000;
q=zeros(N,1);
i=zeros(N,1);
in=zeros(N,1);
q(1)=500e-9;
i(1)=0;
%t=zeros(500,1);
t=(0:h:h*(N-1));
t=t';
func=@(q, i, t) (1/L)*(v_in(t)-(R*i)-(q/C));

for ind = 1:N-1
   [q(ind+1), i(ind+1)]=rukasecond(q(ind), i(ind), t(ind), h, func); 
   in(ind)=v_in(t(ind));
end
%plot(sin(2*pi*109*t));
subplot(3, 1, 1);
plot(t, q);
subplot(3, 1, 2);
plot(t, i);
subplot(3, 1, 3);
plot(t, in);
% plot(t, q, t, i);