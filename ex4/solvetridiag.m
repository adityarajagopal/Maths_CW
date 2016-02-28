function x = solvetridiag(N,a,b,c,vec)
%**************************************************************
% Use this function to solve a tridiagonal matrix equation
% of the form: Tx=v;
% T is a tridiagonal matrix with entries [a, b, c]
% on the three non-zero diagonals
% v is the right hand side vector
% N is the number of equations
%**************************************************************
diag=ones(N,1)*b;
x=zeros(N,1);
for k = 1:N-1
 p=diag(k);
 m=a/p;
 diag(k+1)=diag(k+1)-m*c;
 vec(k+1)=vec(k+1)-m*vec(k);
end
 x(N)=vec(N)/diag(N);
for k = N-1:-1:1
 x(k)=(vec(k)-c*x(k+1))/diag(k);
end

return
