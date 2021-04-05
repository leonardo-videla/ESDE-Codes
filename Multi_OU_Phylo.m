function [t, sol]=Multi_OU_Phylo (x0, A, theta, a, sigma, T, h)

%Solves the multidimensional OU process

t=[0:h:h*floor(T/h)];
N=size(t, 2);
n=size(x0, 1);
sol=zeros(n,N);

%Auxiliary definitions
H=a*eye(n)-A;


for i=1:n
  sol(i, 1)=x0(i);
endfor

for j=2:N
  %This line is Euler method.
  v=randn(n,1);
  sol(:,j)=sol(:, j-1) + h*(-H*sol(:, j-1)+a*theta) +sigma*sqrt(h)*v;   
  
endfor

endfunction