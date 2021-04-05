function l=lambda_dependent(j, x, aver, A, a)
  H=a*eye(j)-A;
  l=1;
  v=x(1:j);
  r=aver;
  l=min(max(0, sqrt(sum ((H*v-a*r(1:j))'*(H*v-a*r(1:j)))))/log(1+j),10);
endfunction
