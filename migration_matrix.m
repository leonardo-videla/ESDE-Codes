function A=migration_matrix (n, m)
  A=zeros(n,n);
  if n==1
      A=m;
  else
      A=-(m/(n-1))*(ones(n,n)-eye(n)) + m*eye(n);
  endif
endfunction
