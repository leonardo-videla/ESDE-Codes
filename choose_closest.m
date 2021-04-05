function ind=choose_closest (n, aver)
  ind=1;
  alpha=20;
  if (n==2)
     return;
  else
    v=zeros(n-1, 1);
    past=aver (1:n-1, n-1);
    present=aver (1:n, n);
    for j=1:n-1 
        v(j,1)=exp(-alpha* abs(present(n)-past(j)));
    endfor
    s=sum(v);
    v=v/s;
    for j=2:n-1
     v(j,1)=v(j,1)+v(j-1,1);
    endfor
    p=rand(1);
    j=1;
    while p>v(j,1)
      j++;
    endwhile 
    ind=j;   
  endif
endfunction
