function proc = ESDE_OU_PHYLO_2 (N, a, m)

  % ESDE with state-dependent rate of branching, extending the paper: Using the Ohrstein-Uhlenbeck process to model 
  %the evolution of interacting populations (Bartoszek, Glemin, Kaj, Lascoux 2017).
  %Migration Matrix
 
  
  %m=0.02;
  %a=0.08;
  vect=5*rand(N, 1);
  
  aver=zeros(N, N);
  for k=1:N
      aver(:, k)= vect;
  endfor
  x0=zeros(N, 1);
  x0(1, 1)=0;
  sigma=0.2*eye (N);
  h=0.005;
  len=zeros(N,1);
  LAMBDA_MAX=10;
  tiempos=zeros(N, 1);
  %theta=0.6*ones(N,1);
 
  %The chain starts with one species. 
  
  
 
  
  j=1;
  T=(1/LAMBDA_MAX)*rande(1);
  [t1, s]=Multi_OU_Phylo (x0(1:1), migration_matrix(1, m), aver(1:1, 1), a, sigma(1:1, 1:1), T, h); 
  len=length(t1);
  proc=zeros(N, len);
  proc(1, 1:len)=s; 
  tiempos (1)=len;
  x0(1)=s(:, len);
 
  if ( rand (1) < lambda_dependent (1,x0, aver(1:1, 1), migration_matrix(1,m), a)/LAMBDA_MAX)  
    j++;
    x0(2)=x0(1)*0.9; 

  endif
  
  while j<=N
    T=(1/LAMBDA_MAX)*rande(1);
    [t1, s]=Multi_OU_Phylo (x0(1:j), migration_matrix (j, m),  aver(1:j, j),a, sigma(1:j, 1:j), T, h);
    len=length(t1);
    aux=zeros(N, len);
    aux(1:j, :)=s;
    proc = [proc aux];
    tiempos (j)=tiempos (j)+len;
    x0(1:j)=s(:, len);
    if ( rand (1) < lambda_dependent (j, x0, aver(1:j, j),migration_matrix(j,m), a)/LAMBDA_MAX)  
      j++;
      ind=floor(rand(1)*(j-1))+1 %USe this in order to choose the ancestor uniformly
      x0(j)=x0(ind)*0.9;
    endif
  endwhile
 
  figure;
  hold on;
  C=colormap(jet(N));
  colormap(C); 
  xlabel("Trait value");
  ylabel("Time");
  vect
  title("Evolution of the traits under speciation.");
  colormap(C);
  q=colorbar();
  ylabel(q, "Number of mutants");
  caxis([0 N]);
  tiempos
  t0=0;
  
  for j=1:N
   abscisa=[t0:1:(t0+tiempos(j)-1)]*h;
   for k=1:j
      plot(proc(k, (t0+1:t0+tiempos(j))), abscisa, "color", C(j, :));
   endfor
   t0=t0+tiempos(j);
  endfor
  
  hold off;
  
  
endfunction

