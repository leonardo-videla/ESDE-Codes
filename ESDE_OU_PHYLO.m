function [n, t, sol]=ESDE_OU_PHYLO (N, a, m, type)

  % ESDE based on the paper: Using the Ohrstein-Uhlenbeck process to model 
  %the evolution of interacting populations (Bartoszek, Glemin, Kaj, Lascoux 2017).
  %Migration Matrix
 
  
  %m=0.02;
  A=m*eye(N);
  %a=0.08;
  
  x0=zeros(N, 1);
  x0(1, 1)=0;
  sigma=0.2*eye (N);
  h=0.005;
  n=ones(N, 1);
  len=zeros(N,1);
 % theta=0.6*ones(N,1);
 
 vect=5*rand(N, 1);
  
  aver=zeros(N, N);
  for k=1:N
      aver(:, k)= vect;
  endfor
 
  %The chain starts with one species. 
  len(1,1)= (1/(mu_phylo(1)+nu_phylo(1)))*rande(1);
  
  for j=2:N
    p=rand;
    ind=n(j-1, 1);
    q=nu_phylo(ind)/(mu_phylo(ind)+nu_phylo(ind));
    if (p<q)
      n(j,1)=ind-1;
    else
       n(j,1)=ind+1;
    endif
    ind=n(j,1);
    len(j)=(1/ (mu_phylo(ind)+nu_phylo(ind)))*rande(1);
  endfor

  TT=sum(len);
  t=[0:h:floor(TT/h)];
  sol=zeros(N, length(t));

  aux=0;
  C=colormap(jet(max(n)-min(n)+1));
  close all;
  figure 1;
  hold on;
  colormap(C);
 
  for j=1:N
    
    [t1, s]=Multi_OU_Phylo (x0(1:n(j)), migration_matrix (n(j), m),  aver(1:n(j), n(j)),a, sigma(1:n(j), 1:n(j)), len(j), h);
    
    nn=length(t1);
    sol(1:n(j), aux+1:aux+nn)=s;
    
    %The last state of the chain. 
    last=s(n(j),nn);
    
    if (j < N)
      %Here we encode the transfer kernel.  
      if (n(j)<n(j+1))
        x0(1:n(j), 1)=s(1:n(j), nn);
        if type==1
            ind=floor(rand(1)*(n(j)-1))+1; %USe this in order to choose the ancestor uniformly.
        else 
            ind=choose_closest(n(j+1), aver);
          %Use this in order to choose the ancestor with probability proportional to the distance between the corrsponding long-term averages
        endif 
        ind       
        x0(n(j+1), 1)= s(ind, nn)*(1+ 0.01*(randn(1)));
      endif  
    endif
    
    %Plotting the tree. 
   
    abscisa=([0:1:nn-1]+(aux)*ones(1, nn))*h;    
    for k=1:n(j)
      plot(s(k, :), abscisa, "color", C(n(j), :));
    endfor
    
    aux=aux+nn;
    
  endfor
  
  
 
  xlabel("Trait value");
  ylabel("Time");
  title("Evolution of the traits under speciation");
  colormap(C);
  q=colorbar();
  ylabel(q, "Number of mutants");
  caxis([0 max(n)]);
  hold off;
  
  
endfunction

