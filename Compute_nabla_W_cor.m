function [nabla_W_cor]=Compute_nabla_W_cor(i,N,x,m,h,rho,nabla_W_cor)
%kernel of SPH as a function
% of the coordinate particle i and j

% input:  %a number initial particle
          %b namber calculated particle 
          %h  blurring radius
          %x coordinate all particle 
% output: W = the force of the particle j effect on the initial i
sumW=0;
gamma=0;
W_cor_tmp=zeros(2,N,N);
Li=zeros(2,2,N);
diad=zeros(2,2);
for j = 1:N
    sumW=sumW+m/rho(1,j)*ComputeW(i,j,x,h);
end

for j = 1:N
    for beta = 1:2
    gamma=gamma+m/rho(1,j)*Compute_nabla_W(i,j,x,h,beta);
    end
end

gamma=gamma/sumW;

for j = 1:N
    for beta = 1:2
    W_cor_tmp(beta,i,j)=(Compute_nabla_W(i,j,x,h,beta)-gamma)/sumW;
    end
end

for j = 1:N
    diad(1,1)=W_cor_tmp(1,i,j)*x(1,j);
    diad(1,2)=W_cor_tmp(1,i,j)*x(2,j);
    diad(2,1)=W_cor_tmp(2,i,j)*x(1,j);
    diad(2,2)=W_cor_tmp(2,i,j)*x(2,j);
    Li(1:2,1:2,j)=Li(1:2,1:2,j)+m/rho(1,j)*diad;
end

for j = 1:N
    Li(1:2,1:2,j)=Li(1:2,1:2,j)';
end

for j = 1:N
    nabla_W_cor(1:2,i,j)=Li(1:2,1:2,j)'*nabla_W_cor(1:2,i,j);
end

nabla_W_cor=nabla_W_cor+W_cor_tmp;





