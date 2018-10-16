function [Ln1]=ComputeL(L,i,j,beta,m,v,rho,nabla_W)

L(1,beta,i)=L(1,beta,i)+(m/rho(1,j)*(v(1,j)-v(1,i))*nabla_W);  
L(2,beta,i)=L(2,beta,i)+(m/rho(1,j)*(v(2,j)-v(2,i))*nabla_W); 

Ln1=L;