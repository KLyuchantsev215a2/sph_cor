function [velocity]=ComputeVelocity(i,j,beta,dt,m,v,rho,SIG,nabla_W0)
   
v(1,i)=v(1,i)+(m/rho(1,j)*v(2,j))*nabla_W0(beta); 
v(2,i)=v(2,i)+(m/rho(1,j)*v(2,j))*nabla_W0(beta);
    
velocity=v;