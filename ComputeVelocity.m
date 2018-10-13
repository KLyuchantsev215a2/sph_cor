function [velocity]=ComputeVelocity(i,j,beta,dt,m,v,rho,SIG,nabla_W0)
   
v(1,i)=v(1,i)-dt*(m*(SIG(1,beta,i)/rho(i)^2+SIG(1,beta,j)/rho(j)^2))*nabla_W0(beta); 
v(2,i)=v(2,i)-dt*(m*(SIG(2,beta,i)/rho(i)^2+SIG(2,beta,j)/rho(j)^2))*nabla_W0(beta);
    
velocity=v;