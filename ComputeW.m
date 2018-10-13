function [W]=ComputeW(i,j,x,h,m,rho)
%kernel of SPH as a function
% of the coordinate particle i and j

% input:  %a number initial particle
          %b namber calculated particle 
          %h  blurring radius
          %x coordinate all particle 
% output: W = the force of the particle j effect on the initial i

W=0;
r=zeros(2);

r(1)=x(1,i)-x(1,j);
r(2)=x(2,i)-x(2,j);
beta=zeros(2);
W_old=0;

 for i = 1:N
     W_old=ComputeW_old(i,j,x,h);
     beta=beta+m/rho(i)*kron(r,r)*W_old*m/rho(i)*r*W_old;
 end
 
q=norm(r,2)/h;
C=1/(pi*h*h);
if (q>=0) && (q<=1)
    W=C*(15 / 7)*(2 / 3 - q*q + 1/2*q*q*q);
elseif (q >= 1) && (q <= 2)
    W = C*(5 / 14)*(2 - q)^3;
end