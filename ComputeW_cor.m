function [W_cor]=ComputeW_cor(i,N,x,m,h,rho,W_cor)
%kernel of SPH as a function
% of the coordinate particle i and j

% input:  %a number initial particle
          %b namber calculated particle 
          %h  blurring radius
          %x coordinate all particle 
% output: W = the force of the particle j effect on the initial i
sumW=0;
W_cor_tmp=zeros(N,N);
 
for j = 1:N
    sumW=sumW+m/rho(1,j)*ComputeW(i,j,x,h);
end

for j = 1:N
    W_cor_tmp(i,j)=ComputeW(i,j,x,h)/sumW;
end

W_cor=W_cor+W_cor_tmp;





