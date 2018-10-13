function [init] =initialization_v(N,sqn)

v=zeros(2,N);
    for i=1:(sqn*sqn)
        v(1,i)= fix((i-1)/sqn);
        v(2,i)=0;
    end
    

init=v;