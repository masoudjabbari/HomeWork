function [Value]=Newton_Cotes(f,a,b,N)
assert(a<b);
assert(floor(N)==N);
assert(N>0);
if(mod(N,2)==1)
    N=N+1;
end
    h=(b-a)/N;
    X=a:h:b;
    w=2*ones(length(X),1);
    evens=2:2:N;
    w(evens)=4;
    w(1)=1;
    w(end)=1;
    w=w*(h/3);
    fv=f(X);
    Value=fv*w;
    
end
