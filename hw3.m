load hw3.mat
numaff=y(:);
cnsttrm=X(:,1);
age=X(:,2);
years_married=X(:,3);
religiousness=X(:,4);
occupation=X(:,5);
selfrate=X(:,6);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n=length(numaff);
facy=(factorial(numaff));
lnfacy=log(facy(:,1));
for i=1:n
    yx(i,:)=y(i,:).*X(i,:);
end
tic
%Likelihood=@(beta) -exp(X*beta)+yx*beta-lnfacy;
Likelihood=@(beta)-sum(-exp(X*beta)+yx*beta-lnfacy);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%Q1 WITH FMINSEARCH

beta0=[7;-0.2;0.5;-1;0.3;-2];
%beta0=[0.5;0.5;0.5;0.5;0.5;0.5];
 betaq1fminsearch=fminsearch(Likelihood,beta0);
 toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%

%  %    Q1 WITH NELDER MEAD
%  tic
toler=1e-6;
beta0=[7;-0.2;0.5;-1;0.3;-2];
%beta0=[0.5;0.5;0.5;0.5;0.5;0.5];
maxim_feval=100;
Likelihood=@(beta)-sum(-exp(X*beta')+yx*beta'-lnfacy); 
[betaq1nelder,minusmaxLikelihood,numoffuncs]=ANMS(Likelihood,beta0,toler,maxim_feval);
toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Q2 WITH FMINUNC
tic
beta0=[7;-0.2;0.5;-1;0.3;-2];
%beta0=[0.5;0.5;0.5;0.5;0.5;0.5];
Likelihood=@(beta)-sum(-exp(X*beta)+yx*beta-lnfacy);
[betaq2,fval,exitflag,output,grad,hessian] = fminunc(Likelihood,beta0)
toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Q2 WITH STEPPEST DESCENT
% Likelihood=@(beta)-sum(-exp(X*beta)+yx*beta-lnfacy);
% 
% beta0=[7;-0.2;0.5;-1;0.3;-2];
% gf=@(beta)-X'*exp(X*beta)+X*y;  
% n=1;  
% while(norm( gf(beta))>0.05)     
%     beta= beta-0.01*(1/n) *gf(beta);    
%     n=n+1;
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Q3
tic
sumresfun=@(beta)sum((y-exp(X*beta)).^2);
[betaq3,resnorm,residual,exitflag,output,lambda,jacobian] = lsqnonlin(sumresfun,beta0);
toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  q4 WITH FMINSEARCH
tic
maxit=100;
beta0=[7;-0.2;0.5;-1;0.3;-2];
%beta0=[0.5;0.5;0.5;0.5;0.5;0.5];

sumres=@(beta)sum((y-exp(X*beta)).^2);
betaq4fminsearch=fminsearch(sumres,beta0);
toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Q4 WITH NELDER MEAD
tic
toler=1e-6;
beta0=[7;-0.2;0.5;-1;0.3;-2];
%beta0=[0.5;0.5;0.5;0.5;0.5;0.5];
maxim_feval=100;
sumres=@(beta)sum((y-exp(X*beta')).^2); 
[betaq4nelder,minusmaxLikelihood,numoffuncs]=ANMS(sumres,beta0,toler,maxim_feval);
toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







