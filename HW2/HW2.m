% %Q2
% 
tic
p = [1; 1];
% you had wrong prices here
fVal = bertrand(p);
iJac = inv(myJac('bertrand', p));     %define Jacbian matrix approximation
%%
%
maxit = 100; 
tol = 1e-6; 
for iter = 1:maxit
    fnorm = norm(fVal);
    fprintf('iter %d: p(1) = %f, p(2) = %f, norm(f(x)) = %.8f\n', iter, p(1), p(2), norm(fVal));
    if norm(fVal) < tol
        break
    end
    d = - (iJac * fVal);
    p = p+d;
    fOld = fVal;
    fVal = bertrand(p);
    u = iJac*(fVal - fOld);
    iJac = iJac + ( (d - u) * (d'*iJac) )/ (d'*u);
end
 v = [2; 2];
   exp_excess=exp(v-p);
   D=exp_excess./(1+sum(exp_excess));
   D0=1-sum(D);
% toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Q3
tic
vq3=[2;2]
f1 = @(paq3,pbq3) [1-paq3+(paq3*exp(vq3(1)-paq3))/(1+exp(vq3(1)-paq3)+exp(vq3(2)-pbq3))];          %equation 1
f2 = @(paq3,pbq3) [1-pbq3+(pbq3*exp(vq3(2)-pbq3))/(1+exp(vq3(1)-paq3)+exp(vq3(2)-pbq3))];           %equation 2
paq3=1.5;                   %set k iteration values
paq3old=1.6;               %set k-1 iteation values
pbq3=1.5;
pbq3old=1.6;               

f1q3old = f1(paq3old,pbq3old);          %find value of functions for initial guesses
f2q3old = f2(paq3old,pbq3old);
% Secant iterations:
tolq3 = 1e-6;
maxitq3 = 100;
for iterq3 =1:maxitq3
    f1Valq3 = f1(paq3,pbq3old);
    f2Valq3 = f2(paq3,pbq3);                %Based on Gauss method immediately uses updated values
   fValq3=[f1Valq3;f2Valq3];
    if norm(fValq3) < tolq3
        break
    else
%         use vector notation like in homework. If you had 100 products,
%         you would get tired writing the update.
        paq3new = paq3 - ( (paq3 - paq3old) ./ (f1Valq3 - f1q3old) ).* f1Valq3;
        pbq3new = pbq3 - ( (pbq3 - pbq3old) ./ (f2Valq3 - f2q3old) ).* f2Valq3;
        paq3old = paq3;
        pbq3old = pbq3;
        paq3    = paq3new;
        pbq3    = pbq3new;
        f1q3old = f1Valq3;
        f2q3old = f2Valq3;
    end
end
pq3=[paq3;pbq3];
  exp_excessq3=exp(vq3-pq3);
   Dq3=exp_excessq3./(1+sum(exp_excessq3));
   D0q3=1-sum(Dq3);
  fprintf('iter %d: paq3 = %.8f, pbq3 = %.8f, f(pq3) = %.8f\n', iterq3, paq3, pbq3, fValq3);
toc

% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Q4
tic
vq4=[2;2]
g = @(pq4) [ (1+exp(vq4(1)-pq4(1))+exp(vq4(2)-pq4(2)))/(1+exp(vq4(2)-pq4(2)));(1+exp(vq4(1)-pq4(1))+exp(vq4(2)-pq4(2)))/(1+exp(vq4(1)-pq4(1)))];     %define the function as a vector based on p^(t+1)=1/1-D(p^t)

pq4 = [1.5;1.5];
tolq4=1e-6;
maxitq4 = 100;
for iterq4 = 1:maxitq4
    nextpq4 = g(pq4);
    if norm(pq4 - nextpq4) < tolq4
        break;
    end
    pq4 = nextpq4;
    disp([iterq4;pq4]);
end

fprintf('iterq4 %d: pq4(1) = %f, pq4(2) = %f, norm(f(x)) = %.8f\n', iterq4, pq4(1), pq4(2), norm(pq4 - nextpq4));
exp_excessq4=exp(vq4-pq4);
   Dq4=exp_excessq4./(1+sum(exp_excessq4));
   D0q4=1-sum(Dq4);
toc
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Q5
tic
for i=1:16              % There are 16 values for v2
    v2=(i-1)/5;        % relating v2 values with i for example v2=1 is the 6th element: i=6
vq5=[2;v2]
f1 = @(paq5,pbq5) [1-paq5+(paq5*exp(vq5(1)-paq5))/(1+exp(vq5(1)-paq5)+exp(vq5(2)-pbq5))];       %all same as question 3
f2 = @(paq5,pbq5) [1-pbq5+(pbq5*exp(vq5(2)-pbq5))/(1+exp(vq5(1)-paq5)+exp(vq5(2)-pbq5))];
paq5=1.5;
paq5old=1.6;
pbq5=1.5;
pbq5old=1.6;

f1q5old = f1(paq5old,pbq5old);
f2q5old = f2(paq5old,pbq5old);
% Secant iterations:
tolq5 = 1e-6;
maxitq5 = 100;
for iterq5 =1:maxitq5
    f1Valq5 = f1(paq5,pbq5old);
    f2Valq5 = f2(paq5,pbq5);
   fValq5=[f1Valq5;f2Valq5];
    if norm(fValq5) < tolq5
        break
    else
        paq5new = paq5 - ( (paq5 - paq5old) ./ (f1Valq5 - f1q5old) ).* f1Valq5;
        pbq5new = pbq5 - ( (pbq5 - pbq5old) ./ (f2Valq5 - f2q5old) ).* f2Valq5;
        paq5old = paq5;
        pbq5old = pbq5;
        paq5    = paq5new;
        pbq5    = pbq5new;
        f1q5old = f1Valq5;
        f2q5old = f2Valq5;
        
    end
end
paq5vec(i)=paq5;       %putting paq5 values in a vector that is building up at each loop to make a vector to be able to plot it vs the vector v2
pbq5vec(i)=pbq5;
v2vec(i)=v2;
end
pq5=[paq5;pbq5];
  fprintf('iter %d: paq3 = %.8f, pbq3 = %.8f, f(pq3) = %.8f\n', iterq5, paq5, pbq5, fValq5);
toc
plot(v2vec,paq5vec)
hold on
plot(v2vec,pbq5vec,'red')
legend('PA','PB')                             
xlabel('VB')




