N1=100;
N2=1000;
N3=10000;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%Q3
x_100=qnwequi(N1,0,1,'N');
y_100=sqrt(1-x_100.^2);
pi3_100=4*sum(y_100)/N1;
seQMC100=(pi3_100-pi)^2;


x_1000=qnwequi(N2,0,1,'N');
y_1000=sqrt(1-x_1000.^2);
pi3_1000=4*sum(y_1000)/N2;
seQMC1000=(pi3_1000-pi)^2;



x_10000=qnwequi(N3,0,1,'N');
y_10000=sqrt(1-x_10000.^2);
pi3_10000=4*sum(y_10000)/N3;
seQMC10000=(pi3_10000-pi)^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%Q4
pi4_100=Newton_Cotes(@(x) 4*sqrt(1-x.^2),0,1,N1);
pi4_1000=Newton_Cotes(@(x) 4*sqrt(1-x.^2),0,1,N2);
pi4_10000=Newton_Cotes(@(x) 4*sqrt(1-x.^2),0,1,N3);
seNC100=(pi4_100-pi)^2;
seNC1000=(pi4_1000-pi)^2;
seNC10000=(pi4_10000-pi)^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%Q5
num_sim=200;
seed=1234567890;
rng(seed);
for i=1:num_sim
    pi5_100(i)=4*sum(sqrt(1-rand(N1,1).^2))/N1;
end
mse5_100=sum((pi5_100-pi).^2)/num_sim;
mean5_100=mean(pi5_100);


for i=1:num_sim
    pi5_1000(i)=4*sum(sqrt(1-rand(N2,1).^2))/N2;
end
mse5_1000=sum((pi5_1000-pi).^2)/num_sim;
mean5_1000=mean(pi5_1000);



for i=1:num_sim
    pi5_10000(i)=4*sum(sqrt(1-rand(N3,1).^2))/N3;
end
mse5_10000=sum((pi5_10000-pi).^2)/num_sim;
mean5_10000=mean(pi5_10000);

