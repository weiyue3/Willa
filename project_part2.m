function project_part2
N_points=50000;
CA_t0=0.186;%CA0=mol/L initial conc. of sugar
CB_t0=0; %CB0=mol/L initial conc. of ethanol
CC_t0=46.18; %CC0=g/L initial conc. of yeast 

trange = [0 300];%unit of minutes
dt = (trange(2)-trange(1))/N_points; %dt in Euler method, in minutes
CA=zeros(N_points,1);
CB=zeros(N_points,1);
CC=zeros(N_points,1);
t=zeros(N_points,1); %minutes
%Set initial conditions
CA(1,1)=CA_t0; %mol/L
CB(1,1)=CB_t0; %mol/L
CC(1,1)=CC_t0; %g/L

for i=1:N_points
    if CA(i,1)>0
        t(i+1,1)=t(i,1)+dt;
        C=[CA(i,1);CB(i,1);CC(i,1)];
        dCdt=diff_eqs(C);
        CA(i+1,1)=CA(i,1)+dCdt(1,1)*dt;
        CB(i+1,1)=CB(i,1)+dCdt(2,1)*dt;
        CC(i+1,1)=CC(i,1)+dCdt(3,1)*dt;
    else
        t(i+1,1)=t(i,1)+dt;
        CA(i,1)=0;
        C=[CA(i,1);CB(i,1);CC(i,1)];
        dCdt=diff_eqs(C);
        CA(i+1,1)=CA(i,1);
        CB(i+1,1)=CB(i,1);
        CC(i+1,1)=CC(i,1)+dCdt(3,1)*dt;
    end
end

figure
hold on
plot(t(:,1),CA(:,1))
plot(t(:,1),CB(:,1))
%plot(t(:,1),CC(:,1))
%title('Concentration Profile of Glucose and ethanol with constant number of yeast cells')
title('Concentration Profile of Glucose and ethanol with growth of yeast cells')
%title('Concentration Profile of yeast cells')
xlabel('time (min)')
ylabel('Concentration (mol/L)')
%ylabel('Concentration of yeast(g/L)')
legend('CA','CB')
%legend('CC')
hold off
end 

function output=diff_eqs(C)
CA=C(1,1);
CB=C(2,1);
CC=C(3,1);
%Define constants
beta=0.5;
gamma=0.95;
u=0.6;
k=0.000782;%mol glucose/(L*min)
n=0;
%If the number of yeast cells is constant
%diff_eq_A=-k*CA^n;
%diff_eq_B=4*k*CA^n;

%If the number of yeast cells is growing
%46.18g/L is the initial concentration of yeast
%0.186 mol/L is the inital concentration of glucose
turnover=-k*0.186^n/46.18;%[mol Glucose/(g yeast * min)]

diff_eq_A=turnover*CC;
diff_eq_B=-4*turnover*CC;
%Concentration of ethanol(CB) must be the same unit of CC: g/L
%Molar mass of ethanol is 46.07 g/mol
diff_eq_C=beta*CC/(1+gamma*CB*46.07+u*CC);

output=[diff_eq_A;diff_eq_B;diff_eq_C];
end

