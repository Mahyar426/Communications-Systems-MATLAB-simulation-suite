%% Mahyar Onsori 309823
%%
clc;
clear all;
close all;
tic;
%% Generating g(D)
% We choose a primitive polynomial of degree 7 and multiply it by (D+1)
% p(D)= D^7+D^3+1
% g(D)=p(D).(D+1)= D^8+D^7+D^4+D^3+D+1
crc_encoder = comm.CRCGenerator('Polynomial','z^8 + z^7 + z^4 + z^3 + z + 1');
%% Choosing k=8 , Encoding with g(D)
k=8;
r=8;
n=r+k;
v=zeros(2^k,k);
for i=1:(2^k)-1
    v(i+1,:)=de2bi(i,k);
    c(i+1,:)=(crc_encoder(v(i+1,:)'))';
end
%% Computing the Distance Profile
for i=1:size(c,1)
    Wh(i,:)=c(1,:)+c(i,:);
end
Wh=mod(Wh,2);
dist_profile=sum(Wh,2);
dist_min=min(dist_profile((2:end),:));
%% Computing Ai
Ai=zeros(n+1,1);
for i=1:size(dist_profile,1)
    Ai(dist_profile(i,1)+1,1)=Ai(dist_profile(i,1)+1,1)+1;
end
%% Computing the Undetected Error Probability for Analytical
powers_analytical=[0.4,0.3,0.2,0.1,0.05,0.01,0.001,0.0001,0.00001];
P_UE_analytical=zeros(1,size(powers_analytical,2));
for Power=1:size(powers_analytical,2)
    p=powers_analytical(1,Power);
    p_not=1-p;
    for i=2:n
        P_UE_analytical(1,Power)=P_UE_analytical(1,Power)+((Ai(i+1,1))*(p^i)*(p_not^(n-i)));
    end
end
%% Computing the Undetected Error Probability for Simulation
crc_detector = comm.CRCDetector('Polynomial','z^8 + z^7 + z^4 + z^3 + z + 1');
powers_Simulation=[0.4,0.3,0.2,0.1];
Max_errors=100;
Num_errors_1=0;
count_error_1=0;
Num_errors_2=0;
count_error_2=0;
Num_errors_3=0;
count_error_3=0;
Num_errors_4=0;
count_error_4=0;
% p=0.4
while(Num_errors_1<Max_errors)
    Rand1=randi(size(c,1),1);
    C_rand_1=c(Rand1,:);
    y1=bsc(C_rand_1,0.4);
    [~,syndrome1]=(crc_detector(y1'));
    if syndrome1==0
        if ~isequal(y1,C_rand_1)
            Num_errors_1=Num_errors_1+1;
        end
    end
    count_error_1=count_error_1+1;
end
P_UE_Simulation(1,1)=Num_errors_1/count_error_1;
% p=0.3
while(Num_errors_2<Max_errors)
    Rand2=randi(size(c,1),1);
    C_rand_2=c(Rand2,:);
    y2=bsc(C_rand_2,0.3);
    [~,syndrome2]=(crc_detector(y2'));
    if syndrome2==0
        if ~isequal(y2,C_rand_2)
            Num_errors_2=Num_errors_2+1;
        end
    end
    count_error_2=count_error_2+1;
end
P_UE_Simulation(1,2)=Num_errors_2/count_error_2;
% p=0.2
while(Num_errors_3<Max_errors)
    Rand3=randi(size(c,1),1);
    C_rand_3=c(Rand3,:);
    y3=bsc(C_rand_3,0.2);
    [~,syndrome3]=(crc_detector(y3'));
    if syndrome3==0
        if ~isequal(y3,C_rand_3)
            Num_errors_3=Num_errors_3+1;
        end
    end
    count_error_3=count_error_3+1;
end
P_UE_Simulation(1,3)=Num_errors_3/count_error_3;
% p=0.1
while(Num_errors_4<Max_errors)
    Rand4=randi(size(c,1),1);
    C_rand_4=c(Rand4,:);
    y4=bsc(C_rand_4,0.1);
    [~,syndrome4]=(crc_detector(y4'));
    if syndrome4==0
        if ~isequal(y4,C_rand_4)
            Num_errors_4=Num_errors_4+1;
        end
    end
    count_error_4=count_error_4+1;
end
P_UE_Simulation(1,4)=Num_errors_4/count_error_4;
%% Plotting
figure,loglog(powers_analytical,P_UE_analytical,'marker','o','color','b');
grid on
hold on
loglog(powers_Simulation,P_UE_Simulation,'marker','x','color','r');
set(gca,'xdir','reverse')
h1=legend('$$Analytic$$','$$Simulation$$');
set(h1,'Interpreter','Latex','FontSize',13)
h2=xlabel('$$ p $$');
set(h2,'Interpreter','Latex','FontSize',13)
h3=ylabel('$$ Pr(UE) $$');
set(h3,'Interpreter','Latex','FontSize',13)

toc;