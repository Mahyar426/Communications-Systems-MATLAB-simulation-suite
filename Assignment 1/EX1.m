%% Mahyar Onsori 309823
%%
clc;
clear all;
close all;
%% Parameters
tic;
n=8;
k=4;
r=n-k;
%% Defining Generator Matrices
G1=[1,0,0,0,0,1,1,1
 0,1,0,0,1,0,1,1
 0,0,1,0,1,1,0,1
 0,0,0,1,1,1,1,0];

G2=[1,0,0,0,0,1,1,0
 0,1,0,0,1,0,1,1
 0,0,1,0,1,1,0,1
 0,0,0,1,1,1,1,0];

%% Generating Non Zero Binary Vectors
Number_of_Vectors=(2^k);
v=zeros(2^k,k);
for i=1:Number_of_Vectors-1
    v(i+1,:)=de2bi(i,k);

end
%% Generating Codeword c1=vG1
for i=1:Number_of_Vectors
    c1(i,:)=v(i,:)*G1;
end
c1=mod(c1,2);
%% Computing Hamming Weight Of Each Codeword C1 and Distance Profile of C2
Wh_1=sum(c1,2);
min_dist_1=min(Wh_1(2:end));
%% Computing Ai_1
Ai_1=zeros(n+1,1);
for i=1:size(Wh_1,1)
    Ai_1(Wh_1(i,1)+1,1)=Ai_1(Wh_1(i,1)+1,1)+1;
end
%% Computing the Undetected Error Probability for First Code
P_UE_1=zeros(1,6);
for Power=1:6
    p=double(10^(-1*Power));
    p_not=double(1-p);
    for i=2:n
        P_UE_1(1,Power)=P_UE_1(1,Power)+((Ai_1(i,1))*(p^(i-1))*(p_not^(n-i+1)));
    end
end
%% Generating Codeword c2=vG2
for i=1:Number_of_Vectors
    c2(i,:)=v(i,:)*G2;
end
c2=mod(c2,2);
%% Computing Hamming Weight Of Each Codeword C2 and Distance Profile of C2
Wh_2=sum(c2,2);
min_dist_2=min(Wh_2(2:end));
%% Computing Ai_2
Ai_2=zeros(n+1,1);
for i=1:size(Wh_2,1)
    Ai_2(Wh_2(i,1)+1,1)=Ai_2(Wh_2(i,1)+1,1)+1;
end
%% Computing the Undetected Error Probability for Second Code
P_UE_2=zeros(1,6);
for Power=1:6
    p=double(10^(-1*Power));
    p_not=double(1-p);
    for i=2:n
        P_UE_2(1,Power)=P_UE_2(1,Power)+((Ai_2(i,1))*(p^(i-1))*(p_not^(n-i+1)));
    end
end
%% Building Parity Matrices H1,H2
H1=zeros(n,r);
H1(1:k,1:k)=G1(1:k,k+1:n);
H1(k+1:n,1:k)=eye(r);
H2=zeros(n,r);
H2(1:k,1:k)=G2(1:k,k+1:n);
H2(k+1:n,1:k)=eye(r);
%% Simulation Part for G1
Max_errors=100;
Num_errors_1_1=0;
count_error_1_1=0;
Num_errors_1_2=0;
count_error_1_2=0;
Num_errors_1_3=0;
count_error_1_3=0;
Num_errors_1_4=0;
count_error_1_4=0;
Num_errors_1_5=0;
count_error_1_5=0;
% p=0.1
while(Num_errors_1_1<Max_errors)
v_rand_1_1=(randi(2,1,k))-1;
C_rand_1_1=v_rand_1_1*G1;
C_rand_1_1=mod(C_rand_1_1,2);
y1_1=bsc(C_rand_1_1,0.1);
syndrome1_1=y1_1*H1;
syndrome1_1=mod(syndrome1_1,2);
if syndrome1_1==0
    if ~isequal(y1_1,C_rand_1_1)
        Num_errors_1_1=Num_errors_1_1+1; 
    end
end
count_error_1_1=count_error_1_1+1;
end
Ratio(1,1)=Num_errors_1_1/count_error_1_1;
% p=0.09
while(Num_errors_1_2<Max_errors)
v_rand_1_2=(randi(2,1,k))-1;
C_rand_1_2=v_rand_1_2*G1;
C_rand_1_2=mod(C_rand_1_2,2);
y1_2=bsc(C_rand_1_2,0.09);
syndrome1_2=y1_2*H1;
syndrome1_2=mod(syndrome1_2,2);
if syndrome1_2==0
    if ~isequal(y1_2,C_rand_1_2)
        Num_errors_1_2=Num_errors_1_2+1; 
    end
end
count_error_1_2=count_error_1_2+1;
end
Ratio(1,2)=Num_errors_1_2/count_error_1_2;
% p=0.08
while(Num_errors_1_3<Max_errors)
v_rand_1_3=(randi(2,1,k))-1;
C_rand_1_3=v_rand_1_3*G1;
C_rand_1_3=mod(C_rand_1_3,2);
y1_3=bsc(C_rand_1_3,0.08);
syndrome1_3=y1_3*H1;
syndrome1_3=mod(syndrome1_3,2);
if syndrome1_3==0
    if ~isequal(y1_3,C_rand_1_3)
        Num_errors_1_3=Num_errors_1_3+1; 
    end
end
count_error_1_3=count_error_1_3+1;
end
Ratio(1,3)=Num_errors_1_3/count_error_1_3;
% p=0.07
while(Num_errors_1_4<Max_errors)
v_rand_1_4=(randi(2,1,k))-1;
C_rand_1_4=v_rand_1_4*G1;
C_rand_1_4=mod(C_rand_1_4,2);
y1_4=bsc(C_rand_1_4,0.07);
syndrome1_4=y1_4*H1;
syndrome1_4=mod(syndrome1_4,2);
if syndrome1_4==0
    if ~isequal(y1_4,C_rand_1_4)
        Num_errors_1_4=Num_errors_1_4+1; 
    end
end
count_error_1_4=count_error_1_4+1;
end
Ratio(1,4)=Num_errors_1_4/count_error_1_4;
% p=0.06
while(Num_errors_1_5<Max_errors)
v_rand_1_5=(randi(2,1,k))-1;
C_rand_1_5=v_rand_1_5*G1;
C_rand_1_5=mod(C_rand_1_5,2);
y1_5=bsc(C_rand_1_5,0.06);
syndrome1_5=y1_5*H1;
syndrome1_5=mod(syndrome1_5,2);
if syndrome1_5==0
    if ~isequal(y1_5,C_rand_1_5)
        Num_errors_1_5=Num_errors_1_5+1; 
    end
end
count_error_1_5=count_error_1_5+1;
end
Ratio(1,5)=Num_errors_1_5/count_error_1_5;
%% Simulation Part for G2
Num_errors_2_1=0;
count_error_2_1=0;
Num_errors_2_2=0;
count_error_2_2=0;
Num_errors_2_3=0;
count_error_2_3=0;
Num_errors_2_4=0;
count_error_2_4=0;
Num_errors_2_5=0;
count_error_2_5=0;
% p=0.1
while(Num_errors_2_1<Max_errors)
v_rand_2_1=(randi(2,1,k))-1;
C_rand_2_1=v_rand_2_1*G2;
C_rand_2_1=mod(C_rand_2_1,2);
y2_1=bsc(C_rand_2_1,0.1);
syndrome2_1=y2_1*H2;
syndrome2_1=mod(syndrome2_1,2);
if syndrome2_1==0
    if ~isequal(y2_1,C_rand_2_1)
        Num_errors_2_1=Num_errors_2_1+1; 
    end
end
count_error_2_1=count_error_2_1+1;
end
Ratio2(1,1)=Num_errors_2_1/count_error_2_1;
% p=0.09
while(Num_errors_2_2<Max_errors)
v_rand_2_2=(randi(2,1,k))-1;
C_rand_2_2=v_rand_2_2*G2;
C_rand_2_2=mod(C_rand_2_2,2);
y2_2=bsc(C_rand_2_2,0.09);
syndrome2_2=y2_2*H2;
syndrome2_2=mod(syndrome2_2,2);
if syndrome2_2==0
    if ~isequal(y2_2,C_rand_2_2)
        Num_errors_2_2=Num_errors_2_2+1; 
    end
end
count_error_2_2=count_error_2_2+1;
end
Ratio2(1,2)=Num_errors_2_2/count_error_2_2;
% p=0.08
while(Num_errors_2_3<Max_errors)
v_rand_2_3=(randi(2,1,k))-1;
C_rand_2_3=v_rand_2_3*G2;
C_rand_2_3=mod(C_rand_2_3,2);
y2_3=bsc(C_rand_2_3,0.08);
syndrome2_3=y2_3*H2;
syndrome2_3=mod(syndrome2_3,2);
if syndrome2_3==0
    if ~isequal(y2_3,C_rand_2_3)
        Num_errors_2_3=Num_errors_2_3+1; 
    end
end
count_error_2_3=count_error_2_3+1;
end
Ratio2(1,3)=Num_errors_2_3/count_error_2_3;
% p=0.07
while(Num_errors_2_4<Max_errors)
v_rand_2_4=(randi(2,1,k))-1;
C_rand_2_4=v_rand_2_4*G2;
C_rand_2_4=mod(C_rand_2_4,2);
y2_4=bsc(C_rand_2_4,0.07);
syndrome2_4=y2_4*H2;
syndrome2_4=mod(syndrome2_4,2);
if syndrome2_4==0
    if ~isequal(y2_4,C_rand_2_4)
        Num_errors_2_4=Num_errors_2_4+1; 
    end
end
count_error_2_4=count_error_2_4+1;
end
Ratio2(1,4)=Num_errors_2_4/count_error_2_4;
% p=0.06
while(Num_errors_2_5<Max_errors)
v_rand_2_5=(randi(2,1,k))-1;
C_rand_2_5=v_rand_2_5*G2;
C_rand_2_5=mod(C_rand_2_5,2);
y2_5=bsc(C_rand_2_5,0.06);
syndrome2_5=y2_5*H2;
syndrome2_5=mod(syndrome2_5,2);
if syndrome2_5==0
    if ~isequal(y2_5,C_rand_2_5)
        Num_errors_2_5=Num_errors_2_5+1; 
    end
end
count_error_2_5=count_error_2_5+1;
end
Ratio2(1,5)=Num_errors_2_5/count_error_2_5;
%% Plotting Part
Powers_Simulation=[0.1,0.09,0.08,0.07,0.06];
xplot=logspace(-1,-6,6);
figure,loglog(xplot,P_UE_1,'marker','o');
grid on
hold on
loglog(xplot,P_UE_2,'marker','o','color','r')
hold on
loglog(Powers_Simulation,Ratio,'marker','square','color','k')
hold on
loglog(Powers_Simulation,Ratio2,'marker','square','color','g')
set(gca,'xdir','reverse');
h1=legend('Code 1 Analytic','Code 2 Analytic','Code 1 Simulation','Code 2 Simulation');
set(h1,'Interpreter','Latex','FontSize',13)
h2=xlabel('$$ p $$');
set(h2,'Interpreter','Latex','FontSize',13)
h3=ylabel('$$ Pr(UE) $$');
set(h3,'Interpreter','Latex','FontSize',13)
toc;