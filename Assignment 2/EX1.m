%% Mahyar Onsori 309823
clc
clear all
close all
%% Parameters and Generator Matrix
tic;
n=8;
k=4;
r=n-k;
G=[1,0,0,0,0,1,1,0
 0,1,0,0,1,0,1,1
 0,0,1,0,1,1,0,1
 0,0,0,1,1,1,1,0];
H=zeros(n,r);
H(1:k,1:k)=G(1:k,k+1:n);
H(k+1:n,1:k)=eye(r);
%% Generating Binary Information Vectors
Number_of_Vectors=(2^k);
v=zeros(2^k,k);
for i=1:Number_of_Vectors-1
    v(i+1,:)=de2bi(i,k);
end
%% Generating codebook c=vG
for i=1:Number_of_Vectors
    c(i,:)=v(i,:)*G;
end
c=mod(c,2);
Wh=sum(c,2);
min_dist=min(Wh(2:end));
t=floor((min_dist-1)/2);
Ai=zeros(n+1,1);
for i=1:size(Wh,1)
    Ai(Wh(i,1)+1,1)=Ai(Wh(i,1)+1,1)+1;
end
Aimin=Ai(min_dist+1);
for i=1:size(c,1)
    for j=1:size(c,2)
        if c(i,j)==0
            cprime(i,j)=-1;
        else
            cprime(i,j)=1;
        end
    end
end
%% Computing the error rate for hard correction using analytical formula
for ebnodb=1:11
    ebno(1,ebnodb)=10^((ebnodb-1)/10);
    p=0.5*erfc(sqrt((k/n)*ebno(1,ebnodb)));
    sum_rate_hard_analytical(1,ebnodb)=0;
    for i=1:t+1
        sum_rate_hard_analytical(1,ebnodb)=sum_rate_hard_analytical(1,ebnodb)+nchoosek(n,i-1)*(p^(i-1))*((1-p)^(n-i+1));
        rate_hard_analytical(1,ebnodb)=1-sum_rate_hard_analytical(1,ebnodb);
    end
end
ebnodb=0:10;
semilogy(ebnodb,rate_hard_analytical,'LineWidth',1);
h2=xlabel('$$ E_b/N_0 $$');
set(h2,'Interpreter','Latex','FontSize',13);
h3=ylabel('$$ P_w (e) $$');
set(h3,'Interpreter','Latex','FontSize',13);
h4=title('Codeword Error Rate');
set(h4,'Interpreter','Latex','FontSize',13);

grid on
hold on
%% Soft Correction Union Bound and Asymptotic Expression
P_W_E_Soft_Asymp=zeros(1,11);
for i=1:11
    P_W_E_Soft_Asymp(1,i)=1/2*Aimin*erfc(sqrt((k/n)*min_dist*ebno(1,i)));
end
semilogy(ebnodb,P_W_E_Soft_Asymp,'LineWidth',1,'LineStyle','--');
P_W_E_Soft_Union=zeros(1,11);
for i=1:11
    for j=2:n+1
        P_W_E_Soft_Union(1,i)=P_W_E_Soft_Union(1,i)+1/2*Ai(j,1)*erfc(sqrt((k/n)*(j-1)*ebno(1,i)));
    end
end
semilogy(ebnodb,P_W_E_Soft_Union,'LineWidth',1);
%% Building possible error vectors and LUT
num_error_corrected=0;
for i=1:t+1
    num_error_corrected=num_error_corrected+nchoosek(n,(i-1));
end
e_lut=[zeros(1,num_error_corrected-1);eye(num_error_corrected-1)];
syndrome_lut=[zeros(1,size(H,2));H];
%% Computing the error rate for hard correction using simulation
for i=1:8
    Sigma2(1,i)=1/(2*(k/n)*(ebno(1,i)));
end
Sigma=sqrt(Sigma2);
Max_errors=500;
Num_errors_0_H=0;count_error_0_H=0;Num_errors_1_H=0;count_error_1_H=0;Num_errors_2_H=0;count_error_2_H=0;Num_errors_3_H=0;count_error_3_H=0;
Num_errors_4_H=0;count_error_4_H=0;Num_errors_5_H=0;count_error_5_H=0;Num_errors_6_H=0;count_error_6_H=0;Num_errors_7_H=0;count_error_7_H=0;
% eb/no=0 db
while(Num_errors_0_H<Max_errors)
    v_rand_0_H=(randi(2,1,k))-1;
    C_rand_0_H=v_rand_0_H*G;
    C_rand_0_H=mod(C_rand_0_H,2);
    for i=1:size(C_rand_0_H,2)
        if C_rand_0_H(1,i)==1
            C_Transmit_0_H(1,i)=1;
        else
            C_Transmit_0_H(1,i)=-1;
        end
    end
    Noise_0_H=randn(1,n)*Sigma(1,1);
    r_0_H=C_Transmit_0_H+Noise_0_H;
    for j=1:size(r_0_H,2)
        if r_0_H(1,j)>0
            y_0_H(1,j)=1;
        else
            y_0_H(1,j)=0;
        end
    end
    Syndrome_0_H=y_0_H*H;
    Syndrome_0_H=mod(Syndrome_0_H,2);
    flag_0_H=0;
    for ii=1:size(syndrome_lut,1)
        if isequal(Syndrome_0_H,syndrome_lut(ii,:))
            C_corrected_0_H=y_0_H+e_lut(ii,:);
            C_corrected_0_H=mod(C_corrected_0_H,2);
            if ~isequal(C_rand_0_H,C_corrected_0_H)
                Num_errors_0_H=Num_errors_0_H+1;
            end
        end
        if ~isequal(Syndrome_0_H,syndrome_lut(ii,:))
            flag_0_H=flag_0_H+1;
        end
        if flag_0_H==size(syndrome_lut,1)
               Num_errors_0_H=Num_errors_0_H+1;
        end
    end
    count_error_0_H=count_error_0_H+1;
end
rate_hard_simulation(1,1)=Num_errors_0_H/count_error_0_H;

% eb/no=1 db
while(Num_errors_1_H<Max_errors)
    v_rand_1_H=(randi(2,1,k))-1;
    C_rand_1_H=v_rand_1_H*G;
    C_rand_1_H=mod(C_rand_1_H,2);
    for i=1:size(C_rand_1_H,2)
        if C_rand_1_H(1,i)==1
            C_Transmit_1_H(1,i)=1;
        else
            C_Transmit_1_H(1,i)=-1;
        end
    end
    Noise_1_H=randn(1,n)*Sigma(1,2);
    r_1_H=C_Transmit_1_H+Noise_1_H;
    for j=1:size(r_1_H,2)
        if r_1_H(1,j)>0
            y_1_H(1,j)=1;
        else
            y_1_H(1,j)=0;
        end
    end
    Syndrome_1_H=y_1_H*H;
    Syndrome_1_H=mod(Syndrome_1_H,2);
    flag_1_H=0;
    for ii=1:size(syndrome_lut,1)
        if isequal(Syndrome_1_H,syndrome_lut(ii,:))
            C_corrected_1_H=y_1_H+e_lut(ii,:);
            C_corrected_1_H=mod(C_corrected_1_H,2);
            if ~isequal(C_rand_1_H,C_corrected_1_H)
                Num_errors_1_H=Num_errors_1_H+1;
            end
        end
        if ~isequal(Syndrome_1_H,syndrome_lut(ii,:))
            flag_1_H=flag_1_H+1;
        end
        if flag_1_H==size(syndrome_lut,1)
               Num_errors_1_H=Num_errors_1_H+1;
        end
    end
    count_error_1_H=count_error_1_H+1;
end
rate_hard_simulation(1,2)=Num_errors_1_H/count_error_1_H;

% eb/no=2 db
while(Num_errors_2_H<Max_errors)
    v_rand_2_H=(randi(2,1,k))-1;
    C_rand_2_H=v_rand_2_H*G;
    C_rand_2_H=mod(C_rand_2_H,2);
    for i=1:size(C_rand_2_H,2)
        if C_rand_2_H(1,i)==1
            C_Transmit_2_H(1,i)=1;
        else
            C_Transmit_2_H(1,i)=-1;
        end
    end
    Noise_2_H=randn(1,n)*Sigma(1,3);
    r_2_H=C_Transmit_2_H+Noise_2_H;
    for j=1:size(r_2_H,2)
        if r_2_H(1,j)>0
            y_2_H(1,j)=1;
        else
            y_2_H(1,j)=0;
        end
    end
    Syndrome_2_H=y_2_H*H;
    Syndrome_2_H=mod(Syndrome_2_H,2);
    flag_2_H=0;
    for ii=1:size(syndrome_lut,1)
        if isequal(Syndrome_2_H,syndrome_lut(ii,:))
            C_corrected_2_H=y_2_H+e_lut(ii,:);
            C_corrected_2_H=mod(C_corrected_2_H,2);
            if ~isequal(C_rand_2_H,C_corrected_2_H)
                Num_errors_2_H=Num_errors_2_H+1;
            end
        end
        if ~isequal(Syndrome_2_H,syndrome_lut(ii,:))
            flag_2_H=flag_2_H+1;
        end
        if flag_2_H==size(syndrome_lut,1)
               Num_errors_2_H=Num_errors_2_H+1;
        end
    end
    count_error_2_H=count_error_2_H+1;
end
rate_hard_simulation(1,3)=Num_errors_2_H/count_error_2_H;

% eb/no=3 db
while(Num_errors_3_H<Max_errors)
    v_rand_3_H=(randi(2,1,k))-1;
    C_rand_3_H=v_rand_3_H*G;
    C_rand_3_H=mod(C_rand_3_H,2);
    for i=1:size(C_rand_3_H,2)
        if C_rand_3_H(1,i)==1
            C_Transmit_3_H(1,i)=1;
        else
            C_Transmit_3_H(1,i)=-1;
        end
    end
    Noise_3_H=randn(1,n)*Sigma(1,4);
    r_3_H=C_Transmit_3_H+Noise_3_H;
    for j=1:size(r_3_H,2)
        if r_3_H(1,j)>0
            y_3_H(1,j)=1;
        else
            y_3_H(1,j)=0;
        end
    end
    Syndrome_3_H=y_3_H*H;
    Syndrome_3_H=mod(Syndrome_3_H,2);
    flag_3_H=0;
    for ii=1:size(syndrome_lut,1)
        if isequal(Syndrome_3_H,syndrome_lut(ii,:))
            C_corrected_3_H=y_3_H+e_lut(ii,:);
            C_corrected_3_H=mod(C_corrected_3_H,2);
            if ~isequal(C_rand_3_H,C_corrected_3_H)
                Num_errors_3_H=Num_errors_3_H+1;
            end
        end
        if ~isequal(Syndrome_3_H,syndrome_lut(ii,:))
            flag_3_H=flag_3_H+1;
        end
        if flag_3_H==size(syndrome_lut,1)
               Num_errors_3_H=Num_errors_3_H+1;
        end
    end
    count_error_3_H=count_error_3_H+1;
end
rate_hard_simulation(1,4)=Num_errors_3_H/count_error_3_H;

% eb/no=4 db
while(Num_errors_4_H<Max_errors)
    v_rand_4_H=(randi(2,1,k))-1;
    C_rand_4_H=v_rand_4_H*G;
    C_rand_4_H=mod(C_rand_4_H,2);
    for i=1:size(C_rand_4_H,2)
        if C_rand_4_H(1,i)==1
            C_Transmit_4_H(1,i)=1;
        else
            C_Transmit_4_H(1,i)=-1;
        end
    end
    Noise_4_H=randn(1,n)*Sigma(1,5);
    r_4_H=C_Transmit_4_H+Noise_4_H;
    for j=1:size(r_4_H,2)
        if r_4_H(1,j)>0
            y_4_H(1,j)=1;
        else
            y_4_H(1,j)=0;
        end
    end
    Syndrome_4_H=y_4_H*H;
    Syndrome_4_H=mod(Syndrome_4_H,2);
    flag_4_H=0;
    for ii=1:size(syndrome_lut,1)
        if isequal(Syndrome_4_H,syndrome_lut(ii,:))
            C_corrected_4_H=y_4_H+e_lut(ii,:);
            C_corrected_4_H=mod(C_corrected_4_H,2);
            if ~isequal(C_rand_4_H,C_corrected_4_H)
                Num_errors_4_H=Num_errors_4_H+1;
            end
        end
        if ~isequal(Syndrome_4_H,syndrome_lut(ii,:))
            flag_4_H=flag_4_H+1;
        end
        if flag_4_H==size(syndrome_lut,1)
               Num_errors_4_H=Num_errors_4_H+1;
        end
    end
    count_error_4_H=count_error_4_H+1;
end
rate_hard_simulation(1,5)=Num_errors_4_H/count_error_4_H;

% eb/no=5 db
while(Num_errors_5_H<Max_errors)
    v_rand_5_H=(randi(2,1,k))-1;
    C_rand_5_H=v_rand_5_H*G;
    C_rand_5_H=mod(C_rand_5_H,2);
    for i=1:size(C_rand_5_H,2)
        if C_rand_5_H(1,i)==1
            C_Transmit_5_H(1,i)=1;
        else
            C_Transmit_5_H(1,i)=-1;
        end
    end
    Noise_5_H=randn(1,n)*Sigma(1,6);
    r_5_H=C_Transmit_5_H+Noise_5_H;
    for j=1:size(r_5_H,2)
        if r_5_H(1,j)>0
            y_5_H(1,j)=1;
        else
            y_5_H(1,j)=0;
        end
    end
    Syndrome_5_H=y_5_H*H;
    Syndrome_5_H=mod(Syndrome_5_H,2);
    flag_5_H=0;
    for ii=1:size(syndrome_lut,1)
        if isequal(Syndrome_5_H,syndrome_lut(ii,:))
            C_corrected_5_H=y_5_H+e_lut(ii,:);
            C_corrected_5_H=mod(C_corrected_5_H,2);
            if ~isequal(C_rand_5_H,C_corrected_5_H)
                Num_errors_5_H=Num_errors_5_H+1;
            end
        end
        if ~isequal(Syndrome_5_H,syndrome_lut(ii,:))
            flag_5_H=flag_5_H+1;
        end
        if flag_5_H==size(syndrome_lut,1)
               Num_errors_5_H=Num_errors_5_H+1;
        end
    end
    count_error_5_H=count_error_5_H+1;
end
rate_hard_simulation(1,6)=Num_errors_5_H/count_error_5_H;

% eb/no=6 db
while(Num_errors_6_H<Max_errors)
    v_rand_6_H=(randi(2,1,k))-1;
    C_rand_6_H=v_rand_6_H*G;
    C_rand_6_H=mod(C_rand_6_H,2);
    for i=1:size(C_rand_6_H,2)
        if C_rand_6_H(1,i)==1
            C_Transmit_6_H(1,i)=1;
        else
            C_Transmit_6_H(1,i)=-1;
        end
    end
    Noise_6_H=randn(1,n)*Sigma(1,7);
    r_6_H=C_Transmit_6_H+Noise_6_H;
    for j=1:size(r_6_H,2)
        if r_6_H(1,j)>0
            y_6_H(1,j)=1;
        else
            y_6_H(1,j)=0;
        end
    end
    Syndrome_6_H=y_6_H*H;
    Syndrome_6_H=mod(Syndrome_6_H,2);
    flag_6_H=0;
    for ii=1:size(syndrome_lut,1)
        if isequal(Syndrome_6_H,syndrome_lut(ii,:))
            C_corrected_6_H=y_6_H+e_lut(ii,:);
            C_corrected_6_H=mod(C_corrected_6_H,2);
            if ~isequal(C_rand_6_H,C_corrected_6_H)
                Num_errors_6_H=Num_errors_6_H+1;
            end
        end
        if ~isequal(Syndrome_6_H,syndrome_lut(ii,:))
            flag_6_H=flag_6_H+1;
        end
        if flag_6_H==size(syndrome_lut,1)
               Num_errors_6_H=Num_errors_6_H+1;
        end
    end
    count_error_6_H=count_error_6_H+1;
end
rate_hard_simulation(1,7)=Num_errors_6_H/count_error_6_H;

% eb/no=7 db
while(Num_errors_7_H<Max_errors)
    v_rand_7_H=(randi(2,1,k))-1;
    C_rand_7_H=v_rand_7_H*G;
    C_rand_7_H=mod(C_rand_7_H,2);
    for i=1:size(C_rand_7_H,2)
        if C_rand_7_H(1,i)==1
            C_Transmit_7_H(1,i)=1;
        else
            C_Transmit_7_H(1,i)=-1;
        end
    end
    Noise_7_H=randn(1,n)*Sigma(1,8);
    r_7_H=C_Transmit_7_H+Noise_7_H;
    for j=1:size(r_7_H,2)
        if r_7_H(1,j)>0
            y_7_H(1,j)=1;
        else
            y_7_H(1,j)=0;
        end
    end
    Syndrome_7_H=y_7_H*H;
    Syndrome_7_H=mod(Syndrome_7_H,2);
    flag_7_H=0;
    for ii=1:size(syndrome_lut,1)
        if isequal(Syndrome_7_H,syndrome_lut(ii,:))
            C_corrected_7_H=y_7_H+e_lut(ii,:);
            C_corrected_7_H=mod(C_corrected_7_H,2);
            if ~isequal(C_rand_7_H,C_corrected_7_H)
                Num_errors_7_H=Num_errors_7_H+1;
            end
        end
        if ~isequal(Syndrome_7_H,syndrome_lut(ii,:))
            flag_7_H=flag_7_H+1;
        end
        if flag_7_H==size(syndrome_lut,1)
               Num_errors_7_H=Num_errors_7_H+1;
        end
    end
    count_error_7_H=count_error_7_H+1;
end
rate_hard_simulation(1,8)=Num_errors_7_H/count_error_7_H;

LOGS=0:7;
semilogy(LOGS,rate_hard_simulation,'LineWidth',1,'Marker','o','Color','m');
%% Computing the error rate for soft correction using simulation
Num_errors_0_S=0;count_error_0_S=0;Num_errors_1_S=0;count_error_1_S=0;Num_errors_2_S=0;count_error_2_S=0;Num_errors_3_S=0;count_error_3_S=0;
Num_errors_4_S=0;count_error_4_S=0;Num_errors_5_S=0;count_error_5_S=0;Num_errors_6_S=0;count_error_6_S=0;Num_errors_7_S=0;count_error_7_S=0;
% Eb/N0=0db
while(Num_errors_0_S<Max_errors)
    v_rand_0_S=(randi(2,1,k))-1;
    C_rand_0_S=v_rand_0_S*G;
    C_rand_0_S=mod(C_rand_0_S,2);
    for i=1:size(C_rand_0_S,2)
        if C_rand_0_S(1,i)==1
            C_Transmit_0_S(1,i)=1;
        else
            C_Transmit_0_S(1,i)=-1;
        end
    end
    Noise_0_S=randn(1,n)*Sigma(1,1);
    r_0_S=C_Transmit_0_S+Noise_0_S;
    for i=1:size(cprime,1)
        Euclidean_0(i,:)=sqrt(sum((r_0_S - cprime(i,:)) .^ 2));
    end
    [~,Index_Min_Euclidean_0]=min(Euclidean_0);
    C_Corrected_0_S=c(Index_Min_Euclidean_0,:);
    if ~isequal(C_Corrected_0_S,C_rand_0_S)
        Num_errors_0_S=Num_errors_0_S+1;
    end
    count_error_0_S=count_error_0_S+1;
end
rate_soft_simulation(1,1)=Num_errors_0_S/count_error_0_S;

% Eb/N0=1db
while(Num_errors_1_S<Max_errors)
    v_rand_1_S=(randi(2,1,k))-1;
    C_rand_1_S=v_rand_1_S*G;
    C_rand_1_S=mod(C_rand_1_S,2);
    for i=1:size(C_rand_1_S,2)
        if C_rand_1_S(1,i)==1
            C_Transmit_1_S(1,i)=1;
        else
            C_Transmit_1_S(1,i)=-1;
        end
    end
    Noise_1_S=randn(1,n)*Sigma(1,2);
    r_1_S=C_Transmit_1_S+Noise_1_S;
    for i=1:size(cprime,1)
        Euclidean_1(i,:)=sqrt(sum((r_1_S - cprime(i,:)) .^ 2));
    end
    [~,Index_Min_Euclidean_1]=min(Euclidean_1);
    C_Corrected_1_S=c(Index_Min_Euclidean_1,:);
    if ~isequal(C_Corrected_1_S,C_rand_1_S)
        Num_errors_1_S=Num_errors_1_S+1;
    end
    count_error_1_S=count_error_1_S+1;
end
rate_soft_simulation(1,2)=Num_errors_1_S/count_error_1_S;

% Eb/N0=2db
while(Num_errors_2_S<Max_errors)
    v_rand_2_S=(randi(2,1,k))-1;
    C_rand_2_S=v_rand_2_S*G;
    C_rand_2_S=mod(C_rand_2_S,2);
    for i=1:size(C_rand_2_S,2)
        if C_rand_2_S(1,i)==1
            C_Transmit_2_S(1,i)=1;
        else
            C_Transmit_2_S(1,i)=-1;
        end
    end
    Noise_2_S=randn(1,n)*Sigma(1,3);
    r_2_S=C_Transmit_2_S+Noise_2_S;
    for i=1:size(cprime,1)
        Euclidean_2(i,:)=sqrt(sum((r_2_S - cprime(i,:)) .^ 2));
    end
    [~,Index_Min_Euclidean_2]=min(Euclidean_2);
    C_Corrected_2_S=c(Index_Min_Euclidean_2,:);
    if ~isequal(C_Corrected_2_S,C_rand_2_S)
        Num_errors_2_S=Num_errors_2_S+1;
    end
    count_error_2_S=count_error_2_S+1;
end
rate_soft_simulation(1,3)=Num_errors_2_S/count_error_2_S;

% Eb/N0=3db
while(Num_errors_3_S<Max_errors)
    v_rand_3_S=(randi(2,1,k))-1;
    C_rand_3_S=v_rand_3_S*G;
    C_rand_3_S=mod(C_rand_3_S,2);
    for i=1:size(C_rand_3_S,2)
        if C_rand_3_S(1,i)==1
            C_Transmit_3_S(1,i)=1;
        else
            C_Transmit_3_S(1,i)=-1;
        end
    end
    Noise_3_S=randn(1,n)*Sigma(1,4);
    r_3_S=C_Transmit_3_S+Noise_3_S;
    for i=1:size(cprime,1)
        Euclidean_3(i,:)=sqrt(sum((r_3_S - cprime(i,:)) .^ 2));
    end
    [~,Index_Min_Euclidean_3]=min(Euclidean_3);
    C_Corrected_3_S=c(Index_Min_Euclidean_3,:);
    if ~isequal(C_Corrected_3_S,C_rand_3_S)
        Num_errors_3_S=Num_errors_3_S+1;
    end
    count_error_3_S=count_error_3_S+1;
end
rate_soft_simulation(1,4)=Num_errors_3_S/count_error_3_S;

% Eb/N0=4db
while(Num_errors_4_S<Max_errors)
    v_rand_4_S=(randi(2,1,k))-1;
    C_rand_4_S=v_rand_4_S*G;
    C_rand_4_S=mod(C_rand_4_S,2);
    for i=1:size(C_rand_4_S,2)
        if C_rand_4_S(1,i)==1
            C_Transmit_4_S(1,i)=1;
        else
            C_Transmit_4_S(1,i)=-1;
        end
    end
    Noise_4_S=randn(1,n)*Sigma(1,5);
    r_4_S=C_Transmit_4_S+Noise_4_S;
    for i=1:size(cprime,1)
        Euclidean_4(i,:)=sqrt(sum((r_4_S - cprime(i,:)) .^ 2));
    end
    [~,Index_Min_Euclidean_4]=min(Euclidean_4);
    C_Corrected_4_S=c(Index_Min_Euclidean_4,:);
    if ~isequal(C_Corrected_4_S,C_rand_4_S)
        Num_errors_4_S=Num_errors_4_S+1;
    end
    count_error_4_S=count_error_4_S+1;
end
rate_soft_simulation(1,5)=Num_errors_4_S/count_error_4_S;

% Eb/N0=5db
while(Num_errors_5_S<Max_errors)
    v_rand_5_S=(randi(2,1,k))-1;
    C_rand_5_S=v_rand_5_S*G;
    C_rand_5_S=mod(C_rand_5_S,2);
    for i=1:size(C_rand_5_S,2)
        if C_rand_5_S(1,i)==1
            C_Transmit_5_S(1,i)=1;
        else
            C_Transmit_5_S(1,i)=-1;
        end
    end
    Noise_5_S=randn(1,n)*Sigma(1,6);
    r_5_S=C_Transmit_5_S+Noise_5_S;
    for i=1:size(cprime,1)
        Euclidean_5(i,:)=sqrt(sum((r_5_S - cprime(i,:)) .^ 2));
    end
    [~,Index_Min_Euclidean_5]=min(Euclidean_5);
    C_Corrected_5_S=c(Index_Min_Euclidean_5,:);
    if ~isequal(C_Corrected_5_S,C_rand_5_S)
        Num_errors_5_S=Num_errors_5_S+1;
    end
    count_error_5_S=count_error_5_S+1;
end
rate_soft_simulation(1,6)=Num_errors_5_S/count_error_5_S;

% Eb/N0=6db
while(Num_errors_6_S<Max_errors)
    v_rand_6_S=(randi(2,1,k))-1;
    C_rand_6_S=v_rand_6_S*G;
    C_rand_6_S=mod(C_rand_6_S,2);
    for i=1:size(C_rand_6_S,2)
        if C_rand_6_S(1,i)==1
            C_Transmit_6_S(1,i)=1;
        else
            C_Transmit_6_S(1,i)=-1;
        end
    end
    Noise_6_S=randn(1,n)*Sigma(1,7);
    r_6_S=C_Transmit_6_S+Noise_6_S;
    for i=1:size(cprime,1)
        Euclidean_6(i,:)=sqrt(sum((r_6_S - cprime(i,:)) .^ 2));
    end
    [~,Index_Min_Euclidean_6]=min(Euclidean_6);
    C_Corrected_6_S=c(Index_Min_Euclidean_6,:);
    if ~isequal(C_Corrected_6_S,C_rand_6_S)
        Num_errors_6_S=Num_errors_6_S+1;
    end
    count_error_6_S=count_error_6_S+1;
end
rate_soft_simulation(1,7)=Num_errors_6_S/count_error_6_S;

% Eb/N0=7db
while(Num_errors_7_S<Max_errors)
    v_rand_7_S=(randi(2,1,k))-1;
    C_rand_7_S=v_rand_7_S*G;
    C_rand_7_S=mod(C_rand_7_S,2);
    for i=1:size(C_rand_7_S,2)
        if C_rand_7_S(1,i)==1
            C_Transmit_7_S(1,i)=1;
        else
            C_Transmit_7_S(1,i)=-1;
        end
    end
    Noise_7_S=randn(1,n)*Sigma(1,8);
    r_7_S=C_Transmit_7_S+Noise_7_S;
    for i=1:size(cprime,1)
        Euclidean_7(i,:)=sqrt(sum((r_7_S - cprime(i,:)) .^ 2));
    end
    [~,Index_Min_Euclidean_7]=min(Euclidean_7);
    C_Corrected_7_S=c(Index_Min_Euclidean_7,:);
    if ~isequal(C_Corrected_7_S,C_rand_7_S)
        Num_errors_7_S=Num_errors_7_S+1;
    end
    count_error_7_S=count_error_7_S+1;
end 
rate_soft_simulation(1,8)=Num_errors_7_S/count_error_7_S;
LOGS_H=0:7;
semilogy(LOGS,rate_soft_simulation,'LineWidth',1,'Marker','o')
toc; 
%% Additional
Syndrome_lut_full=zeros(2^r,r);
e_lut_full=zeros(2^r,n);
for i=1:2^r-1
    Syndrome_lut_full(i+1,:)=de2bi(i,r,'left-msb');
    for j=1:size(syndrome_lut,1)
        if isequal(Syndrome_lut_full(i,:),syndrome_lut(j,:))
            e_lut_full(i,:)=e_lut(j,:);
        end
    end
end
weight_2=[1 1 0 0 0 0 0 0];
for i=1:100
    e2=weight_2(randperm(8));
    s2=e2*H;
    s2=mod(s2,2);
    for j=1:size(Syndrome_lut_full,1)
        if isequal(Syndrome_lut_full(j,:),s2)
            if(e_lut_full(j,:)==0)
                e_lut_full(j,:)=e2;
            end
        end
    end
end
% eb/no=0 db
Num_errors_0_A=0;
count_error_0_A=0;
while(Num_errors_0_A<Max_errors)
    v_rand_0_A=(randi(2,1,k))-1;
    C_rand_0_A=v_rand_0_A*G;
    C_rand_0_A=mod(C_rand_0_A,2);
    for i=1:size(C_rand_0_A,2)
        if C_rand_0_A(1,i)==1
            C_Transmit_0_A(1,i)=1;
        else
            C_Transmit_0_A(1,i)=-1;
        end
    end
    Noise_0_A=randn(1,n)*Sigma(1,1);
    r_0_A=C_Transmit_0_A+Noise_0_A;
    for j=1:size(r_0_A,2)
        if r_0_A(1,j)>0
            y_0_A(1,j)=1;
        else
            y_0_A(1,j)=0;
        end
    end
    Syndrome_0_A=y_0_A*H;
    Syndrome_0_A=mod(Syndrome_0_A,2);
    flag_0_A=0;
    for ii=1:size(Syndrome_lut_full,1)
        if isequal(Syndrome_0_A,Syndrome_lut_full(ii,:))
            C_corrected_0_A=y_0_A+e_lut_full(ii,:);
            C_corrected_0_A=mod(C_corrected_0_A,2);
            if ~isequal(C_rand_0_A,C_corrected_0_A)
                Num_errors_0_A=Num_errors_0_A+1;
            end
        end
        if ~isequal(Syndrome_0_A,Syndrome_lut_full(ii,:))
            flag_0_A=flag_0_A+1;
        end
        if flag_0_A==size(Syndrome_lut_full,1)
               Num_errors_0_A=Num_errors_0_A+1;
        end
    end
    count_error_0_A=count_error_0_A+1;
end
rate_hard_simulation_Additional(1,1)=Num_errors_0_A/count_error_0_A;
% eb/no=1 db
Num_errors_1_A=0;
count_error_1_A=0;
while(Num_errors_1_A<Max_errors)
    v_rand_1_A=(randi(2,1,k))-1;
    C_rand_1_A=v_rand_1_A*G;
    C_rand_1_A=mod(C_rand_1_A,2);
    for i=1:size(C_rand_1_A,2)
        if C_rand_1_A(1,i)==1
            C_Transmit_1_A(1,i)=1;
        else
            C_Transmit_1_A(1,i)=-1;
        end
    end
    Noise_1_A=randn(1,n)*Sigma(1,2);
    r_1_A=C_Transmit_1_A+Noise_1_A;
    for j=1:size(r_1_A,2)
        if r_1_A(1,j)>0
            y_1_A(1,j)=1;
        else
            y_1_A(1,j)=0;
        end
    end
    Syndrome_1_A=y_1_A*H;
    Syndrome_1_A=mod(Syndrome_1_A,2);
    flag_1_A=0;
    for ii=1:size(Syndrome_lut_full,1)
        if isequal(Syndrome_1_A,Syndrome_lut_full(ii,:))
            C_corrected_1_A=y_1_A+e_lut_full(ii,:);
            C_corrected_1_A=mod(C_corrected_1_A,2);
            if ~isequal(C_rand_1_A,C_corrected_1_A)
                Num_errors_1_A=Num_errors_1_A+1;
            end
        end
        if ~isequal(Syndrome_1_A,Syndrome_lut_full(ii,:))
            flag_1_A=flag_1_A+1;
        end
        if flag_1_A==size(Syndrome_lut_full,1)
               Num_errors_1_A=Num_errors_1_A+1;
        end
    end
    count_error_1_A=count_error_1_A+1;
end
rate_hard_simulation_Additional(1,2)=Num_errors_1_A/count_error_1_A;
% eb/no=2 db
Num_errors_2_A=0;
count_error_2_A=0;
while(Num_errors_2_A<Max_errors)
    v_rand_2_A=(randi(2,1,k))-1;
    C_rand_2_A=v_rand_2_A*G;
    C_rand_2_A=mod(C_rand_2_A,2);
    for i=1:size(C_rand_2_A,2)
        if C_rand_2_A(1,i)==1
            C_Transmit_2_A(1,i)=1;
        else
            C_Transmit_2_A(1,i)=-1;
        end
    end
    Noise_2_A=randn(1,n)*Sigma(1,3);
    r_2_A=C_Transmit_2_A+Noise_2_A;
    for j=1:size(r_2_A,2)
        if r_2_A(1,j)>0
            y_2_A(1,j)=1;
        else
            y_2_A(1,j)=0;
        end
    end
    Syndrome_2_A=y_2_A*H;
    Syndrome_2_A=mod(Syndrome_2_A,2);
    flag_2_A=0;
    for ii=1:size(Syndrome_lut_full,1)
        if isequal(Syndrome_2_A,Syndrome_lut_full(ii,:))
            C_corrected_2_A=y_2_A+e_lut_full(ii,:);
            C_corrected_2_A=mod(C_corrected_2_A,2);
            if ~isequal(C_rand_2_A,C_corrected_2_A)
                Num_errors_2_A=Num_errors_2_A+1;
            end
        end
        if ~isequal(Syndrome_2_A,Syndrome_lut_full(ii,:))
            flag_2_A=flag_2_A+1;
        end
        if flag_2_A==size(Syndrome_lut_full,1)
               Num_errors_2_A=Num_errors_2_A+1;
        end
    end
    count_error_2_A=count_error_2_A+1;
end
rate_hard_simulation_Additional(1,3)=Num_errors_2_A/count_error_2_A;
% eb/no=3 db
Num_errors_3_A=0;
count_error_3_A=0;
while(Num_errors_3_A<Max_errors)
    v_rand_3_A=(randi(2,1,k))-1;
    C_rand_3_A=v_rand_3_A*G;
    C_rand_3_A=mod(C_rand_3_A,2);
    for i=1:size(C_rand_3_A,2)
        if C_rand_3_A(1,i)==1
            C_Transmit_3_A(1,i)=1;
        else
            C_Transmit_3_A(1,i)=-1;
        end
    end
    Noise_3_A=randn(1,n)*Sigma(1,4);
    r_3_A=C_Transmit_3_A+Noise_3_A;
    for j=1:size(r_3_A,2)
        if r_3_A(1,j)>0
            y_3_A(1,j)=1;
        else
            y_3_A(1,j)=0;
        end
    end
    Syndrome_3_A=y_3_A*H;
    Syndrome_3_A=mod(Syndrome_3_A,2);
    flag_3_A=0;
    for ii=1:size(Syndrome_lut_full,1)
        if isequal(Syndrome_3_A,Syndrome_lut_full(ii,:))
            C_corrected_3_A=y_3_A+e_lut_full(ii,:);
            C_corrected_3_A=mod(C_corrected_3_A,2);
            if ~isequal(C_rand_3_A,C_corrected_3_A)
                Num_errors_3_A=Num_errors_3_A+1;
            end
        end
        if ~isequal(Syndrome_3_A,Syndrome_lut_full(ii,:))
            flag_3_A=flag_3_A+1;
        end
        if flag_3_A==size(Syndrome_lut_full,1)
               Num_errors_3_A=Num_errors_3_A+1;
        end
    end
    count_error_3_A=count_error_3_A+1;
end
rate_hard_simulation_Additional(1,4)=Num_errors_3_A/count_error_3_A;
% eb/no=4 db
Num_errors_4_A=0;
count_error_4_A=0;
while(Num_errors_4_A<Max_errors)
    v_rand_4_A=(randi(2,1,k))-1;
    C_rand_4_A=v_rand_4_A*G;
    C_rand_4_A=mod(C_rand_4_A,2);
    for i=1:size(C_rand_4_A,2)
        if C_rand_4_A(1,i)==1
            C_Transmit_4_A(1,i)=1;
        else
            C_Transmit_4_A(1,i)=-1;
        end
    end
    Noise_4_A=randn(1,n)*Sigma(1,5);
    r_4_A=C_Transmit_4_A+Noise_4_A;
    for j=1:size(r_4_A,2)
        if r_4_A(1,j)>0
            y_4_A(1,j)=1;
        else
            y_4_A(1,j)=0;
        end
    end
    Syndrome_4_A=y_4_A*H;
    Syndrome_4_A=mod(Syndrome_4_A,2);
    flag_4_A=0;
    for ii=1:size(Syndrome_lut_full,1)
        if isequal(Syndrome_4_A,Syndrome_lut_full(ii,:))
            C_corrected_4_A=y_4_A+e_lut_full(ii,:);
            C_corrected_4_A=mod(C_corrected_4_A,2);
            if ~isequal(C_rand_4_A,C_corrected_4_A)
                Num_errors_4_A=Num_errors_4_A+1;
            end
        end
        if ~isequal(Syndrome_4_A,Syndrome_lut_full(ii,:))
            flag_4_A=flag_4_A+1;
        end
        if flag_4_A==size(Syndrome_lut_full,1)
               Num_errors_4_A=Num_errors_4_A+1;
        end
    end
    count_error_4_A=count_error_4_A+1;
end
rate_hard_simulation_Additional(1,5)=Num_errors_4_A/count_error_4_A;
% eb/no=5 db
Num_errors_5_A=0;
count_error_5_A=0;
while(Num_errors_5_A<Max_errors)
    v_rand_5_A=(randi(2,1,k))-1;
    C_rand_5_A=v_rand_5_A*G;
    C_rand_5_A=mod(C_rand_5_A,2);
    for i=1:size(C_rand_5_A,2)
        if C_rand_5_A(1,i)==1
            C_Transmit_5_A(1,i)=1;
        else
            C_Transmit_5_A(1,i)=-1;
        end
    end
    Noise_5_A=randn(1,n)*Sigma(1,6);
    r_5_A=C_Transmit_5_A+Noise_5_A;
    for j=1:size(r_5_A,2)
        if r_5_A(1,j)>0
            y_5_A(1,j)=1;
        else
            y_5_A(1,j)=0;
        end
    end
    Syndrome_5_A=y_5_A*H;
    Syndrome_5_A=mod(Syndrome_5_A,2);
    flag_5_A=0;
    for ii=1:size(Syndrome_lut_full,1)
        if isequal(Syndrome_5_A,Syndrome_lut_full(ii,:))
            C_corrected_5_A=y_5_A+e_lut_full(ii,:);
            C_corrected_5_A=mod(C_corrected_5_A,2);
            if ~isequal(C_rand_5_A,C_corrected_5_A)
                Num_errors_5_A=Num_errors_5_A+1;
            end
        end
        if ~isequal(Syndrome_5_A,Syndrome_lut_full(ii,:))
            flag_5_A=flag_5_A+1;
        end
        if flag_5_A==size(Syndrome_lut_full,1)
               Num_errors_5_A=Num_errors_5_A+1;
        end
    end
    count_error_5_A=count_error_5_A+1;
end
rate_hard_simulation_Additional(1,6)=Num_errors_5_A/count_error_5_A;
% eb/no=6 db
Num_errors_6_A=0;
count_error_6_A=0;
while(Num_errors_6_A<Max_errors)
    v_rand_6_A=(randi(2,1,k))-1;
    C_rand_6_A=v_rand_6_A*G;
    C_rand_6_A=mod(C_rand_6_A,2);
    for i=1:size(C_rand_6_A,2)
        if C_rand_6_A(1,i)==1
            C_Transmit_6_A(1,i)=1;
        else
            C_Transmit_6_A(1,i)=-1;
        end
    end
    Noise_6_A=randn(1,n)*Sigma(1,7);
    r_6_A=C_Transmit_6_A+Noise_6_A;
    for j=1:size(r_6_A,2)
        if r_6_A(1,j)>0
            y_6_A(1,j)=1;
        else
            y_6_A(1,j)=0;
        end
    end
    Syndrome_6_A=y_6_A*H;
    Syndrome_6_A=mod(Syndrome_6_A,2);
    flag_6_A=0;
    for ii=1:size(Syndrome_lut_full,1)
        if isequal(Syndrome_6_A,Syndrome_lut_full(ii,:))
            C_corrected_6_A=y_6_A+e_lut_full(ii,:);
            C_corrected_6_A=mod(C_corrected_6_A,2);
            if ~isequal(C_rand_6_A,C_corrected_6_A)
                Num_errors_6_A=Num_errors_6_A+1;
            end
        end
        if ~isequal(Syndrome_6_A,Syndrome_lut_full(ii,:))
            flag_6_A=flag_6_A+1;
        end
        if flag_6_A==size(Syndrome_lut_full,1)
               Num_errors_6_A=Num_errors_6_A+1;
        end
    end
    count_error_6_A=count_error_6_A+1;
end
rate_hard_simulation_Additional(1,7)=Num_errors_6_A/count_error_6_A;
% eb/no=7 db
Num_errors_7_A=0;
count_error_7_A=0;
while(Num_errors_7_A<Max_errors)
    v_rand_7_A=(randi(2,1,k))-1;
    C_rand_7_A=v_rand_7_A*G;
    C_rand_7_A=mod(C_rand_7_A,2);
    for i=1:size(C_rand_7_A,2)
        if C_rand_7_A(1,i)==1
            C_Transmit_7_A(1,i)=1;
        else
            C_Transmit_7_A(1,i)=-1;
        end
    end
    Noise_7_A=randn(1,n)*Sigma(1,8);
    r_7_A=C_Transmit_7_A+Noise_7_A;
    for j=1:size(r_7_A,2)
        if r_7_A(1,j)>0
            y_7_A(1,j)=1;
        else
            y_7_A(1,j)=0;
        end
    end
    Syndrome_7_A=y_7_A*H;
    Syndrome_7_A=mod(Syndrome_7_A,2);
    flag_7_A=0;
    for ii=1:size(Syndrome_lut_full,1)
        if isequal(Syndrome_7_A,Syndrome_lut_full(ii,:))
            C_corrected_7_A=y_7_A+e_lut_full(ii,:);
            C_corrected_7_A=mod(C_corrected_7_A,2);
            if ~isequal(C_rand_7_A,C_corrected_7_A)
                Num_errors_7_A=Num_errors_7_A+1;
            end
        end
        if ~isequal(Syndrome_7_A,Syndrome_lut_full(ii,:))
            flag_7_A=flag_7_A+1;
        end
        if flag_7_A==size(Syndrome_lut_full,1)
               Num_errors_7_A=Num_errors_7_A+1;
        end
    end
    count_error_7_A=count_error_7_A+1;
end
rate_hard_simulation_Additional(1,8)=Num_errors_7_A/count_error_7_A;
semilogy(LOGS,rate_hard_simulation_Additional,'LineWidth',1,'Marker','o','Color','k')
h1=legend('CER HARD analytic','CER SOFT asympt. app.','CER SOFT union bound','CER HARD simulation','CER SOFT simulation','CER HARD full LUT','Location','best');
set(h1,'Interpreter','Latex','FontSize',9)