%% Mahyar Onsori 309823
clc;
clear all;
close all;
%% Initializing
tic;
N=16;
p=randi(2,1,N)-1;
pprime=2*p-1;
SNRdB=-10;
Ns=100000;
Ne=30;
%% Plotting T0 and T1
SNR=10^((SNRdB)/10);
Sigma=sqrt(1/SNR);
for j=1:Ns
    n=randn(1,N)*Sigma;
    r0=n;
    r1=n+pprime;
    T0=r0.*pprime;
    T0sum(j)=sum(T0);
    T1=r1.*pprime;
    T1sum(j)=sum(T1);
end
T0mean=mean(T0sum)*ones(1,Ns);
T1mean=mean(T1sum)*ones(1,Ns);
figure,plot(1:Ns,T0sum,'r');
hold on
plot(1:Ns,T0mean,'w');
title('Soft correlation: simulation outcomes. SNR=-10 dB')
xlabel('simulation number');
h1=ylabel('$$ T_0 $$ ');
set(h1,'Interpreter','Latex','FontSize',13);
axis([1 Ns -60 60]);
figure,plot(1:Ns,T1sum,1:Ns,T1mean);
title('Soft correlation: simulation outcomes. SNR=-10 dB')
xlabel('simulation number');
h2=ylabel('$$ T_1 $$ ');
set(h2,'Interpreter','Latex','FontSize',13);
axis([1 Ns -60 60]);
%% Plotting the Histograms
figure;
hist1=histogram(T0sum);
hold on
hist2=histogram(T1sum);
title('Soft correlation: PDF of T under H_0 and H_1. SNR=-10 dB')
%% Plotting Pfa and Pmd against threshold
minT0=round(min(T0sum));
maxT1=round(max(T1sum));
for thresh=minT0+abs(minT0)+1:maxT1+abs(minT0)+1
    fa(thresh)=numel(find(T0sum>thresh-abs(minT0)+1))/Ns;
    md(thresh)=numel(find(T1sum<thresh-abs(minT0)+1))/Ns;
end
d=1-md;
for i=1:length(fa)
    if log10(fa(i))<-3
        fa(i)=0;
    end
    if log10(md(i))<-3
        md(i)=0;
    end
end
theresholds=minT0:maxT1;
figure,semilogy(theresholds,fa,'ro',theresholds,md,'bo')
legend('P_{fa}','P_{md}','Location','best');
title('Soft correlation: P_{fa} and P_{md} vs. threshold. SNR=-10 dB');
h4=xlabel('threshold t');
set(h4,'Interpreter','Latex','FontSize',13);
h3=ylabel('$$ P_{fa} , P_{md} $$');
set(h3,'Interpreter','Latex','FontSize',13);
axis([minT0 maxT1 10^-8 10^0])
grid on
%% Plotting the ROC
figure,loglog(fa,d,'ko');
title('ROC curve. SNR=-10 dB , N=32');
h5=xlabel('$$P_{fa}$$');
set(h5,'Interpreter','Latex','FontSize',13);
h6=ylabel('$$ P_d $$');
set(h6,'Interpreter','Latex','FontSize',13);
axis([10^-8 10^0 10^-8 10^0])
grid on
%% Analytic PFA and PMD vs. Thresholds
PFA_Analytic=1/2*erfc(((theresholds-T0mean(1,1))/sqrt(N*Sigma^2))/sqrt(2));
PMD_Analytic=1-(1/2*erfc((((theresholds-T1mean(1,1))/sqrt(N*Sigma^2)))/sqrt(2)));
figure,semilogy(theresholds,fa,'ro',theresholds,md,'bo',theresholds,PFA_Analytic,theresholds,PMD_Analytic);
legend('P_{fa,s}','P_{md,s}','P_{fa,a}','P_{md,a}','Location','best');
title('Soft correlation: P_{fa} and P_{md} vs. threshold. SNR=-10 dB');
h7=xlabel('threshold t');
set(h7,'Interpreter','Latex','FontSize',13);
h8=ylabel('$$ P_{fa} , P_{md} $$');
set(h8,'Interpreter','Latex','FontSize',13);
axis([minT0 maxT1 10^-8 10^0])
grid on
%% Analytic ROC
PD_Analytic=1-PMD_Analytic;
figure,loglog(fa,d,'ko',PFA_Analytic,PD_Analytic,'K--');
legend('Simulation','Analytical','Location','best');
axis([10^-8 10^0 10^-8 10^0]);
grid on
title('ROC curve. SNR=-10 dB , N=32');
h9=xlabel('$$P_{fa}$$');
set(h9,'Interpreter','Latex','FontSize',13);
h10=ylabel('$$ P_d $$');
set(h10,'Interpreter','Latex','FontSize',13);
%% Pmd vs. SNR
for i=1:19
    SNRnew(i)=10^((i-11)/10);
    Sigmanew(i)=sqrt(1/SNRnew(i));
    Gamma2(i)=sqrt(Sigmanew(i)^2*N)*sqrt(2)*erfinv(1-2*0.01);
    PMD2(i)=1-(1/2*erfc(((Gamma2(i)-T1mean(1,1))/sqrt(Sigmanew(i)^2*N))/sqrt(2)));
    Gamma4(i)=sqrt(Sigmanew(i)^2*N)*sqrt(2)*erfinv(1-2*0.0001);
    PMD4(i)=1-(1/2*erfc(((Gamma4(i)-T1mean(1,1))/sqrt(Sigmanew(i)^2*N))/sqrt(2)));
end
SNRPLOT=-10:8;
figure,semilogy(SNRPLOT,PMD2,'r-o',SNRPLOT,PMD4,'b-o')
legend('P_{fa}=1e2','P_{fa}=1e4','Location','best');
grid on
title('P_{md} for P_{fa}=1e2 and P_{fa}=1e4'); 
h11=xlabel('$$ SNR[dB] $$');
set(h11,'Interpreter','Latex','FontSize',13);
h12=ylabel('$$ P_{md} $$');
set(h12,'Interpreter','Latex','FontSize',13);
toc;
