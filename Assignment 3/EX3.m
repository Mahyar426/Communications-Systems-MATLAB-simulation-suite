%% Mahyar Onsori 309823
clc;
clear all;
close all;
%% Parameters
Delta=1;
N=128;
a1=1;a2=0.5;a3=0.9;a4=0.3;
A1=rand*a1*(2*(randi(2,1,1)-1)-1);A2=rand*a2*(2*(randi(2,1,1)-1)-1);A3=rand*a3*(2*(randi(2,1,1)-1)-1);A4=rand*a4*(2*(randi(2,1,1)-1)-1);
B1=sqrt(a1^2-A1^2);B2=sqrt(a2^2-A2^2);B3=sqrt(a3^2-A3^2);B4=sqrt(a4^2-A4^2);
C1=complex(A1,B1);C2=complex(A2,B2);C3=complex(A3,B3);C4=complex(A4,B4);
D1=1;D2=1.01;D3=1.015;D4=1.02;
%% Attenuation Profile
for f=1:128
    H(f)=C1*exp(-i*2*pi*f*D1)+C2*exp(-i*2*pi*f*D2)+C3*exp(-i*2*pi*f*D3)+C4*exp(-i*2*pi*f*D4);
    Att(f)=1/(abs(H(f))^2);
end
fplot=1:128;
Attdb=10*log10(Att);
figure,subplot(4,1,1),plot(fplot,Attdb,'b');
title('Attenuation')
xlabel('f');
ylabel('Att[dB]');
Maxatt=max(Attdb);
Minatt=min(Attdb);
Maxattplot=max(Maxatt,abs(Minatt));
maxplotdb=Maxattplot-mod(Maxattplot,5)+5;
minplotdb=-maxplotdb;
axis([0 128 minplotdb maxplotdb])
grid on
%% Hartogs Hughes Algorithm
PTx=median(Att)*N;
PTx2=PTx;
[Sorted,Indices]=sort(Att);
Combined=[];
K=20;
for i=1:K
    Combined=[Combined,i*Sorted];
end
[SortedCombined,IndicesCombined]=sort(Combined);
Tones=zeros(1,N);
for i=1:N*K
    if PTx-sum(SortedCombined(1:i))<0
        Stop=i;
        break;
    end
end
for i=1:Stop-1
        PTx=PTx-SortedCombined(i);
        if(mod(IndicesCombined(i),128)~=0)
            Tones(Indices(mod(IndicesCombined(i),128)))=Tones(Indices(mod(IndicesCombined(i),128)))+1;
        else
            Tones(Indices(128))=Tones(Indices(128))+1;
        end
        if (PTx<0)
            break;
        end
end
NumActiveTones=0;
ActiveTones=zeros(1,128);
TotalBits=0;
for i=1:N
    if Tones(i)>0
        NumActiveTones=NumActiveTones+1;
        TotalBits=Tones(i)+TotalBits;
        ActiveTones(i)=1;
    end
end
subplot(4,1,2),plot(fplot,Tones,'ro');
title(['HH, Active tones= ',num2str(NumActiveTones),', Total Bits= ',num2str(TotalBits)]);
xlabel('f');
ylabel('bit');
grid on
Max=max(Tones);
axis([0 128 0 Max+1])
%% Uniform on Active Tones
PTxmean=(PTx2/NumActiveTones)*ActiveTones;
Tones2=zeros(1,N);
for i=1:N
    j=1;
    while(PTxmean(i)>0)
        PTxmean(i)=PTxmean(i)-j*Att(i);
        Tones2(i)=Tones2(i)+1;
        j=j+1;
    end
    if PTxmean(i)<0
        Tones2(i)=Tones2(i)-1;
    end
end
NumActiveTones2=0;
ActiveTones2=zeros(1,128);
TotalBits2=0;
for i=1:N
    if Tones2(i)>0
        NumActiveTones2=NumActiveTones2+1;
        TotalBits2=Tones2(i)+TotalBits2;
        ActiveTones2(i)=1;
    end
end
subplot(4,1,3),plot(fplot,Tones2,'go');
title(['Uniform on HH, Active tones= ',num2str(NumActiveTones2),', Total Bits= ',num2str(TotalBits2)]);
xlabel('f');
ylabel('bit');
axis([0 128 0 Max+1])
grid on
%% Uniform on All Tones
PTxmean2=(PTx2/N)*ones(1,N);
Tones3=zeros(1,N);
for i=1:N
    j=1;
    while(PTxmean2(i)>0)
        PTxmean2(i)=PTxmean2(i)-j*Att(i);
        Tones3(i)=Tones3(i)+1;
        j=j+1;
    end
    if PTxmean2(i)<0
        Tones3(i)=Tones3(i)-1;
    end
end
NumActiveTones3=0;
ActiveTones3=zeros(1,128);
TotalBits3=0;
for i=1:N
    if Tones3(i)>0
        NumActiveTones3=NumActiveTones3+1;
        TotalBits3=Tones3(i)+TotalBits3;
        ActiveTones3(i)=1;
    end
end
subplot(4,1,4),plot(fplot,Tones3,'co');
title(['Uniform on HH, Active tones= ',num2str(NumActiveTones3),', Total Bits= ',num2str(TotalBits3)]);
xlabel('f');
ylabel('bit');
grid on
axis([0 128 0 Max+1])