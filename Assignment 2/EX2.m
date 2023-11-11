%% Mahyar Onsori 309823
clc;
clear all;
close all;
%% Initializing
tic;
m=7;
Nb=2^m-1;
pnSequence = comm.PNSequence('Polynomial',[7 4 0],'SamplesPerFrame',Nb,'InitialConditions',[1 1 1 0 1 1 0]);
vn = pnSequence()';
L=length(vn);
%% Generate and plot the 127-symbol bipolar sequence
vnprime=2*vn-1;
xplot=1:127;
plot(xplot,vnprime);
axis([0 128 -2 2])
h1=title('$$ 5G PSS m-sequence, 127 symbols, N_1=64, N_0=63, N_T=64, N_{NT}=63$$');
set(h1,'Interpreter','Latex','FontSize',10);
%% Creating a table with N0 and N1 and NT and NNT
N0=0;N1=0;NT=0;NNT=0;
for i=1:L
    if vn(i)==1
        N1=N1+1;
    else
        N0=N0+1;
    end
    if i>1
        if vn(i)~=vn(i-1)
            NT=NT+1;
        else
            NNT=NNT+1;
        end
    end
    if i==1
        if vn(i)~=vn(L)
            NT=NT+1;
        else
            NNT=NNT+1;
        end
    end
end
N1N0=N1/N0;
NTNNT=NT/NNT;
%% Runs
cnt1=1;
cnt2=1;
NR0=ones(1,L);
NR1=ones(1,L);
while(cnt1<L+1)
    while(vn(cnt2)==vn(cnt2+1))
        if vn(cnt2)~=vn(cnt2-1)
            Flag=cnt2;
        end
        if vn(cnt2)==0
            NR0(Flag)=NR0(Flag)+1;
        end
        if vn(cnt2)==1
            NR1(Flag)=NR1(Flag)+1;
        end
        cnt2=cnt2+1;
    end
    cnt1=cnt1+1;
    cnt2=cnt1;
    Flag=cnt1;
    if cnt1==L
        if vn(L)~=vn(L-1)
            if vn(length(vn))==0
                NR0(L)=1;
                break;
            else
                NR1(L)=1;
                break;
            end
        end
    end
end
for i=1:L-1
    if NR1(i)>1
        J=NR1(i);
        NR1(i+1:i+J)=0;
        NR0(i:i+J-1)=0;
    end
    if NR0(i)>1
        K=NR0(i);
        NR0(i+1:i+K)=0;
        NR1(i:i+K-1)=0;
    end
    if NR1(i)<NR1(i+1)
        NR1(i)=0;
    end
    if NR0(i)<NR0(i+1)
        NR0(i)=0;
    end
    if NR0(i)==1
        cnt3=i;
        while(NR0(cnt3)==1)
            if(NR0((cnt3)+1)==1)
                NR0(cnt3+1)=0;
                cnt3=cnt3+1;
            end
        end
    end
    if NR1(i)==1
        cnt4=i;
        while(NR1(cnt4)==1)
            if(NR1((cnt4)+1)==1)
                NR1(cnt4+1)=0;
                cnt4=cnt4+1;
            end
        end
    end

end
NR0RUN=zeros(1,round(log2(L)));
NR1RUN=zeros(1,round(log2(L)));
NR=0;
for i=1:L
    if NR0(i)>0
        NR=NR+1;
    end
    if(NR1(i)>0)
        NR=NR+1;
    end
end

for i=1:L
    if(NR0(i)>0)
        NR0RUN(NR0(i))=NR0RUN(NR0(i))+1;
    end
end
for i=1:L
    if(NR1(i)>0)
        NR1RUN(NR1(i))=NR1RUN(NR1(i))+1;
    end
end
NR0RUNsum=sum(NR0RUN);
NR1RUNsum=sum(NR1RUN);
figure,bar(xplot,NR0);
hold on
grid on
bar(xplot,NR1)
h2=xlabel('Starting position in the sequence');
h3=ylabel('Run length');
h4=title('5G PSS m-sequence, 127 symbols');
%% Autocorrelation
R=ifft(fft(vnprime).*conj(fft(vnprime)));
R=fftshift((R));
fshift = (-63:63);
figure,plot(fshift,R)
grid on
%% Truncated
RTruncated=ifft(fft(vnprime(1:117)).*conj(fft(vnprime(1:117))));
MPSL=max(RTruncated(1,2:end));
pnSequence2 = comm.PNSequence('Polynomial',[7 4 0],'SamplesPerFrame',Nb,'InitialConditions',[1 0 1 0 0 0 1]);
vn2 = pnSequence2()';
vnprime2=2*vn2-1;
RTruncated2=ifft(fft(vnprime2(1:117)).*conj(fft(vnprime2(1:117))));
MPSL2=max(RTruncated2(1,2:end));
toc;