%% Mahyar Onsori 309823
clc
clear all;
close all;
%%
tic;
m=7;
L=2^m-1;
goldseq1 = comm.GoldSequence('FirstPolynomial','x^7+x^4+1',...
'SecondPolynomial','x^7+x+1',...
'FirstInitialConditions',[1 1 0 0 0 0 1],...
'SecondInitialConditions',[1 1 0 0 0 0 1],...
'Index',26,'SamplesPerFrame',L);
x1 = goldseq1()';
x1b=1-2*x1;
%%
R=ifft(fft(x1b).*conj(fft(x1b)));
R=fftshift((R));
fshift = (-63:63);
figure,plot(fshift,R)
title('Cyclic Autocorrelation');
grid on
%%
goldseq2 = comm.GoldSequence('FirstPolynomial','x^7+x^4+1',...
'SecondPolynomial','x^7+x+1',...
'FirstInitialConditions',[0 0 1 0 0 1 1],...
'SecondInitialConditions',[0 0 1 0 0 1 1],...
'Index',4,'SamplesPerFrame',L);
x2 = goldseq2()';
x2b=1-2*x2;
R2=ifft(fft(x1b).*conj(fft(x2b)));
R2=fftshift((R2));
figure,plot(fshift,R2)
title('Cross Correlation of two Gold sequences with N=127');
grid on
%%
N=127;
welch1_127=zeros(1,200);
for k=1:200
    for s=1:10
        d1=1/(k*N-1);
        d2=k*N;
        d3=nchoosek(N+s-1,s);
        d4=d1*((d2/d3)-1);
        if d4<0
            d4=0;
        end
        if ceil(N*d4^(1/(2*s)))>welch1_127(k)
            welch1_127(k)=ceil(N*d4^(1/(2*s)));
        end
    end
end
for k=1:200
    welch2_127(k)=ceil(N*(sqrt((k-1)/(k*N-1))));
end
for k=1:200
    welch3_127(k)=ceil(sqrt(N));
end
%%
sidelnikov127=zeros(1,200);
for k=1:200
    for s=0:52
        d5=(2*s+1)*(N-s);
        d6=s*(s+1)/2;
        d7=(2^s)*(N^(2*s+1));
        d8=k*factorial(2*s)*nchoosek(N,s);
        d9=d5+d6-(d7/d8);
        if d9<0
            d9=0;
        end
        if ceil(sqrt(d9))>sidelnikov127(k)
            sidelnikov127(k)=ceil(sqrt(d9));
        end
    end
end
%%
Ra=0;
Rs=0;
Rm=0;
Counter=1;
for i=1:129
    goldseq3 = comm.GoldSequence('FirstPolynomial','x^7+x^4+1',...
    'SecondPolynomial','x^7+x+1',...
    'FirstInitialConditions',[1 1 0 0 0 0 1],...
    'SecondInitialConditions',[1 1 0 0 0 0 1],...
    'Index',i-3,'SamplesPerFrame',127);
    x3 = goldseq3()';
    x3b=1-2*x3;
    Auto1=abs(ifft(fft(x3b).*conj(fft(x3b))));
    SortedAuto1=sort(Auto1);
    Maxauto1=SortedAuto1(L-1);
    if Maxauto1>Ra
        Ra=Maxauto1;
    end
    for j=i+1:129
            goldseq4 = comm.GoldSequence('FirstPolynomial','x^7+x^4+1',...
            'SecondPolynomial','x^7+x+1',...
            'FirstInitialConditions',[1 1 0 0 0 0 1],...
            'SecondInitialConditions',[1 1 0 0 0 0 1],...
            'Index',j-3,'SamplesPerFrame',127);
            x4 = goldseq4()';
            x4b=1-2*x4;
            Cross1=abs(ifft(fft(x3b).*conj(fft(x4b))));
            Maxcross1=max(max(Cross1));
            if Maxcross1>Rs
                Rs=Maxcross1;
            end
        if Ra>Rs
            Rm=Ra;
        else
            Rm=Rs;
        end
        Counter=Counter+1;
    end
end
kplot=1:200;
K=129;
figure,plot(kplot,welch1_127,'bo',kplot,welch2_127,'c',kplot,welch3_127,'k',kplot,sidelnikov127,'r',K,Rm,'kx','Markersize',12);
legend('Original Welch','First Approximation of Welch','Second Approximation of Welch','Sidelnikov','r_M computed','Location','best')
grid on
xlabel('K')
ylabel('bound')
title('N=127')
%% Truncated
Ra2=0;
Rs2=0;
Rm2=0;
Counter2=1;
for i=1:129
    goldseq5 = comm.GoldSequence('FirstPolynomial','x^7+x^4+1',...
    'SecondPolynomial','x^7+x+1',...
    'FirstInitialConditions',[1 1 0 0 0 0 1],...
    'SecondInitialConditions',[1 1 0 0 0 0 1],...
    'Index',i-3,'SamplesPerFrame',127);
    x5 = goldseq5()';
    x5b=1-2*x5;
    x6b=x5b(1:126);
    Auto2=abs(ifft(fft(x6b).*conj(fft(x6b))));
    SortedAuto2=sort(Auto2);
    Maxauto2=SortedAuto1(125);
    if Maxauto2>Ra2
        Ra2=Maxauto2;
    end
    for j=i+1:129
            goldseq7 = comm.GoldSequence('FirstPolynomial','x^7+x^4+1',...
            'SecondPolynomial','x^7+x+1',...
            'FirstInitialConditions',[1 1 0 0 0 0 1],...
            'SecondInitialConditions',[1 1 0 0 0 0 1],...
            'Index',j-3,'SamplesPerFrame',127);
            x7 = goldseq7()';
            x7b=1-2*x7;
            x8b=x7b(1:126);
            Cross2=abs(ifft(fft(x6b).*conj(fft(x8b))));
            Maxcross2=max(max(Cross2));
            if Maxcross2>Rs2
                Rs2=Maxcross2;
            end
        if Ra2>Rs2
            Rm2=Ra2;
        else
            Rm2=Rs2;
        end
        Counter2=Counter2+1;
    end
end
N2=126;
welch1_126=zeros(1,200);
for k=1:200
    for s=1:10
        d11=1/(k*N2-1);
        d22=k*N2;
        d33=nchoosek(N2+s-1,s);
        d44=d11*((d22/d33)-1);
        if d44<0
            d44=0;
        end
        if ceil(N2*d44^(1/(2*s)))>welch1_126(k)
            welch1_126(k)=ceil(N2*d44^(1/(2*s)));
        end
    end
end
for k=1:200
    welch2_126(k)=ceil(N2*(sqrt((k-1)/(k*N2-1))));
end
for k=1:200
    welch3_126(k)=ceil(sqrt(N2));
end
sidelnikov126=zeros(1,200);
for k=1:200
    for s=1:10
        d55=(2*s+1)*(N2-s);
        d66=s*(s+1)/2;
        d77=(2^s)*(N2^(2*s+1));
        d88=k*factorial(2*s)*nchoosek(N2,s);
        d99=d55+d66-(d77/d88);
        if d99<0
            d99=0;
        end
        if ceil(sqrt(d99))>sidelnikov126(k)
            sidelnikov126(k)=ceil(sqrt(d99));
        end
    end
end
figure,plot(kplot,welch1_126,'bo',kplot,welch2_126,'c',kplot,welch3_126,'k',kplot,sidelnikov126,'r',K,Rm2,'kx','Markersize',12);
legend('Original Welch','First Approximation of Welch','Second Approximation of Welch','Sidelnikov','r_M computed','Location','best')
grid on
xlabel('K')
ylabel('bound')
title('N=126')
toc;