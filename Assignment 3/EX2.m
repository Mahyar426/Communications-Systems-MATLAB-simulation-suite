%% Mahyar Onsori 309823
clc;
clear all;
close all;
%% One User Scenario
N=8;
% Analytical
for ebnodb=1:21
        ebno(ebnodb)=10^((ebnodb-1)/10);
        P(ebnodb)=1/2*erfc(sqrt(ebno(ebnodb)));
end
% Simulation
c=randi(2,1,N)-1;
cprime=2*c-1;
for ebnodb=1:8
    ebno(ebnodb)=10^((ebnodb-1)/10);
    Sigma2(ebnodb)= N/(2*ebno(ebnodb));
    Num_errors=0;
    count_error=0;
    while(Num_errors<500)
        v=randi(2,1,1)-1;
        vprime=2*v-1;
        Noise=randn(1,N)*sqrt(Sigma2(ebnodb));
        r=cprime*vprime+Noise;
        S=sum(r.*cprime);
        if S>0
            vg=1;
        else
            vg=-1;
        end
            if ~isequal(vprime,vg)
                Num_errors=Num_errors+1;
            end
        count_error=count_error+1; 
        rate(ebnodb)=Num_errors/count_error;
    end
end

%% Two User Scenario
c1=randi(2,1,N)-1;
c2=randi(2,1,N)-1;
cprime1=2*c1-1;
cprime2=2*c2-1;
Projection=sum(cprime1.*cprime2);
Positive=((N+Projection)/N)^2;
Negative=((N-Projection)/N)^2;
% Analytical
for ebnodb=1:21
        Pb1e=1/2*erfc(sqrt(Positive*ebno(ebnodb)));
        Pb2e=1/2*erfc(sqrt(Negative*ebno(ebnodb)));
        Pe(ebnodb)=(Pb1e+Pb2e)/2;
end
% Simulation
for ebnodb=1:8
    Num_errors2=0;
    count_error2=0;
    while(Num_errors2<500)
        v1=randi(2,1,1)-1;
        v2=randi(2,1,1)-1;
        vprime1=2*v1-1;
        vprime2=2*v2-1;
        Noise=randn(1,N)*sqrt(Sigma2(ebnodb));
        rr=cprime1*vprime1+cprime2*vprime2+Noise;
        S1=sum(rr.*cprime1);
        if S1>0
            vg1=1;
        else
            vg1=-1;
        end
        if ~isequal(vprime1,vg1)
            Num_errors2=Num_errors2+1;
        end
        count_error2=count_error2+1; 
    end
    rate2(ebnodb)=Num_errors2/count_error2;
end
%% Two user scenario with cyclic shift
% Simulation
for ebnodb=1:8
    for Shift=1:8
        cprime3=circshift(cprime2,Shift);
        Projection2(Shift)=sum(cprime1.*cprime3);
        Num_errors3=0;
        count_error3=0;
        while(Num_errors3<500)
            v3=randi(2,1,1)-1;
            v4=randi(2,1,1)-1;
            vprime3=2*v3-1;
            vprime4=2*v4-1;
            Noise=randn(1,N)*sqrt(Sigma2(ebnodb));
            rrr=cprime1*vprime3+cprime3*vprime4+Noise;
            S2=sum(rrr.*cprime1);
            if S2>0
                vg2=1;
            else
                vg2=-1;
            end
            if ~isequal(vprime3,vg2)
                Num_errors3=Num_errors3+1;
            end
            count_error3=count_error3+1; 
        end
        rate3(Shift)=Num_errors3/count_error3;
    end
    AverageProjection=sum(abs(Projection2))/8;
    rate4(ebnodb)=sum(rate3)/8;
end
%% Plotting
Ebnoplot=0:20;
Ebnoplot2=0:7;
figure,semilogy(Ebnoplot,P,'g',Ebnoplot,Pe,'r',Ebnoplot2,rate,'bx',Ebnoplot2,rate2,'ko');
legend('Analytical 2-PAM','Analytical 2 Users','Simulation 2-PAM','Simulation 2 Users')
title(['CDMA with 2 users and p=' , num2str(Projection),' and N=',num2str(N)]);
xlabel('E_{b}/N_{0}')
ylabel('P_{b}(e)')
ylim([10^-10,1])
grid on
figure,semilogy(Ebnoplot,P,'g',Ebnoplot,Pe,'r',Ebnoplot2,rate,'bx',Ebnoplot2,rate2,'ko',Ebnoplot2,rate4,'m*');
legend('Analytical 2-PAM','Analytical 2 Users','Simulation 2-PAM','Simulation 2 Users','Simulation 2 Users with shift')
xlabel('E_{b}/N_{0}')
ylabel('P_{b}(e)')
ylim([10^-10,1])
grid on