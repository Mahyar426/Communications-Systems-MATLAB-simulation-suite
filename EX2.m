%% Mahyar Onsori 309823
clc;
clear all;
close all;
% g(D)=p(D).(D+1)= D^4+D^3+D^2+1
tic;
crc_encoder = comm.CRCGenerator('Polynomial','z^4 + z^3 + z^2 + 1');
k2=2;
k3=3;
k4=4;
k5=5;
k6=6;
k7=7;
k8=8;
[dist_profile2,dist2,Ai2]=DistanceProfile(k2,crc_encoder);
[dist_profile3,dist3,Ai3]=DistanceProfile(k3,crc_encoder);
[dist_profile4,dist4,Ai4]=DistanceProfile(k4,crc_encoder);
[dist_profile5,dist5,Ai5]=DistanceProfile(k5,crc_encoder);
[dist_profile6,dist6,Ai6]=DistanceProfile(k6,crc_encoder);
[dist_profile7,dist7,Ai7]=DistanceProfile(k7,crc_encoder);
[dist_profile8,dist8,Ai8]=DistanceProfile(k8,crc_encoder);
toc;