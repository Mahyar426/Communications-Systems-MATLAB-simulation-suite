%% Mahyar Onsori 309823
clc;
clear all;
close all;
%% Creating a Scenario
starttime=datetime(2023,01,01,00,00,0);
stoptime=starttime+days(7);
sampletime=60;
sc=satelliteScenario(starttime,stoptime,sampletime);
%% Creating a Satelite
semiMajorAxis=earthRadius+500e3;
inclination=97.4;
eccentrcity=0;
rightAscensionOfAscendingNode=0;
argumentOfPeriapsis=0;
trueAnomaly=0;
name='EO Sat';
Satellite=satellite(sc,semiMajorAxis,eccentrcity,inclination,rightAscensionOfAscendingNode,argumentOfPeriapsis,trueAnomaly,'Name',name);
%% Adding Gimbal To the Satellite
gimbaltxSat = gimbal(Satellite);%Correct
%% Adding a Transmitter To the Gimbal of the Satellite
frequency_txsat= 8.212e9;
power_txsat=3;
bitRate_txsat= 120; 
systemLoss_txsat= 1;
txSat = transmitter(gimbaltxSat,Name="Satellite Transmitter",Frequency=frequency_txsat,power=power_txsat,BitRate=bitRate_txsat,SystemLoss=systemLoss_txsat);
%% Adding Gaussian Antenna To the Satellite
dishDiameter_sat = 0.1818;
apertureEfficiency_sat = 0.65;
gaussianAntenna(txSat,DishDiameter=dishDiameter_sat,ApertureEfficiency=apertureEfficiency_sat);
%% Creating Ground Station #1
lat1=68.35;
lon1=133.72;
alt1=15;
elevationdegree1=10;
gs1 = groundStation(sc,'Name','Inuvik','Latitude',lat1,'Longitude',lon1,'Altitude',alt1,'MinElevationAngle',elevationdegree1);
%% Creating Ground Station #2
lat2=78.22;
lon2=15.38;
alt2=440;
elevationdegree2=10;
gs2 = groundStation(sc,'Name','Svalbard','Latitude',lat2,'Longitude',lon2,'Altitude',alt2,'MinElevationAngle',elevationdegree2);
%% Creating Ground Station #3
lat3=-46.52;
lon3=168.48;
alt3=0;
elevationdegree3=10;
gs3 = groundStation(sc,'Name','Awarua','Latitude',lat3,'Longitude',lon3,'Altitude',alt3,'MinElevationAngle',elevationdegree3);
%% Creating Ground Station #4
lat4=-72.01;
lon4=2.53;
alt4=1200;
elevationdegree4=10;
gs4 = groundStation(sc,'Name','Troll','Latitude',lat4,'Longitude',lon4,'Altitude',alt4,'MinElevationAngle',elevationdegree4);
%% Adding Gimbal, Receiver, and Antenna to Ground Station #1
gimbalgs1= gimbal(gs1,'MountingAngles', [0; 180;0]);
gainToNoiseTemperatureRatio_rxgs1=25;                                       
systemLoss_rxgs1=1;
Requiredebno_rxgs1=18;
Prereceiverloss_rxgs1=0;
rxGs1=receiver(gimbalgs1,Name="Ground Station 1 Receiver",GainToNoiseTemperatureRatio=gainToNoiseTemperatureRatio_rxgs1,SystemLoss=systemLoss_rxgs1,PreReceiverLoss=Prereceiverloss_rxgs1,RequiredEbNo=Requiredebno_rxgs1);
dishDiameter_gs = 3.7;
apertureEfficiency_gs=0.65;
gaussianAntenna(rxGs1,DishDiameter=dishDiameter_gs,ApertureEfficiency=apertureEfficiency_gs);
%% Adding Gimbal, Receiver, and Antenna to Ground Station #2
gimbalgs2= gimbal(gs2,'MountingAngles', [0; 180;0]);
gainToNoiseTemperatureRatio_rxgs2=25;                                      
systemLoss_rxgs2=1;
Requiredebno_rxgs2=18;
Prereceiverloss_rxgs2=0;
rxGs2=receiver(gimbalgs2,Name="Ground Station 2 Receiver",GainToNoiseTemperatureRatio=gainToNoiseTemperatureRatio_rxgs2,SystemLoss=systemLoss_rxgs2,PreReceiverLoss=Prereceiverloss_rxgs2,RequiredEbNo=Requiredebno_rxgs2);
gaussianAntenna(rxGs2,DishDiameter=dishDiameter_gs,ApertureEfficiency=apertureEfficiency_gs);
%% Adding Gimbal, Receiver, and Antenna to Ground Station #3
gimbalgs3= gimbal(gs3,'MountingAngles', [0; 180;0]);
gainToNoiseTemperatureRatio_rxgs3=25;                                      
systemLoss_rxgs3=1;
Requiredebno_rxgs3=18;
Prereceiverloss_rxgs3=0;
rxGs3=receiver(gimbalgs3,Name="Ground Station 3 Receiver",GainToNoiseTemperatureRatio=gainToNoiseTemperatureRatio_rxgs3,SystemLoss=systemLoss_rxgs3,PreReceiverLoss=Prereceiverloss_rxgs3,RequiredEbNo=Requiredebno_rxgs3);
gaussianAntenna(rxGs3,DishDiameter=dishDiameter_gs,ApertureEfficiency=apertureEfficiency_gs);
%% Adding Gimbal, Receiver, and Antenna to Ground Station #4
gimbalgs4= gimbal(gs4,'MountingAngles', [0; 180;0]);
gainToNoiseTemperatureRatio_rxgs4=25;                                      
systemLoss_rxgs4=1;
Requiredebno_rxgs4=18;
Prereceiverloss_rxgs4=0;
rxGs4=receiver(gimbalgs4,Name="Ground Station 4 Receiver",GainToNoiseTemperatureRatio=gainToNoiseTemperatureRatio_rxgs4,SystemLoss=systemLoss_rxgs4,PreReceiverLoss=Prereceiverloss_rxgs4,RequiredEbNo=Requiredebno_rxgs4);
gaussianAntenna(rxGs4,DishDiameter=dishDiameter_gs,ApertureEfficiency=apertureEfficiency_gs);
%% Pointing the Gimbals at the Satellite
pointAt(gimbalgs1,Satellite);
pointAt(gimbalgs2,Satellite);
pointAt(gimbalgs3,Satellite);
pointAt(gimbalgs4,Satellite);
%pointAt(gimbaltxSat,gs1);
%pointAt(gimbaltxSat,gs2);
pointAt(gimbaltxSat,gs3);
%pointAt(gimbaltxSat,gs4);
%% Links
lnk1 = link(txSat,rxGs1);
lnk2 = link(txSat,rxGs2);
lnk3 = link(txSat,rxGs3);
lnk4 = link(txSat,rxGs4);
intvls_lnk1=linkIntervals(lnk1);
intvls_lnk2=linkIntervals(lnk2);
intvls_lnk3=linkIntervals(lnk3);
intvls_lnk4=linkIntervals(lnk4);
%% Accesses
ac1 = access(Satellite,gs1);
ac2 = access(Satellite,gs2);
ac3 = access(Satellite,gs3);
ac4 = access(Satellite,gs4);
intvls_ac1 = accessIntervals(ac1);
intvls_ac2 = accessIntervals(ac2);
intvls_ac3 = accessIntervals(ac3);
intvls_ac4 = accessIntervals(ac4);
Total_Connectivity1=sum(intvls_ac1.Duration);
Total_Connectivity2=sum(intvls_ac2.Duration);
Total_Connectivity3=sum(intvls_ac3.Duration);
Total_Connectivity4=sum(intvls_ac4.Duration);
Total_Connectivity=Total_Connectivity1+Total_Connectivity2+Total_Connectivity3+Total_Connectivity4;
Mean_Daily_Connectivity=Total_Connectivity/7;
%% EbNo
[e1,t1]=ebno(lnk1);
e1=isfinite(e1).*e1;
[e2,t2]=ebno(lnk2);
e2=isfinite(e2).*e2;
[e3,t3]=ebno(lnk3);
e3=isfinite(e3).*e3;
[e4,t4]=ebno(lnk4);
e4=isfinite(e4).*e4;
Max1=max(e1);
Min1=min(e1);
Max2=max(e2);
Min2=min(e2);
Max3=max(e3);
Min3=min(e3);
Max4=max(e4);
Min4=min(e4);
%%
PLcfgP618_1 = p618Config('Frequency', 8.212e9 , 'ElevationAngle', 10, 'Latitude', lat1 ,'Longitude', lon1 ,'TotalAnnualExceedance', 0.1 ,'AntennaDiameter', dishDiameter_gs,'AntennaEfficiency', apertureEfficiency_gs);
[PL_1, ~, TSKY_1] = p618PropagationLosses(PLcfgP618_1);
PLcfgP618_2 = p618Config('Frequency', 8.212e9 , 'ElevationAngle', 10, 'Latitude', lat2 ,'Longitude', lon2 ,'TotalAnnualExceedance', 0.1 ,'AntennaDiameter', dishDiameter_gs,'AntennaEfficiency', apertureEfficiency_gs);
[PL_2, ~, TSKY_2] = p618PropagationLosses(PLcfgP618_2);
PLcfgP618_3 = p618Config('Frequency', 8.212e9 , 'ElevationAngle', 10, 'Latitude', lat3 ,'Longitude', lon3 ,'TotalAnnualExceedance', 0.1 ,'AntennaDiameter', dishDiameter_gs,'AntennaEfficiency', apertureEfficiency_gs);
[PL_3, ~, TSKY_3] = p618PropagationLosses(PLcfgP618_3);
PLcfgP618_4 = p618Config('Frequency', 8.212e9 , 'ElevationAngle', 10, 'Latitude', lat4 ,'Longitude', lon4 ,'TotalAnnualExceedance', 0.1 ,'AntennaDiameter', dishDiameter_gs,'AntennaEfficiency', apertureEfficiency_gs);
[PL_4, ~, TSKY_4] = p618PropagationLosses(PLcfgP618_4);
%% Creating Rainy Receivers
gainToNoiseTemperatureRatio_rxgs1_rainy=gainToNoiseTemperatureRatio_rxgs1-10*log10(TSKY_1);   %correct                                       
Requiredebno_rxgs1_rainy=3;
rxGs1_rainy=receiver(gimbalgs1,Name="Ground Station 1 Rainy Receiver",GainToNoiseTemperatureRatio=gainToNoiseTemperatureRatio_rxgs1_rainy,SystemLoss=systemLoss_rxgs3,PreReceiverLoss=Prereceiverloss_rxgs3,RequiredEbNo=Requiredebno_rxgs1_rainy);
gaussianAntenna(rxGs1,DishDiameter=dishDiameter_gs,ApertureEfficiency=apertureEfficiency_gs);%Correct

gainToNoiseTemperatureRatio_rxgs2_rainy=gainToNoiseTemperatureRatio_rxgs2-10*log10(TSKY_2);   %correct                                       
Requiredebno_rxgs2_rainy=3;
rxGs2_rainy=receiver(gimbalgs2,Name="Ground Station 2 Rainy Receiver",GainToNoiseTemperatureRatio=gainToNoiseTemperatureRatio_rxgs2_rainy,SystemLoss=systemLoss_rxgs3,PreReceiverLoss=Prereceiverloss_rxgs3,RequiredEbNo=Requiredebno_rxgs2_rainy);
gaussianAntenna(rxGs2,DishDiameter=dishDiameter_gs,ApertureEfficiency=apertureEfficiency_gs);%Correct

gainToNoiseTemperatureRatio_rxgs3_rainy=gainToNoiseTemperatureRatio_rxgs3-10*log10(TSKY_3);   %correct                                       
Requiredebno_rxgs3_rainy=3;
rxGs3_rainy=receiver(gimbalgs3,Name="Ground Station 3 Rainy Receiver",GainToNoiseTemperatureRatio=gainToNoiseTemperatureRatio_rxgs3_rainy,SystemLoss=systemLoss_rxgs3,PreReceiverLoss=Prereceiverloss_rxgs3,RequiredEbNo=Requiredebno_rxgs3_rainy);
gaussianAntenna(rxGs3,DishDiameter=dishDiameter_gs,ApertureEfficiency=apertureEfficiency_gs);%Correct

gainToNoiseTemperatureRatio_rxgs4_rainy=gainToNoiseTemperatureRatio_rxgs4-10*log10(TSKY_4);   %correct                                       
Requiredebno_rxgs4_rainy=3;
rxGs4_rainy=receiver(gimbalgs4,Name="Ground Station 4 Rainy Receiver",GainToNoiseTemperatureRatio=gainToNoiseTemperatureRatio_rxgs4_rainy,SystemLoss=systemLoss_rxgs3,PreReceiverLoss=Prereceiverloss_rxgs3,RequiredEbNo=Requiredebno_rxgs4_rainy);
gaussianAntenna(rxGs4,DishDiameter=dishDiameter_gs,ApertureEfficiency=apertureEfficiency_gs);%Correct
%% Rainy Links
lnk1_rainy = link(txSat,rxGs1_rainy);
lnk2_rainy = link(txSat,rxGs2_rainy);
lnk3_rainy = link(txSat,rxGs3_rainy);
lnk4_rainy = link(txSat,rxGs4_rainy);
intvls_lnk1_rainy=linkIntervals(lnk1_rainy);
intvls_lnk2_rainy=linkIntervals(lnk2_rainy);
intvls_lnk3_rainy=linkIntervals(lnk3_rainy);
intvls_lnk4_rainy=linkIntervals(lnk4_rainy);
[e1_rainy,~]=ebno(lnk1_rainy);
e1_rainy=isfinite(e1_rainy).*e1_rainy;
e1_rainy=e1_rainy-PL_1.At;
[e2_rainy,~]=ebno(lnk2_rainy);
e2_rainy=isfinite(e2_rainy).*e2_rainy;
e2_rainy=e2_rainy-PL_2.At;
[e3_rainy,~]=ebno(lnk3_rainy);
e3_rainy=isfinite(e3_rainy).*e3_rainy;
e3_rainy=e3_rainy-PL_3.At;
[e4_rainy,~]=ebno(lnk4_rainy);
e4_rainy=isfinite(e4_rainy).*e4_rainy;
e4_rainy=e4_rainy-PL_4.At;
Max1_rainy=max(e1_rainy);
Min1_rainy=min(e1_rainy);
Max2_rainy=max(e2_rainy);
Min2_rainy=min(e2_rainy);
Max3_rainy=max(e3_rainy);
Min3_rainy=min(e3_rainy);
Max4_rainy=max(e4_rainy);
Min4_rainy=min(e4_rainy);
Total_Connectivity1_rainy=sum(intvls_lnk1_rainy.Duration);
Total_Connectivity2_rainy=sum(intvls_lnk2_rainy.Duration);
Total_Connectivity3_rainy=sum(intvls_lnk3_rainy.Duration);
Total_Connectivity4_rainy=sum(intvls_lnk4_rainy.Duration);
Total_Connectivity_rainy=Total_Connectivity1_rainy+Total_Connectivity2_rainy+Total_Connectivity3_rainy+Total_Connectivity4_rainy;
Mean_Daily_Connectivity_rainy=Total_Connectivity_rainy/7;
%% Initialize
LDRS.planes=3;
LDRS.satsPerPlane=6;
LDRS.numSats= LDRS.planes*LDRS.satsPerPlane;
LDRS.eccentricity=0;
LDRS.inclination=97.4;
LDRS.altitude=1000e3;
LDRS.semiMajorAxis= earthRadius ('meters') + LDRS.altitude;
LDRS.argOfPeriapsis=0;
%% Walker
LDRS.planePhase= 360/LDRS.planes;
LDRS.satsPhase=360/ LDRS.satsPerPlane;
LDRS.RAAN = zeros(1,LDRS.numSats);
LDRS.trueanomaly= zeros(1, LDRS.numSats);
idx=0;
for plane = 1:LDRS.planes
    for sat=1:LDRS.satsPerPlane
        idx= idx + 1;
        LDRS.RAAN(idx)= plane * LDRS.planePhase;
        LDRS.trueanomaly(idx) = sat * LDRS.satsPhase;
    end
end
LDRS.Satellites=satellite(sc,LDRS.semiMajorAxis * ones(1, LDRS.numSats ),LDRS.eccentricity* ones(1, LDRS.numSats ),LDRS.inclination* ones(1, LDRS.numSats ),LDRS.RAAN,LDRS.argOfPeriapsis* ones(1, LDRS.numSats ),LDRS.trueanomaly);
% Add a conical sensor to gimbal
LCTGimbal=gimbal(Satellite);
LCT(1)=conicalSensor (LCTGimbal, 'MaxViewAngle', 90);
% To compute when the LDRS is in visibility of the satellite
ac_sat_sats=access(Satellite,LDRS.Satellites(1));
ac_sat_gs1=access(gs1,LDRS.Satellites(1));
ac_sat_gs2=access(gs2,LDRS.Satellites(1));
ac_sat_gs3=access(gs3,LDRS.Satellites(1));
ac_sat_gs4=access(gs4,LDRS.Satellites(1));
intvls_ac_sat_sats=accessIntervals(ac_sat_sats);
intvls_ac_sat_gs1=accessIntervals(ac_sat_gs1);
intvls_ac_sat_gs2=accessIntervals(ac_sat_gs2);
intvls_ac_sat_gs3=accessIntervals(ac_sat_gs3);
intvls_ac_sat_gs4=accessIntervals(ac_sat_gs4);
Total_Connectivity1_LDRS=sum(intvls_ac_sat_gs1.Duration);
Total_Connectivity2_LDRS=sum(intvls_ac_sat_gs2.Duration);
Total_Connectivity3_LDRS=sum(intvls_ac_sat_gs3.Duration);
Total_Connectivity4_LDRS=sum(intvls_ac_sat_gs4.Duration);
Total_Connectivitysat_LDRS=sum(intvls_ac_sat_sats.Duration);
Total_Connectivity_LDRS=Total_Connectivity1_LDRS+Total_Connectivity2_LDRS+Total_Connectivity3_LDRS+Total_Connectivity4_LDRS;
Mean_Daily_Connectivity_LDRS=Total_Connectivity_LDRS/7;
%% Playing Scenario
%play(sc);