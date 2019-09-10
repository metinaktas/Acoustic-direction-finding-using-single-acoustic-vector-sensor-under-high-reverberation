clc
clear all
close all
% Add required paths
addpath '../Simulation' '..\Utility'
%% SIMULATION PARAMETERS
acousticParameters = AdjustSimulationParameters;
multipathNumber = 1;
multipathRangesList = [0.5,1,3,5,10,20,30,40,50,100];
multipathMinAngleDifference = 30;
trialNum = 100;
AngleType = 'Free';
saveDirectory = '..\Scenarios\MultipathRangePerformance\';
saveName = sprintf('%sMultiPathNum%dAngleDif%dTrial-',AngleType,multipathNumber,multipathMinAngleDifference);
%% TRANSMIT SIGNAL
pathname = '..\Sound Samples\Gunshot Sounds\';
soundSampleName = '9mm.wav';
[signalOrj, fs] = audioread([pathname,soundSampleName],[4000,7000]);
signalOrj = signalOrj * exp(acousticParameters.Environment.Attenuation * acousticParameters.Environment.Distance);
%% RECEIVED SIGNAL PARAMETERS
deadZoneAzimuth = [0,90,180,270];
deadZoneElevation = [-80,0,80];
deadZoneRange = 10;
len = round(length(signalOrj));
%% GENERATE SCENARIOS
for trIndx = 1:trialNum
    AzimuthActualInDeadZone = mod(deadZoneAzimuth(randi(length(deadZoneAzimuth))) + (rand(1)-0.5)*2*deadZoneRange,360);
    ElevationActualInDeadZone = mod(deadZoneElevation(2) + (rand(1)-0.5)*2*deadZoneRange + 90,180)-90;
    
    indxAzStart = randi(length(deadZoneAzimuth)-1);
    startAz = deadZoneAzimuth(indxAzStart) + deadZoneRange;
    endAz = deadZoneAzimuth(indxAzStart+1) - deadZoneRange;
    if endAz < startAz
        error('endAz < startAz')
    end
    AzimuthActualOutofDeadZone = mod(rand(1) * (endAz - startAz) + startAz, 360);
    
    indxElevStart = randi(length(deadZoneElevation)-1);
    startElev = deadZoneElevation(indxElevStart) + deadZoneRange;
    endElev = deadZoneElevation(indxElevStart+1) - deadZoneRange;
    if endElev < startElev
        error('endElev < startElev')
    end
    ElevationActualOutofDeadZone = mod(rand(1) * (endElev - startElev) + startElev + 90, 180) - 90;
    switch AngleType
        case 'SafeAzSafeEl'
            AzimuthActual = mod(AzimuthActualOutofDeadZone + [0,(rand(1,multipathNumber)-0.5)*180],360);
            ElevationActual = mod(ElevationActualOutofDeadZone + [0,(rand(1,multipathNumber)-0.5)*180]+90,180)-90;
        case 'SafeAzDeadEl'
            AzimuthActual = mod(AzimuthActualOutofDeadZone + [0,(rand(1,multipathNumber)-0.5)*180],360);
            ElevationActual = mod(ElevationActualInDeadZone + [0,(rand(1,multipathNumber)-0.5)*180]+90,180)-90;        
        case 'DeadAzSafeEl'
            AzimuthActual = mod(AzimuthActualInDeadZone + [0,(rand(1,multipathNumber)-0.5)*180],360);
            ElevationActual = mod(ElevationActualOutofDeadZone + [0,(rand(1,multipathNumber)-0.5)*180]+90,180)-90;
        case 'DeadAzDeadEl'
            AzimuthActual = mod(AzimuthActualInDeadZone + [0,(rand(1,multipathNumber)-0.5)*180],360);
            ElevationActual = mod(ElevationActualInDeadZone + [0,(rand(1,multipathNumber)-0.5)*180]+90,180)-90;
        case 'Free'
            AzimuthActual = rand(1,multipathNumber+1)*360;
            ElevationActual = [(rand(1)-0.5)*2*75,(rand(1,multipathNumber)-0.5)*2*90];
        otherwise
            error('AngleType is not defined properly')
    end
    indx = find(abs(AzimuthActual - AzimuthActual(1)) <= multipathMinAngleDifference);
    if length(indx) > 1
        indx = indx(2:end);
        AzimuthActual(indx) = mod(AzimuthActual(1) + multipathMinAngleDifference * sign(round(AzimuthActual(indx) - AzimuthActual(1)) + 0.5),360);
    end
    indx = find(abs(ElevationActual - ElevationActual(1)) <= multipathMinAngleDifference);
    if length(indx) > 1
        indx = indx(2:end);
        ElevationActual(indx) = mod(ElevationActual(1) + multipathMinAngleDifference * sign(round(ElevationActual(indx) - ElevationActual(1)) + 0.5) + 90,180)-90;
    end
    for rIndx = 1:length(multipathRangesList)
        if multipathRangesList(rIndx) > 50
            multipathRanges = 0;            
        else
            multipathRanges = [0,multipathRangesList(rIndx)];
        end
        
        [~, avsDataNoiseless{rIndx}, noiseOrg{rIndx}] = SimulateAvsDataOutput(acousticParameters, signalOrj, len, AzimuthActual(1:length(multipathRanges)), ElevationActual(1:length(multipathRanges)), fs, multipathRanges);
    end
    save([saveDirectory,saveName,sprintf('%d',trIndx)],'acousticParameters','noiseOrg','avsDataNoiseless','fs','AzimuthActual','ElevationActual','multipathRangesList');
    disp(sprintf('Trial %d / %d was saved successfully\n',trIndx,trialNum));
end
