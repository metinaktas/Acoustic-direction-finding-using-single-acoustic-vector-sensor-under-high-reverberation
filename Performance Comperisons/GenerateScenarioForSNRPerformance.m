clc
clear all
close all
% Add required paths
addpath '../Simulation' '..\Utility'
%% SIMULATION PARAMETERS
acousticParameters = AdjustSimulationParameters;
multipathMinRange = 5;
multipathMaxRange = 50;
multipathNumber = 10;
multipathMinAngleDifference = 30;
trialNum = 100;
AngleType = 'Free';
saveDirectory = '..\Scenarios\SNRPerformance\';
saveName = sprintf('%sMultiPathRange%d_%dNum%dAngleDif%dTrial-',AngleType,multipathMinRange,multipathMaxRange,multipathNumber,multipathMinAngleDifference);
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
    multipathRanges = [0,rand(1,multipathNumber)*(multipathMaxRange-multipathMinRange)+multipathMinRange];
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
    [~, avsDataNoiseless, noiseOrg] = SimulateAvsDataOutput(acousticParameters, signalOrj, len, AzimuthActual, ElevationActual, fs, multipathRanges);
    save([saveDirectory,saveName,sprintf('%d',trIndx)],'acousticParameters','noiseOrg','avsDataNoiseless','fs','AzimuthActual','ElevationActual','multipathRanges');
    disp(sprintf('Trial %d / %d was saved successfully\n',trIndx,trialNum));
end
