clc
clear all
close all
%% SIMULATION PARAMETERS
SNR = 20;
scenarioDirectory = '..\Scenarios\MultipathRangePerformance\';
scenarioName = 'FreeMultiPathNum1AngleDif30Trial-';
saveResultDirectory = '..\Results\MultipathRangePerformanceResults\';
saveResultName = [sprintf('Results_MultipathRange_Metin_SNR%d',SNR),scenarioName];saveResultName(end-5:end) = [];
trialRunList = 1:100;
trialNumRun = length(trialRunList);
%% ALGORITHM PARAMETERS
AlgParameters.DOAEst.TailThresholdScaleGlobal = 0.3;
AlgParameters.DOAEst.TailThresholdScaleLocal = 0.3;

AlgParameters.Stft.WindowSize = 128;                            % The short-time fourier transform window size (sample)
AlgParameters.Stft.Win = hann(AlgParameters.Stft.WindowSize);   % The window for short-time fourier transform
AlgParameters.Stft.OverlapRatio = 0.25;                         % The short-time fourier transform overlapping ratio (0 < ratio < 1)
AlgParameters.Stft.FrqNumber = 512;

AlgParameters.CoherenceTest.ConsistencyCheckThreshold = 0.05;
AlgParameters.CoherenceTest.ClusteringAngleThreshold = 10;
%% MULTIPATH RANGE PERFORMANCE MEASURE
load([scenarioDirectory,scenarioName,sprintf('%d',1)])
AzimuthActual_List = realmax*ones(trialNumRun,length(multipathRangesList));
ElevationActual_List = realmax*ones(trialNumRun,length(multipathRangesList));
AzimuthEstimated_MetinLowComplexity = realmax*ones(trialNumRun,length(multipathRangesList));
ElevationEstimated_MetinLowComplexity = realmax*ones(trialNumRun,length(multipathRangesList));
AzimuthEstimated_Mohan = realmax*ones(trialNumRun,length(multipathRangesList));
ElevationEstimated_Mohan = realmax*ones(trialNumRun,length(multipathRangesList));
totalElapsedTime_MetinLowComplexity = 0;
totalElapsedTime_Mohan = 0;
for trIndxRun = 1:trialNumRun
    load([scenarioDirectory,scenarioName,sprintf('%d',trialRunList(trIndxRun))])
    snrCorrection = sqrt(10^((acousticParameters.Environment.SNR - SNR)/10));
    for rIndx = 1:length(multipathRangesList)
        noise = noiseOrg{rIndx} * snrCorrection;
        avsData = avsDataNoiseless{rIndx} + noise;
        avsDataInTime.p = avsData(1,:);
        avsDataInTime.vx = avsData(2,:);
        avsDataInTime.vy = avsData(3,:);
        avsDataInTime.vz = avsData(4,:);
        %% DIRECTION FINDIG   
        startTime = tic;
        [Azimuth_MetinLowComplexity, Elevation_MetinLowComplexity] = AcousticDirectionFinding_MetinLowComplexity(avsDataInTime, fs, AlgParameters);
        totalElapsedTime_MetinLowComplexity = totalElapsedTime_MetinLowComplexity + toc(startTime);

        startTime = tic;
        [Azimuth_Mohan, Elevation_Mohan, DirectionalTimeFrqMap_Mohan] = AcousticDirectionFinding_Mohan(avsDataInTime, fs, []);
        totalElapsedTime_Mohan = totalElapsedTime_Mohan + toc(startTime);
        %% RESULTS
        AzimuthActual_List(trIndxRun,rIndx) = AzimuthActual(1);
        ElevationActual_List(trIndxRun,rIndx) = ElevationActual(1);

        AzimuthEstimated_MetinLowComplexity(trIndxRun,rIndx) = Azimuth_MetinLowComplexity(1);
        ElevationEstimated_MetinLowComplexity(trIndxRun,rIndx) = Elevation_MetinLowComplexity(1);        
        AzimuthEstimated_Mohan(trIndxRun,rIndx) = Azimuth_Mohan(1);
        ElevationEstimated_Mohan(trIndxRun,rIndx) = Elevation_Mohan(1);        
    end
    disp(sprintf('%d / %d was finished for Multipath Range Performance with MULTIPATH and SAFEZone Azimuth and SAFEZone Elevation',trIndxRun,trialNumRun))
end
nanListMetinLowComplexity = isnan(AzimuthEstimated_MetinLowComplexity) | isnan(ElevationEstimated_MetinLowComplexity);
indxMetinLowComplexity = find(nanListMetinLowComplexity == 1);
AzimuthEstimated_MetinLowComplexity(indxMetinLowComplexity) = AzimuthActual_List(indxMetinLowComplexity);
ElevationEstimated_MetinLowComplexity(indxMetinLowComplexity) = ElevationActual_List(indxMetinLowComplexity);
effectiveTrial_MetinLowComplexity = sum(not(nanListMetinLowComplexity));
MseAzimuth_MetinLowComplexityTot = zeros(1,length(multipathRangesList));
for trIndxRun = 1:trialNumRun    
    MseAzimuth_MetinLowComplexityTot = MseAzimuth_MetinLowComplexityTot + min([abs(AzimuthEstimated_MetinLowComplexity(trIndxRun,:) - AzimuthActual_List(trIndxRun,:));360-abs(AzimuthEstimated_MetinLowComplexity(trIndxRun,:) - AzimuthActual_List(trIndxRun,:))]).^2;
end
MseAzimuth_MetinLowComplexity = sqrt(MseAzimuth_MetinLowComplexityTot./effectiveTrial_MetinLowComplexity);
MseElevation_MetinLowComplexity = sqrt(sum(abs(ElevationEstimated_MetinLowComplexity - ElevationActual_List).^2)./effectiveTrial_MetinLowComplexity);
ElapsedTime_MetinLowComplexity = totalElapsedTime_MetinLowComplexity / trialNumRun / length(multipathRangesList);

nanListMohan = isnan(AzimuthEstimated_Mohan) | isnan(ElevationEstimated_Mohan);
indxMohan = find(nanListMohan == 1);
AzimuthEstimated_Mohan(indxMohan) = AzimuthActual_List(indxMohan);
ElevationEstimated_Mohan(indxMohan) = ElevationActual_List(indxMohan);
effectiveTrial_Mohan = sum(not(nanListMohan));
MseAzimuth_MohanTot = zeros(1,length(multipathRangesList));
for trIndxRun = 1:trialNumRun    
    MseAzimuth_MohanTot = MseAzimuth_MohanTot + min([abs(AzimuthEstimated_Mohan(trIndxRun,:) - AzimuthActual_List(trIndxRun,:));360-abs(AzimuthEstimated_Mohan(trIndxRun,:) - AzimuthActual_List(trIndxRun,:))]).^2;
end
MseAzimuth_Mohan = sqrt(MseAzimuth_MohanTot./effectiveTrial_Mohan);
MseElevation_Mohan = sqrt(sum(abs(ElevationEstimated_Mohan - ElevationActual_List).^2)./effectiveTrial_Mohan);
ElapsedTime_Mohan = totalElapsedTime_Mohan / trialNumRun / length(multipathRangesList);

save([saveResultDirectory,saveResultName,sprintf('_TrialNum%d',trialNumRun)])

figure(3)
hold off
semilogy(multipathRangesList, MseAzimuth_MetinLowComplexity)
hold on
semilogy(multipathRangesList, MseAzimuth_Mohan,'r')
legend('Metin Low Complexity','Mohan')
xlabel('Multipath Range, meters')
ylabel('MSE, dB')
title('Azimuth Estimation')

figure(4)
hold off
semilogy(multipathRangesList, MseElevation_MetinLowComplexity)
hold on
semilogy(multipathRangesList, MseElevation_Mohan,'r')
legend('Metin Low Complexity','Mohan')
xlabel('Multipath Range, meters')
ylabel('MSE, dB')
title('Elevation Estimation')