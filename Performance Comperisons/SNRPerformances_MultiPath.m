clc
clear all
close all
% Add required paths
addpath '..\Localization' '..\Utility'
%% SIMULATION PARAMETERS
scenarioDirectory = '..\Scenarios\SNRPerformance\';
scenarioName = 'FreeMultiPathRange5_50Num10AngleDif30Trial-';
saveResultDirectory = '..\Results\SNRPerformance\';
saveResultName = ['Results_SNR_Metin_',scenarioName];saveResultName(end-5:end) = [];
trialRunList = 1:100;
trialNumRun = length(trialRunList);
%% ALGORITHM PARAMETERS
AlgParameters.DOAEst.TailThresholdScaleGlobal = 0.4;
AlgParameters.DOAEst.TailThresholdScaleLocal = 0.4;

AlgParameters.Stft.WindowSize = 128;                            % The short-time fourier transform window size (sample)
AlgParameters.Stft.Win = hann(AlgParameters.Stft.WindowSize);   % The window for short-time fourier transform
AlgParameters.Stft.OverlapRatio = 0.25;                         % The short-time fourier transform overlapping ratio (0 < ratio < 1)
AlgParameters.Stft.FrqNumber = 512;

AlgParameters.CoherenceTest.ConsistencyCheckThreshold = 0.05;
AlgParameters.CoherenceTest.ClusteringAngleThreshold = 5;
%% SNR PERFORMANCE MEASURE
SNRList = 0:5:40;
AzimuthActual_List = realmax*ones(trialNumRun,length(SNRList));
ElevationActual_List = realmax*ones(trialNumRun,length(SNRList));
AzimuthEstimated_MetinLowComplexity = realmax*ones(trialNumRun,length(SNRList));
ElevationEstimated_MetinLowComplexity = realmax*ones(trialNumRun,length(SNRList));
AzimuthEstimated_Mohan = realmax*ones(trialNumRun,length(SNRList));
ElevationEstimated_Mohan = realmax*ones(trialNumRun,length(SNRList));
totalElapsedTime_MetinLowComplexity = 0;
totalElapsedTime_Metin = 0;
totalElapsedTime_Mohan = 0;
for trIndxRun = 1:trialNumRun
    load([scenarioDirectory,scenarioName,sprintf('%d',trialRunList(trIndxRun))])
    for snrIndx = 1:length(SNRList)
        snrCorrection = sqrt(10^((acousticParameters.Environment.SNR - SNRList(snrIndx))/10));
        noise = noiseOrg * snrCorrection;
        avsData = avsDataNoiseless + noise;
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
        AzimuthActual_List(trIndxRun,snrIndx) = AzimuthActual(1);
        ElevationActual_List(trIndxRun,snrIndx) = ElevationActual(1);

        AzimuthEstimated_MetinLowComplexity(trIndxRun,snrIndx) = Azimuth_MetinLowComplexity(1);
        ElevationEstimated_MetinLowComplexity(trIndxRun,snrIndx) = Elevation_MetinLowComplexity(1);        
        AzimuthEstimated_Mohan(trIndxRun,snrIndx) = Azimuth_Mohan(1);
        ElevationEstimated_Mohan(trIndxRun,snrIndx) = Elevation_Mohan(1);        
    end
    disp(sprintf('%d / %d was finished for SNR Performance with MULTIPATH and SAFEZone Azimuth and SAFEZone Elevation',trIndxRun,trialNumRun))
end
nanListMetinLowComplexity = isnan(AzimuthEstimated_MetinLowComplexity) | isnan(ElevationEstimated_MetinLowComplexity);
indxMetinLowComplexity = find(nanListMetinLowComplexity == 1);
AzimuthEstimated_MetinLowComplexity(indxMetinLowComplexity) = AzimuthActual_List(indxMetinLowComplexity);
ElevationEstimated_MetinLowComplexity(indxMetinLowComplexity) = ElevationActual_List(indxMetinLowComplexity);
effectiveTrial_MetinLowComplexity = sum(not(nanListMetinLowComplexity));
MseAzimuth_MetinLowComplexityTot = zeros(1,length(SNRList));
for trIndxRun = 1:trialNumRun    
    MseAzimuth_MetinLowComplexityTot = MseAzimuth_MetinLowComplexityTot + min([abs(AzimuthEstimated_MetinLowComplexity(trIndxRun,:) - AzimuthActual_List(trIndxRun,:));360-abs(AzimuthEstimated_MetinLowComplexity(trIndxRun,:) - AzimuthActual_List(trIndxRun,:))]).^2;
end
MseAzimuth_MetinLowComplexity = sqrt(MseAzimuth_MetinLowComplexityTot./effectiveTrial_MetinLowComplexity);
MseElevation_MetinLowComplexity = sqrt(sum(abs(ElevationEstimated_MetinLowComplexity - ElevationActual_List).^2)./effectiveTrial_MetinLowComplexity);
ElapsedTime_MetinLowComplexity = totalElapsedTime_MetinLowComplexity / trialNumRun / length(SNRList);

nanListMohan = isnan(AzimuthEstimated_Mohan) | isnan(ElevationEstimated_Mohan);
indxMohan = find(nanListMohan == 1);
AzimuthEstimated_Mohan(indxMohan) = AzimuthActual_List(indxMohan);
ElevationEstimated_Mohan(indxMohan) = ElevationActual_List(indxMohan);
effectiveTrial_Mohan = sum(not(nanListMohan));
MseAzimuth_MohanTot = zeros(1,length(SNRList));
for trIndxRun = 1:trialNumRun    
    MseAzimuth_MohanTot = MseAzimuth_MohanTot + min([abs(AzimuthEstimated_Mohan(trIndxRun,:) - AzimuthActual_List(trIndxRun,:));360-abs(AzimuthEstimated_Mohan(trIndxRun,:) - AzimuthActual_List(trIndxRun,:))]).^2;
end
MseAzimuth_Mohan = sqrt(MseAzimuth_MohanTot./effectiveTrial_Mohan);
MseElevation_Mohan = sqrt(sum(abs(ElevationEstimated_Mohan - ElevationActual_List).^2)./effectiveTrial_Mohan);
ElapsedTime_Mohan = totalElapsedTime_Mohan / trialNumRun / length(SNRList);

save([saveResultDirectory,saveResultName,sprintf('_TrialNum%d',trialNumRun)])

figure(1)
hold off
semilogy(SNRList, MseAzimuth_MetinLowComplexity)
hold on
semilogy(SNRList, MseAzimuth_Mohan,'r')
legend('Metin Low Complexity','Mohan')
xlabel('SNR, dB')
ylabel('MSE, dB')
title('Azimuth Estimation')

figure(2)
hold off
semilogy(SNRList, MseElevation_MetinLowComplexity)
hold on
semilogy(SNRList, MseElevation_Mohan,'r')
legend('Metin Low Complexity','Mohan')
xlabel('SNR, dB')
ylabel('MSE, dB')
title('Elevation Estimation')