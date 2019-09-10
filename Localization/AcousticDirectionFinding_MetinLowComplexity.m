% The function for determining acoustic source direction

% Metin AKTAÞ
% 29.04.2013

function varargout = AcousticDirectionFinding_MetinLowComplexity(avsDataInTime, fs, AlgParameters)
%% DETERMINE PARAMETERS
% The short-time fourier transform overlapping size (ms)
stftOverlapSize = min(round(AlgParameters.Stft.WindowSize*AlgParameters.Stft.OverlapRatio),AlgParameters.Stft.WindowSize-1);
%% TRANSFORMATION TO TIME-FREQUENCY
avsDataInFrq = struct('p',[],'vx',[], 'vy',[], 'vz',[]);
frqVector = struct('p',[],'vx',[], 'vy',[], 'vz',[]);
timeVector = struct('p',[],'vx',[], 'vy',[], 'vz',[]);
sampleLength = length(avsDataInTime.p);
[avsDataInFrq.p,frqVector.p,timeVector.p] = stft(avsDataInTime.p(1:sampleLength),AlgParameters.Stft.Win,stftOverlapSize,AlgParameters.Stft.FrqNumber,fs);
[avsDataInFrq.vx,frqVector.vx,timeVector.vx] = stft(avsDataInTime.vx(1:sampleLength),AlgParameters.Stft.Win,stftOverlapSize,AlgParameters.Stft.FrqNumber,fs);
[avsDataInFrq.vy,frqVector.vy,timeVector.vy] = stft(avsDataInTime.vy(1:sampleLength),AlgParameters.Stft.Win,stftOverlapSize,AlgParameters.Stft.FrqNumber,fs);
[avsDataInFrq.vz,frqVector.vz,timeVector.vz] = stft(avsDataInTime.vz(1:sampleLength),AlgParameters.Stft.Win,stftOverlapSize,AlgParameters.Stft.FrqNumber,fs);
%% COHERENCE TEST WITH DOA ESTIMATION
% Optain covariance matrices for each time-frequenc bin and test them for
% the coherency
[azimuthMap, elevationMap, energyMap, errorMap] = CoherenceTest_WithDOAEst(avsDataInFrq, AlgParameters.CoherenceTest);
%% Mapping Elimination
%% Global Cluster Elimination
clusteringMapGlobal = ClusteringMap(azimuthMap, elevationMap, AlgParameters.CoherenceTest.ClusteringAngleThreshold, 'global');
clustNum = max(max(clusteringMapGlobal));
totEliminationGlobal = 0;
if clustNum >= 1
    Nclust = zeros(1,clustNum);
    for i = 1:clustNum
        Nclust(i) = length(find(clusteringMapGlobal == i));
    end
    indxSmallClust = find(Nclust < (max(Nclust)*AlgParameters.DOAEst.TailThresholdScaleGlobal));    
    for i = 1:length(indxSmallClust)
        totEliminationGlobal = totEliminationGlobal + length(find(clusteringMapGlobal == indxSmallClust(i)));
        indxG = find(clusteringMapGlobal == indxSmallClust(i));
        azimuthMap(indxG) = NaN;
        elevationMap(indxG) = NaN;
        clusteringMapGlobal(indxG) = NaN;        
    end    
end
%% Local Cluster Elimination
clusteringMapLocal = ClusteringMap(azimuthMap, elevationMap, AlgParameters.CoherenceTest.ClusteringAngleThreshold, 'local');
totEliminationLocal = 0;
for t = 1:size(clusteringMapLocal,2)
    clustNum = max(clusteringMapLocal(:,t));
    if clustNum >= 1
        Nclust = zeros(1,clustNum);
        for i = 1:clustNum
            Nclust(i) = length(find(clusteringMapLocal(:,t) == i));
        end
        indxSmallClust = find(Nclust < (max(Nclust)*AlgParameters.DOAEst.TailThresholdScaleLocal));
        for i = 1:length(indxSmallClust)
            totEliminationLocal = totEliminationLocal + length(find(clusteringMapLocal(:,t) == indxSmallClust(i)));
            indxL = find(clusteringMapLocal(:,t) == indxSmallClust(i));
            azimuthMap(indxL,t) = NaN;
            elevationMap(indxL,t) = NaN;
            clusteringMapLocal(indxL,t) = NaN;
        end
        temp = clusteringMapLocal(find(isnan(clusteringMapLocal(:,t)) == false),t);
        if sum(temp)/length(temp) ~= temp(1)
            azimuthMap(:,t) = NaN;
            elevationMap(:,t) = NaN;
            totEliminationLocal = totEliminationLocal + length(temp);
        end
    end
end
%%
clusteringMapGlobal = ClusteringMap(azimuthMap, elevationMap, AlgParameters.CoherenceTest.ClusteringAngleThreshold, 'global');
clustNum = max(max(clusteringMapGlobal));
azResults = NaN;
elResults = NaN;
if clustNum >= 1
    overalSelectedClustIndx = [];
    for i = 1:clustNum
        selectedClustIndx = find(clusteringMapGlobal(:) == i);
        minAz = min(azimuthMap(selectedClustIndx));
        if minAz < 2*AlgParameters.CoherenceTest.ClusteringAngleThreshold
            outOfBoundIndx = find(azimuthMap(selectedClustIndx) > 360 - 2*AlgParameters.CoherenceTest.ClusteringAngleThreshold);
            if isempty(outOfBoundIndx) == false
                az = azimuthMap(selectedClustIndx);
                az(outOfBoundIndx) = az(outOfBoundIndx) - 360;
                azimuthMap(selectedClustIndx) = az;
            end
        end
        azResults = sum(azimuthMap(selectedClustIndx).*(energyMap(selectedClustIndx)./errorMap(selectedClustIndx)))/sum(energyMap(selectedClustIndx)./errorMap(selectedClustIndx));
        elResults = sum(elevationMap(selectedClustIndx).*(energyMap(selectedClustIndx)./errorMap(selectedClustIndx)))/sum(energyMap(selectedClustIndx)./errorMap(selectedClustIndx));
        if i == 1
            angInitial = [azResults;elResults];
            overalSelectedClustIndx = [overalSelectedClustIndx,selectedClustIndx];
        else
            angDif = min(sqrt(min([abs(angInitial(1,:) - azResults);360-abs(angInitial(1,:) - azResults)]).^2 + abs(angInitial(2,:) - elResults).^2));
            if angDif < 2*AlgParameters.CoherenceTest.ClusteringAngleThreshold
                overalSelectedClustIndx = [overalSelectedClustIndx;selectedClustIndx];
            end
        end
    end
    minAz = min(azimuthMap(overalSelectedClustIndx));
    if minAz < 2*AlgParameters.CoherenceTest.ClusteringAngleThreshold
        outOfBoundIndx = find(azimuthMap(overalSelectedClustIndx) > 360 - 2*AlgParameters.CoherenceTest.ClusteringAngleThreshold);
        if isempty(outOfBoundIndx) == false
            az = azimuthMap(overalSelectedClustIndx);
            az(outOfBoundIndx) = az(outOfBoundIndx) - 360;
            azimuthMap(overalSelectedClustIndx) = az;
        end
    end
    azResults = sum(azimuthMap(overalSelectedClustIndx).*(energyMap(overalSelectedClustIndx)./errorMap(overalSelectedClustIndx)))/sum(energyMap(overalSelectedClustIndx)./errorMap(overalSelectedClustIndx));
    elResults = sum(elevationMap(overalSelectedClustIndx).*(energyMap(overalSelectedClustIndx)./errorMap(overalSelectedClustIndx)))/sum(energyMap(overalSelectedClustIndx)./errorMap(overalSelectedClustIndx));
end
%% OUTPUTS
varargout{1} = azResults(:).';
varargout{2} = elResults(:).';
