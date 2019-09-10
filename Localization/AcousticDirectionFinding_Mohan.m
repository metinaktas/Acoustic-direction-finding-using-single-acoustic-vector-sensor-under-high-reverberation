% The function for determining acoustic source direction

% Metin AKTAÞ
% 29.04.2013

function varargout = AcousticDirectionFinding_Mohan(avsDataInTime, fs, calibrationMap)
%% DETERMINE PARAMETERS
% The time length for covariance matrix estimation (s)
timeLengthForCovariance = 30e-3;
% The time length for DOA search
timeLengthForDOASearch = max(timeLengthForCovariance,1);
% The sample length for covariance matrix estimation (s)
sampleLengthForCovariance = 30;
% The short-time fourier transform overlapping ratio (0 < ratio < 1)
stftOverlapRatio = 0.25;
% The short-time fourier transform rectengular window size (ms)
stftWindowSize = 128;%256;
% The short-time fourier transform overlapping size (ms)
stftOverlapSize = min(round(stftWindowSize*stftOverlapRatio),stftWindowSize-1);
% The number of frequency points in short time fourier transform
if isempty(calibrationMap) == true
    stftFrqNumber = 256;%512;
else
    stftFrqNumber = (size(calibrationMap.SteeringVector,3)-1)*2;
end
% Magnitude Squared Coherences threshold for coherence test
mscThreshold = 0.9;
% Phase threshold for coherence test
phaseThreshold = 0.01; % Radian
%% TRANSFORMATION TO TIME-FREQUENCY
avsDataInFrq = struct('p',[],'vx',[], 'vy',[], 'vz',[]);
frqVector = struct('p',[],'vx',[], 'vy',[], 'vz',[]);
timeVector = struct('p',[],'vx',[], 'vy',[], 'vz',[]);
% sampleLength = min([length(avsDataInTime.p),ceil(timeLengthForDOASearch * fs)]);
sampleLength = length(avsDataInTime.p);
win = ones(stftWindowSize,1);
[avsDataInFrq.p,frqVector.p,timeVector.p] = stft(avsDataInTime.p(1:sampleLength),win,stftOverlapSize,stftFrqNumber,fs);
[avsDataInFrq.vx,frqVector.vx,timeVector.vx] = stft(avsDataInTime.vx(1:sampleLength),win,stftOverlapSize,stftFrqNumber,fs);
[avsDataInFrq.vy,frqVector.vy,timeVector.vy] = stft(avsDataInTime.vy(1:sampleLength),win,stftOverlapSize,stftFrqNumber,fs);
[avsDataInFrq.vz,frqVector.vz,timeVector.vz] = stft(avsDataInTime.vz(1:sampleLength),win,stftOverlapSize,stftFrqNumber,fs);
frqLen = round(stftFrqNumber/2)+1;
timeLen = length(timeVector.p);
%% COHERENCE TEST
% Optain covariance matrices for each time-frequenc bin and test them for
% the coherency
[coherenceTest, covarianceMatrices] = CoherenceTest_Mohan(avsDataInFrq, sampleLengthForCovariance, mscThreshold, phaseThreshold);
%% Apply coherence test result
avsDataInFrq.p = avsDataInFrq.p .* coherenceTest;
avsDataInFrq.vx = avsDataInFrq.vx .* coherenceTest;
avsDataInFrq.vy = avsDataInFrq.vy .* coherenceTest;
avsDataInFrq.vz = avsDataInFrq.vz .* coherenceTest;
%% DIRECTIONAL SPECTRUM ESTIMATION
% Determine the search space for the doa angles
if isempty(calibrationMap) == true
    azimuthSearchDeg = [0:1:359];
    azimuthSearchRad = azimuthSearchDeg/180*pi;
    elevationSearchDeg = [-90:1:90];
    elevationSearchRad = elevationSearchDeg/180*pi;
    AngleList = [kron(ones(1,length(elevationSearchRad)), azimuthSearchRad);kron(elevationSearchRad,ones(1,length(azimuthSearchRad)))];
    arrayResponsePerBin = [ones(1,size(AngleList,2));...                    % Pressure
                           cos(AngleList(2,:)) .* cos(AngleList(1,:));...          % Velocity x
                           cos(AngleList(2,:)) .* sin(AngleList(1,:));...          % Velocity y
                           sin(AngleList(2,:))];
else
    AngleList = calibrationMap.AngleList;
end
% Estimate the directional spectrum for all time-frequency bin that
% satisfies the coherence test
i = kron(ones(1,timeLen),[1:frqLen]);
j = kron([1:timeLen],ones(1,frqLen));
% BLOCK COMPUTATION IS ERRENEOUS
[V,eigenVectMatrix,eigVal] = EigenvalueDecomposition4x4Block(covarianceMatrices,frqLen,timeLen);

indx = kron([0:16:(frqLen*timeLen-1)*16].',ones(3,1)) + kron(ones(frqLen*timeLen,1),[5;9;13]);
nullVectH = conj([eigenVectMatrix(indx),eigenVectMatrix(indx+1),eigenVectMatrix(indx+2),eigenVectMatrix(indx+3)]);
% SpectraFullMatrix = abs((nullVectH) * arrayResponsePerBin).^2;
cnt = 0;
cntq = 0;
cntSvd = 0;
SpectraFull = zeros(size(arrayResponsePerBin,2),frqLen*timeLen);
for timeIndx = 1:timeLen
    for frqIndx = 1:frqLen
        if (coherenceTest(frqIndx,timeIndx) == 1)
            SpectraFullMatrix = abs(nullVectH(cnt+1:cnt+3,:) * arrayResponsePerBin).^2;
            R = covarianceMatrices(:,cntSvd+1:cntSvd+4);
            [U,S,V] = svd(R);
            SpectraFull(:,cntq+1) = 1./sum(abs(V(:,2:end)' * arrayResponsePerBin).^2);
        end
        cnt = cnt + 3;
        cntq = cntq + 1;
        cntSvd = cntSvd + 4;
    end
end

MapStr = struct('SpectraFull',SpectraFull, 'frqLen',frqLen, 'timeLen',timeLen, 'azimuthSearchRad', azimuthSearchRad, 'elevationSearchRad', elevationSearchRad);
threshold = sum(sum(SpectraFull))/sum(sum(coherenceTest))/10;
DirectionalTimeFrqMap = ExtractDirectionalTimeFrequencyMap_Mohan(MapStr, threshold, 10);
if isempty(DirectionalTimeFrqMap.MeanAngle) == true
    DirectionalTimeFrqMap = ExtractDirectionalTimeFrequencyMap_Mohan(MapStr, 0, 10);
end
if isempty(DirectionalTimeFrqMap.MeanAngle) == false
    Azimuth = DirectionalTimeFrqMap.MeanAngle(1,:);
    Elevation = DirectionalTimeFrqMap.MeanAngle(2,:);
else
    Azimuth = NaN;
    Elevation = NaN;
end
%% OUTPUT
varargout{1} = Azimuth;
varargout{2} = Elevation;
varargout{3} = DirectionalTimeFrqMap;