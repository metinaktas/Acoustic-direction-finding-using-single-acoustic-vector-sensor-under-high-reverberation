function varargout = SimulateAvsDataOutput(parameters, transmitSignal, receivedLength, azimuthDeg, elevationDeg, fs, multipathRanges, varargin)
if nargin > 7
    frqBinNumber = varargin{1};
else
    frqBinNumber = [];
end
%% INITIALIZATIONS
azimuthRad = azimuthDeg / 180 * pi;
elevationRad = elevationDeg / 180 * pi;
%% Filter the transmitted signal
if parameters.LPFilter.FrqStop > fs/2
    parameters.LPFilter.FrqStop = fs/2;
    parameters.LPFilter.FrqPass = min(parameters.LPFilter.FrqStop - 2000, 6000);
end
Hd = LowpassFilter(fs, parameters.LPFilter.FrqPass, parameters.LPFilter.FrqStop);
transmitSignal = conv(Hd.Numerator, transmitSignal);
%% Generate received signal at the position of acoustic sensors
vs = parameters.Environment.SoundSpeed;
timeDifSampleNum = (parameters.Environment.Distance + (rand-0.5)/50) / vs * fs;
timeDifSampleNumInt = round(timeDifSampleNum);
pastSampleNum = receivedLength - length(transmitSignal) - timeDifSampleNumInt;
t = [0:length(transmitSignal)-1];
cs = spline(t,transmitSignal);
tinterp = [0:length(transmitSignal)-1] + (timeDifSampleNum - timeDifSampleNumInt);
transmitSignalInterp = ppval(cs,tinterp);
if pastSampleNum <= 0
    lastIndx = length(transmitSignal) + pastSampleNum;
    undistortedReceivedSignal = [zeros(1,timeDifSampleNumInt), transmitSignalInterp(1:lastIndx)];
else
    undistortedReceivedSignal = [zeros(1,timeDifSampleNumInt), transmitSignalInterp, zeros(1,pastSampleNum)];
end
rangeDif = multipathRanges - multipathRanges(1);
range = parameters.Environment.Distance + rangeDif;
timeShiftMultipath = rangeDif / vs;
sampleShiftMultipath = timeShiftMultipath * fs;
intSampleShiftMultipath = round(timeShiftMultipath * fs);
receivedSignal = zeros(length(multipathRanges),length(undistortedReceivedSignal) + length(parameters.Environment.RoomImp)-1+max(intSampleShiftMultipath));
filteredSignal = conv(undistortedReceivedSignal, parameters.Environment.RoomImp) * exp(-parameters.Environment.Attenuation * parameters.Environment.Distance);
receivedSignal(1,1:length(filteredSignal)) = filteredSignal;
for i = 2:length(multipathRanges)
    subSampleShift = sampleShiftMultipath(i) - intSampleShiftMultipath(i);
    cs = spline([-1, [0:length(filteredSignal)]],[0,filteredSignal,0]);
    tinterp = [0:length(filteredSignal)-1] + subSampleShift;
    interpedSignal = ppval(cs,tinterp);
    receivedSignal(i,1:length(interpedSignal)+intSampleShiftMultipath(i)) = [zeros(1,intSampleShiftMultipath(i)),interpedSignal * exp(-parameters.Environment.Attenuation * range(i))];
end
clear i
%% Generate the AVS array response
avsDataNoiseless = 0;
for i = 1:length(azimuthRad)
    % Construct gain distortion matrix
    gainError = diag([parameters.Gain.p, parameters.Gain.x, parameters.Gain.y, parameters.Gain.z]);
    indxAzimuth = find(parameters.RandValueAzimuthList <= azimuthDeg(i));indxAzimuth = indxAzimuth(end);
    indxElevation = find(parameters.RandValueElevationList <= elevationDeg(i));indxElevation = indxElevation(end);
    randValue = zeros(1,4);
    for randIndx = 1:4
        randValue(randIndx) = parameters.RandGainValue(indxAzimuth,indxElevation,randIndx);
    end
    GainDistortionMatrix = gainError*diag(randValue) + eye(4);
    % Construct orientation distortion matrix
    orientErrAzXRad = parameters.Orientation.x.Azimuth / 180 * pi;
    orientErrAzYRad = parameters.Orientation.y.Azimuth / 180 * pi;
    orientErrAzZRad = parameters.Orientation.z.Azimuth / 180 * pi;
    orientErrElXRad = parameters.Orientation.x.Elevation / 180 * pi;
    orientErrElYRad = parameters.Orientation.y.Elevation / 180 * pi;
    orientErrElZRad = parameters.Orientation.z.Elevation / 180 * pi;
    orientDistortionMatrix = diag([1,...
        cos(orientErrAzXRad) * cos(orientErrElXRad),...
        cos(orientErrAzYRad) * cos(orientErrElYRad),...
        cos(orientErrAzZRad) * cos(orientErrElZRad)]);
    % Construct array manifold error vector
    manifoldErrorVector = [0;...
        -cos(elevationRad(i))*sin(azimuthRad(i))*cos(orientErrElXRad)*sin(orientErrAzXRad) + sin(elevationRad(i))*sin(orientErrElXRad);...
        cos(elevationRad(i))*cos(azimuthRad(i))*cos(orientErrElYRad)*sin(orientErrAzYRad) - sin(elevationRad(i))*sin(orientErrElYRad);....
        -cos(elevationRad(i))*cos(azimuthRad(i))*cos(orientErrElZRad)*sin(orientErrAzZRad) + cos(elevationRad(i))*sin(azimuthRad(i))*sin(orientErrElZRad)];
    % Generate noisy sound sequence
    idealArrayResponse = [1;...                    % Pressure
        cos(elevationRad(i)) * cos(azimuthRad(i));...          % Velocity x
        cos(elevationRad(i)) * sin(azimuthRad(i));...          % Velocity y
        sin(elevationRad(i))];                           % Velocity z
    % Displacement error in meter
    unitVector = [cos(elevationRad(i)) * cos(azimuthRad(i)); cos(elevationRad(i)) * sin(azimuthRad(i)); sin(elevationRad(i))];
    azRad = parameters.Displacement.p.Azimuth/180*pi;
    elRad = parameters.Displacement.p.Elevation/180*pi;
    displacementErrP = parameters.Displacement.p.Range/1000 * [cos(azRad)*cos(elRad);sin(azRad)*cos(elRad);sin(elRad)].' * unitVector;
    azRad = parameters.Displacement.x.Azimuth/180*pi;
    elRad = parameters.Displacement.x.Elevation/180*pi;
    displacementErrX = parameters.Displacement.x.Range/1000 * [cos(azRad)*cos(elRad);sin(azRad)*cos(elRad);sin(elRad)].' * unitVector;
    azRad = parameters.Displacement.y.Azimuth/180*pi;
    elRad = parameters.Displacement.y.Elevation/180*pi;
    displacementErrY = parameters.Displacement.y.Range/1000 * [cos(azRad)*cos(elRad);sin(azRad)*cos(elRad);sin(elRad)].' * unitVector;
    azRad = parameters.Displacement.z.Azimuth/180*pi;
    elRad = parameters.Displacement.z.Elevation/180*pi;
    displacementErrZ = parameters.Displacement.z.Range/1000 * [cos(azRad)*cos(elRad);sin(azRad)*cos(elRad);sin(elRad)].' * unitVector;
    displacementErr = [displacementErrP;displacementErrX;displacementErrY;displacementErrZ];
    timeShift = displacementErr / vs;
    sampleShift = timeShift * fs;
    sampleShiftInteger = round(sampleShift);
    sampleShiftRemaining = sampleShift - sampleShiftInteger;
    arrayGainResponse = GainDistortionMatrix * (orientDistortionMatrix * idealArrayResponse + manifoldErrorVector);
    if isempty(frqBinNumber) == false
        arrayPhaseResponse = exp(-sqrt(-1)*2*pi*timeShift*linspace(0,fs/2,frqBinNumber+1));
        arrayResponse = diag(arrayGainResponse) * arrayPhaseResponse;
    else
        arrayResponse = [];
    end
    %% Generate AVS data
    % Interpolated and time shifted received signal
    tinterp = zeros(length(sampleShiftRemaining),length(receivedSignal(i,:)));
    torg = [0:1:length(receivedSignal(i,:))-1]/fs;
    signalAtSensors = zeros(length(sampleShiftRemaining),length(receivedSignal(i,:)));
    for s = 1:length(sampleShift)
        cs = spline([torg(1) - max(abs(timeShift(s)),1/fs), torg, torg(end)+max(abs(timeShift(s)),1/fs)],[0,receivedSignal(i,:),0]);
        tinterp(s,:) = torg + timeShift(s);
        signalAtSensors(s,:) = ppval(cs,tinterp(s,:)).';
    end
    % Generate avs data
    avsDataNoiseless = avsDataNoiseless + diag(arrayGainResponse) * signalAtSensors;
end
% Generate noise sequence for each sensor in AVS
signalPower = mean(diag(receivedSignal * receivedSignal' / size(receivedSignal,2)));
noise = randn(4,length(receivedSignal(i,:)));
noisePower = diag(noise * noise' / size(noise,2));
snrCorrection = sqrt(10^(-parameters.Environment.SNR/10) * signalPower ./ noisePower);
noise = diag(snrCorrection) * noise;
% Generate noisy avs data
avsData = avsDataNoiseless + noise;
%% OUTPUT
varargout{1} = avsData;
varargout{2} = avsDataNoiseless;
varargout{3} = noise;