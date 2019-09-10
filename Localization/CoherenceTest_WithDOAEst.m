% The function for optaining covariance matrices for each time-frequenc bin
% and test them for the coherency and estimate DOA for each time-frequency
% bin that satisfies single source condition

function varargout = CoherenceTest_WithDOAEst(avsDataInFrq, parameters)

%% ALGORITHM
frqLen = round(size(avsDataInFrq.p,1)/2)+1;
timeLen = size(avsDataInFrq.p,2);
result = zeros(frqLen,timeLen);
azimuthEstimation = NaN*zeros(frqLen,timeLen);
elevationEstimation = NaN*zeros(frqLen,timeLen);
eneryMap = NaN*zeros(frqLen,timeLen);
errorMap = NaN*zeros(frqLen,timeLen);
%% Analize each time-frequency bin
for covTimeIndx = 1:timeLen
    for covFrqIndx = 1:frqLen
        %% Data in time-frequency bin (covFrqIndx, covTimeIndx)
        timeFrqBin = zeros(4,1);
        timeFrqBin(1) = avsDataInFrq.p(covFrqIndx,covTimeIndx);
        timeFrqBin(2) = avsDataInFrq.vx(covFrqIndx,covTimeIndx);
        timeFrqBin(3) = avsDataInFrq.vy(covFrqIndx,covTimeIndx);
        timeFrqBin(4) = avsDataInFrq.vz(covFrqIndx,covTimeIndx);
        %% DOA Estimation
        Xscaled = timeFrqBin / timeFrqBin(1);
        az = atan(real(timeFrqBin(3)/timeFrqBin(2))) + pi*(sign(real(Xscaled(2))) < 0);
        el = atan(real(Xscaled(4)/sqrt(abs(Xscaled(2))^2 + abs(Xscaled(3))^2)));%atan2(real(Xscaled(4)),sqrt(Xscaled(2)^2 + Xscaled(3)^2))
        %% Single Source Test
        X = [1;cos(az)*cos(el);cos(el)*sin(az);sin(el)];
        if sqrt(sum(abs(Xscaled - X).^2) / 4) < parameters.ConsistencyCheckThreshold
            Azimuth = mod(az * 180 / pi,360);
            Elevation = el * 180 / pi;
            azimuthEstimation(covFrqIndx,covTimeIndx)   = Azimuth;
            elevationEstimation(covFrqIndx,covTimeIndx) = Elevation;
            result(covFrqIndx,covTimeIndx) = 1;
            eneryMap(covFrqIndx,covTimeIndx) = abs(timeFrqBin(1,1)).^2;
            errorMap(covFrqIndx,covTimeIndx) = sqrt(sum(abs(Xscaled - X).^2) / 4);
        end
    end
end
result = [result;flipud(result(2:end-1,:))];
%% OUTPUTS
varargout{1} = azimuthEstimation;
varargout{2} = elevationEstimation;
varargout{3} = eneryMap;
varargout{4} = errorMap;
