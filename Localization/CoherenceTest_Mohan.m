% The function for optaining covariance matrices for each time-frequenc bin
% and test them for the coherency

function varargout = CoherenceTest_Mohan(avsDataInFrq, sampleLengthForCovariance, mscThreshold, phaseThreshold)

frqLen = round(size(avsDataInFrq.p,1)/2)+1;
timeLen = size(avsDataInFrq.p,2);
covarianceMatrices = zeros(4,4*timeLen*frqLen);
result = ones(frqLen,timeLen);
timeStart = ones(frqLen,1);
testResult = ones(frqLen,1);
scaleForR = ones(frqLen,1);
cnt = 0;
for covTimeIndx = 1:timeLen
    timeScale = max(1,covTimeIndx-sampleLengthForCovariance+1):covTimeIndx;
    timeScaleMat = ones(frqLen,length(timeScale));
    if size(timeScaleMat,2) > 1
        scaleForR = sum(timeScaleMat.').';
    end
    timeFrqBin = zeros(frqLen*4,length(timeScale));
    timeFrqBin(1:4:end,:) = avsDataInFrq.p(1:frqLen,timeScale) .* timeScaleMat ./ sqrt(scaleForR*ones(1,size(timeScaleMat,2)));
    timeFrqBin(2:4:end,:) = avsDataInFrq.vx(1:frqLen,timeScale) .* timeScaleMat ./ sqrt(scaleForR*ones(1,size(timeScaleMat,2)));
    timeFrqBin(3:4:end,:) = avsDataInFrq.vy(1:frqLen,timeScale) .* timeScaleMat ./ sqrt(scaleForR*ones(1,size(timeScaleMat,2)));
    timeFrqBin(4:4:end,:) = avsDataInFrq.vz(1:frqLen,timeScale) .* timeScaleMat ./ sqrt(scaleForR*ones(1,size(timeScaleMat,2)));
    R = zeros(4,frqLen*4);
    cntR = 0;
    for covFrqIndx = 1:frqLen
        R(:,cntR+1:cntR+4) = timeFrqBin(cntR+1:cntR+4,:) * timeFrqBin(cntR+1:cntR+4,:)';
        cntR = cntR + 4;
    end    
    for covFrqIndx = 1:frqLen
        tempR = R(:,(covFrqIndx-1)*4+1:covFrqIndx*4);
        covarianceMatrices(:,cnt+1:cnt+4) = tempR;
        cnt = cnt + 4;
    end
    diagR = zeros(1,size(R,2));
    diagRMat = zeros(4,frqLen);
    for i = 1:4
        diagR(i:4:end) = R(i,i:4:end);
        diagRMat(i,:) = R(i,i:4:end);
    end
    mscIndx = [1,1,1,2,2,3;2,3,4,3,4,4];
    phaseIndx = [1,1,1,2;2,2,3,3;3,4,4,4];
    offsetIndicesForMsc = kron(0:16:prod(size(R))-16,ones(1,size(mscIndx,2)));
    i = mscIndx(1,:);
    j = mscIndx(2,:);
    indxij =  kron(ones(1,frqLen),i + (j-1)*size(R,1)) + offsetIndicesForMsc;
    indxii =  kron(ones(1,frqLen),i + (i-1)*size(R,1)) + offsetIndicesForMsc;
    indxjj =  kron(ones(1,frqLen),j + (j-1)*size(R,1)) + offsetIndicesForMsc;
    mscVect = (abs(R(indxij)).^2 ./ (R(indxii).*R(indxjj)) >= mscThreshold);
    mscMat = zeros(size(mscIndx,2),frqLen);
    for i = 1:size(mscIndx,2)
        mscMat(i,:) = mscVect(i:size(mscIndx,2):end);
    end
    mscTest = prod(mscMat);
    
    offsetIndicesForPhase = kron(0:16:prod(size(R))-16,ones(1,size(phaseIndx,2)));
    i = phaseIndx(1,:);
    j = phaseIndx(2,:);
    k = phaseIndx(3,:);
    indxij =  kron(ones(1,frqLen),i + (j-1)*size(R,1)) + offsetIndicesForPhase;
    indxik =  kron(ones(1,frqLen),i + (k-1)*size(R,1)) + offsetIndicesForPhase;
    indxjk =  kron(ones(1,frqLen),j + (k-1)*size(R,1)) + offsetIndicesForPhase;    
    phaseVect = (mod(abs(angle(R(indxij)) - angle(R(indxik)) + angle(R(indxjk)))+phaseThreshold,2*pi) <= (2*phaseThreshold));
    phaseMat = zeros(size(phaseIndx,2),frqLen);
    for i = 1:size(phaseIndx,2)
        phaseMat(i,:) = phaseVect(i:size(phaseIndx,2):end);
    end
    phaseTest = prod(phaseMat);

    testResult = prod([mscTest;phaseTest]);
    
    result(:,covTimeIndx) = testResult.';
end
result = [result;flipud(result(2:end-1,:))];
%% OUTPUTS
varargout{1} = result;
varargout{2} = covarianceMatrices;