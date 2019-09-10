% The function for short time fourier transform
function varargout = stft(dataInTime, win, stftOverlapSize, stftFrqNumber, fs, varargin)

if nargin > 5
    Delta = [0,varargin{1}];
end
stftWindowSize = length(win);
dataInTime = dataInTime(:);
sampleLength = length(dataInTime);

dataInFrq = zeros(stftFrqNumber, floor((sampleLength  - stftWindowSize) / (stftWindowSize - stftOverlapSize))+1);
indx = 1;
for i = 1:size(dataInFrq,2)
    x = dataInTime(indx:indx + stftWindowSize - 1) .* win;
    indx = indx + stftWindowSize - stftOverlapSize;
    xf = fft(x,stftFrqNumber);
    if size(xf,2) > 1
        dataInFrq(:,i) = sum((xf .* exp(-sqrt(-1)*[0:stftFrqNumber-1].' * Delta)).').';
    else
        dataInFrq(:,i) = xf;
    end
end
%% OUTPUTS
varargout{1} = dataInFrq;
varargout{2} = linspace(-fs/2,fs/2,size(dataInFrq,1));
varargout{3} = (round(stftWindowSize/2) + [0:size(dataInFrq,2)-1]*(stftWindowSize - stftOverlapSize))/fs;