% The function for finding the peaks in 2D MUSIC spectra
% INPUTS
% S: The 2D MUSIC spectra with rows for azimuth and columns for elevation
% angles
% thresholdRatio: The ratio scaled with the maximum value on the MUSIC
% spectra
% azimuth: Azimuth angles for S
% elevation: Elevation angles for S
function varargout = PeakSearch2DVectorized(S, thresholdRatio, azimuth, elevation, AngleAz)

% Determine the threshold and reset the pixels whose values are smaller
% than the threshold
threshold = max(max(S)) * thresholdRatio;
S(S < threshold) = 0;

% AngleAz = [kron(ones(1,length(elevation)), azimuth);kron(elevation,ones(1,length(azimuth)))];
SElev = Vector2Matrix(S,length(elevation),length(azimuth),'column');
SElev = SElev(:);

[SbinaryAz, SpectraValAz] = PeakBinary(S);
[SbinaryElev, SpectraValElev] = PeakBinary(SElev);

SbinaryElev2D = Vector2Matrix(SbinaryElev,length(elevation),length(azimuth),'row').';
Sbinary = SbinaryAz.*SbinaryElev2D(:);
indx = find(Sbinary == 1);
Angles = AngleAz(:,indx);
SpectraVal = SpectraValAz(indx);

%% OUTPUTS
varargout{1} = Angles;
varargout{2} = SpectraVal;
varargout{3} = Sbinary;



function varargout = PeakBinary(S)

% Define the output parameters
Sbinary = zeros(size(S));
SpectraVal = zeros(size(S));

% Extend the MUSIC spectra to avoids missing the peaks on the edges
% For azimuth 360+[1,2,3,...], [1,2,3,...] will be used
% For azimuth 0-[1,2,3,...], [359,358,357...] will be used
% For elevation 90+[1,2,3,...], [89,88,87,...] will be used
% For elevation -90-[1,2,3,...], [-89,-88,-87,...] will be used
S = [S(end-1);S;S(2)];

% Evaluate the 2D derivative for each pixel of original 2D MUSIC spectra
% and generate the binary image showing the peaks as 1
SMat = [S(1:end-2).';S(2:end-1).';S(3:end).'];
indx = find(max(SMat) - SMat(2,:) == 0 & SMat(2,:) > 0);
Sbinary(indx) = 1;
SpectraVal(indx) = SMat(2,indx);

%% OUTPUTS
varargout{1} = Sbinary;
varargout{2} = SpectraVal;