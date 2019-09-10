function acousticParameters = AdjustSimulationParameters()

acousticParameters.Gain.p = 0;%0.1;
acousticParameters.Gain.x = 0;%0.1;
acousticParameters.Gain.y = 0;%0.1;
acousticParameters.Gain.z = 0;%0.1;

acousticParameters.Displacement.p.Range = 0;
acousticParameters.Displacement.p.Azimuth = 0;
acousticParameters.Displacement.p.Elevation = 0;

acousticParameters.Displacement.x.Range = 0;%5;
acousticParameters.Displacement.x.Azimuth = 0;%30;
acousticParameters.Displacement.x.Elevation = 0;

acousticParameters.Displacement.y.Range = 0;%5;
acousticParameters.Displacement.y.Azimuth = 0;%150;
acousticParameters.Displacement.y.Elevation = 0;

acousticParameters.Displacement.z.Range = 0;%10;
acousticParameters.Displacement.z.Azimuth = 0;%45;
acousticParameters.Displacement.z.Elevation = 0;%-45;

acousticParameters.Orientation.x.Azimuth = 0;%0.1;
acousticParameters.Orientation.x.Elevation = 0;%0.1;

acousticParameters.Orientation.y.Azimuth = 0;%-0.1;
acousticParameters.Orientation.y.Elevation = 0;%0.1;

acousticParameters.Orientation.z.Azimuth = 0;%0.1;
acousticParameters.Orientation.z.Elevation = 0;%-0.1;

acousticParameters.Environment.SNR = 30;
acousticParameters.Environment.Distance = 10;
acousticParameters.Environment.RoomImp = 1;
acousticParameters.Environment.Attenuation = 0.0448;
acousticParameters.Environment.SoundSpeed = 330;

acousticParameters.RandValueAzimuthList = 0:1:359;
acousticParameters.RandValueElevationList = -90:1:90;
lenAz = length(acousticParameters.RandValueAzimuthList);
lenEl = length(acousticParameters.RandValueElevationList);
lowPass = ones([lenAz,lenEl]);
acousticParameters.RandGainValue = zeros(lenAz,lenEl,4);
shiftAz = round(size(lowPass,1)*3/2);
shiftEl = round(size(lowPass,2)*3/2);
for i = 1:4
    unfilteredNoise = (rand(lenAz+size(lowPass,1),lenEl+size(lowPass,2)) - 0.5)*2;
    filteredNoise = conv2(unfilteredNoise,lowPass);
    filteredNoise = filteredNoise(shiftAz:lenAz+shiftAz-1,shiftEl:lenEl+shiftEl-1);
    filteredNoise = filteredNoise - min(min(filteredNoise));
    filteredNoise = filteredNoise / max(max(filteredNoise)) * 2 - 1;
    acousticParameters.RandGainValue(:,:,i) = filteredNoise;
end

acousticParameters.LPFilter.FrqPass = 16000;
acousticParameters.LPFilter.FrqStop = 18000;
