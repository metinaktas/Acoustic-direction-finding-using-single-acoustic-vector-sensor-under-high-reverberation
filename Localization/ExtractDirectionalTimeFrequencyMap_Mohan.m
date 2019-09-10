function DirectionalTimeFrqMap = ExtractDirectionalTimeFrequencyMap_Mohan(MapStr, threshold, angleResolution)
% Initialize parameters
azimuthSearchDeg = MapStr.azimuthSearchRad / pi * 180;
elevationSearchDeg = MapStr.elevationSearchRad / pi * 180;
AngleList = [];
AngleCntList = [];
AngleTotList = [];
SpectraValList = [];
MapValList = [];
MapValIncrement = 10;
AngleAz = [kron(ones(1,length(elevationSearchDeg)), azimuthSearchDeg);kron(elevationSearchDeg,ones(1,length(azimuthSearchDeg)))];
% Initialize map
map = zeros(MapStr.frqLen, MapStr.timeLen);
% Analyze each time-frequency bin
cnt = 1;
for timeIndx = 1:MapStr.timeLen
    for frqIndx = 1:MapStr.frqLen
        S = MapStr.SpectraFull(:,cnt);
        if max(S) > threshold
            [Angles, SpectraVal] = PeakSearch2DVectorized(S, 0.01, azimuthSearchDeg, elevationSearchDeg, AngleAz);
            [maxSpectraVal,indxAng] = max(SpectraVal);
            if isempty(Angles) == false
                Azimuth = Angles(1,indxAng);
                Elevation = Angles(2,indxAng);
                if isempty(AngleList) == true
                    MapValList = MapValIncrement;
                    AngleList = [Azimuth;Elevation];
                    map(frqIndx,timeIndx) = MapValIncrement;
                    SpectraValList = maxSpectraVal;
                    AngleCntList = 1;
                    AngleTotList = AngleList;
                else
                    [val,indx] = min(mod(sqrt(sum(abs(AngleList - [Azimuth;Elevation]*ones(1,size(AngleList,2))).^2)),360));
                    if val > angleResolution
                        MapValList = [MapValList,MapValList(end) + MapValIncrement];
                        AngleList = [AngleList, [Azimuth;Elevation]];
                        map(frqIndx,timeIndx) = MapValList(end);
                        SpectraValList = [SpectraValList, maxSpectraVal];
                        AngleCntList = [AngleCntList, 1];
                        AngleTotList = [AngleTotList, [Azimuth;Elevation]];
                    else
                        map(frqIndx,timeIndx) = MapValList(indx);
                        SpectraValList(indx) = SpectraValList(indx) + maxSpectraVal;
                        AngleCntList(indx) = AngleCntList(indx) + 1;
                        AngleTotList(:,indx) = AngleTotList(:,indx) + [Azimuth;Elevation];
                    end
                end
            end
        end
        cnt = cnt + 1;
    end
end
indxDelete = find(AngleCntList < max(AngleCntList)*0.1);
AngleCntList(indxDelete) = [];
AngleTotList(:,indxDelete) = [];
for i = 1:length(indxDelete)
    map(find(map == MapValList(indxDelete(i)))) = 0;
end
MapValList(indxDelete) = [];
SpectraValList(indxDelete) = [];
%% OUTPUTS
if isempty(AngleTotList) == false
    MeanAngle = AngleTotList ./ (ones(2,1)*AngleCntList);
else
    MeanAngle = [];
end
try
    DirectionalTimeFrqMap = struct('map',map, 'MeanAngle', MeanAngle, 'MapVal',MapValList, 'SpectraVal',SpectraValList, 'AngleCnt',AngleCntList);
catch
    AngleTotList
end