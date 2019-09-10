function varargout = ClusteringMap(azimuthMap, elevationMap, threshold, type)
frqLen = size(azimuthMap,1);
timeLen = size(azimuthMap,2);

clusteringMap = NaN*zeros(frqLen,timeLen);
clusteringCntMap = zeros(frqLen,timeLen);
AngleList = [];
AngleCntList = [];
AngleTotList = [];
MapValList = [];
MapValIncrement = 1;
%% Analize each time-frequency bin
for covTimeIndx = 1:timeLen
    if strcmp(type, 'local') == 1
        AngleList = [];
        AngleCntList = [];
        AngleTotList = [];
        MapValList = [];
        MapValIncrement = 1;
    end
    for covFrqIndx = 1:frqLen
        Azimuth = azimuthMap(covFrqIndx,covTimeIndx);
        Elevation = elevationMap(covFrqIndx,covTimeIndx);
        if isnan(Azimuth) == false && isnan(Elevation) == false
            if isempty(AngleList) == true
                MapValList = MapValIncrement;
                AngleList = [Azimuth;Elevation];
                clusteringMap(covFrqIndx,covTimeIndx) = MapValIncrement;
                AngleCntList = 1;
                AngleTotList = AngleList;
            else
                [val,indx] = min(sqrt(min([abs(AngleList(1,:) - Azimuth);360-abs(AngleList(1,:) - Azimuth)]).^2 + abs(AngleList(2,:) - Elevation).^2));
                if val > threshold
                    MapValList = [MapValList,MapValList(end) + MapValIncrement];
                    AngleList = [AngleList, [Azimuth;Elevation]];
                    clusteringMap(covFrqIndx,covTimeIndx) = MapValList(end);
                    AngleCntList = [AngleCntList, 1];
                    AngleTotList = [AngleTotList, [Azimuth;Elevation]];
                else
                    clusteringMap(covFrqIndx,covTimeIndx) = MapValList(indx);
                    clusteringCntMap(covFrqIndx,covTimeIndx) = clusteringCntMap(covFrqIndx,covTimeIndx) + 1;
                    AngleCntList(indx) = AngleCntList(indx) + 1;
                    AngleTotList(:,indx) = AngleTotList(:,indx) + [Azimuth;Elevation];
                end
            end
        end
    end
end
%% OUTPUTS
varargout{1} = clusteringMap;