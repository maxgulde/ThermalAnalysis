% thermale Analyse ERNST 2.3-te Näherung
% Autor: Max Gulde
% Stand: 2017-09-20

% Annahmen
% - Energietransport zwischen internen Komponenten über Strahlung vernachlässigbar
%   - konservativer Ansatz, Strahlungstransport entspannt zusätzlich
% - Perfekte Wärmeleitung innerhalb eines Bauteils
%   - Zeitskala für Tempertaturänderungen langsam

% todo
% - ViewFactors einbetten für Aussenflächen in ThermoSim
% - Implement Earth view factor determination (self-shadow)

% Änderungen
% 2.4
% - Changed albedo to fitted function of Lockheed Martin's Data (includes
%   view factor)
% - Added view factors for internal radiation
% 2.3
% - kompatibel gemacht zun Sidos Berechnungen für Validierung via Thermal Desktop
% - Sol_Flux (hot case), fixer Albedowert, Celsius als Einheit, keine Zugriffszeiten, 
% - Plotunterschiede zwischen hot und cold
% - Fix: Albedokorrektur Subsolarwinkel
% - Fix: Maximaler Albedowinkel eingeführt
% - Fix: Zeitkonstante korrigiert für thermal coupling

%% Einstellungen
measureTime = 1;
if measureTime == 1
    tic
end

% % % Orbit
Inc = 98;   % [deg]
Alt = 700;  % [km]
RAA = 10;   % [deg]

% % % Satellit
satName = 'ERNST';              % Name des Satelliten
initTempC = 10;                 % Initial temperature in degrees C
T_Start = 273.15+initTempC;          % Starttemperatur des Satelliten
Sat_RadEffArea = 1.36;          % effektive Fläche des Radiators, Pyramide
Sat_RadName = 'Radiator';
Sat_CellEff = 0.34;             % Effizienz Solarzellen
Sat_CellName = 'SolarCells';
Sat_SurfNum = 1;                % Anzahl der Oberflächen wie in Simulation berechnet

% % % Darstellung
f_FigNum = 1;           % Nummer der figure
f_DrawOnTop = 0;        % drüberzeichnen
f_PlotGenPower = 0;
f_DrawCaseIndex = 2;    % 1: mean, 2: hot, 3: cold, 4: both
f_Verbose = 0;
f_UseMeanLoads = 0;
f_TemperatureLimits = [253 273 293 308];
f_TemperatureLimitsC = 'bgr';

% % % Optionen zum Zu- und Abschalten bestimmter Effekt
f_Sun = 1;      % Sonneneinstrahlung
f_Cmp = 0;      % Abwärme Komponenten
f_Alb = 1;      % Albedo
f_EIR = 1;      % Erde IR
f_Emi = 1;      % Emission Oberflächen
f_TCo = 1;      % Thermische Kopplung
f_XLo = 0;      % extra Loads, die nur kurzzeitig anfallen
f_IncludedParts = 0;    % Indices simulierte Strukturteile, 0 = alle
% f_DrawParts = [1:10 11];

% % % Simulation
t_Res = 60;                     % [s] zeitliche Auflösung
t_ResAcc = t_Res/6;             % [s] zeitliche Auflösung Zugriff, nur wenn f_UseMeanLoads == 0
t_Range = [0 1];                % simulierte Zeit [start ende], [0 1] voll
t_Step = 1;                     % Schrittweite
t_IntLim = [1/60 5];            % Grenzen der Zeitkonstanten für die Simulation der thermischen Kopplung

% % % Paths %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Testing simple cube sat
d_Par = 'Data\';
d_Cmp = [d_Par '_Komponenten.txt'];
d_XLo = [d_Par '_XLoads.txt'];
%d_TEx = [d_Par 'Temperatur'];
% d_Suff = ' - short';
d_Suff = '';

% % % Paths
d_Bas = sprintf('Data_%s_i%i_a%i_r%i_t%i',satName,Inc,Alt,RAA,t_Res); % Base path
% EXPERIMENT
d_Bas = 'Functional_Cube';
d_StrFolder = fullfile(d_Bas,'Structure');

% % % Data Files
d_SubSol = fullfile(d_Bas,'SunAngles.csv');
d_EarthAngles = fullfile(d_Bas,'EarthAngles.csv');
d_AreaS = fullfile(d_Bas,'Out_AreaSunView.txt');
d_Mat = fullfile(d_Bas,'_materials.txt');        % Materials file
d_Str = fullfile(d_StrFolder,'_structure.txt');      % Structure file
d_TCo = fullfile(d_Bas,'_conduction.txt');    % Thermal conductivity
d_IntRad = fullfile(d_Bas,'matrix.txt'); % Internal viewfactors

% Daten
fmt_Base = sprintf('_i%.0f_a%.0f_r%02.0f_t%03.0f',Inc,Alt,RAA,t_Res);
%d_TEx = [d_TEx fmt_Base '.txt'];
d_Power = sprintf('Out_Power_i98_1929_a700_r10_t120.txt');
d_Access = sprintf('SimpleSat AccessTimes.csv');

% Formatstrings
fstrAreas = ['%d %s %d %12s' repmat(';%f',1,Sat_SurfNum)];
fstrPower = '%d %s %d %12s;%f';
fstrAccess = '%f,%f %s %f %f:%f:%f,%f %s %f %f:%f:%f,%f';
fstrDateOut = '%.0f %s %.0f %02.0f:%02.0f:%02.0f';
fstrCmp = '%s %d %f %f';
fstrXLo = '%s %f %f %f';

% % % Format Strings
fstrMat = '%s %f %f %f'; % Name Abs Emi HCap
fstrStr = '%s %s %s %f %f %s %d'; % Name Optical Bulk TotArea Weight Component Internal
fstrStrSurf = '%s %f %f %f %f %d'; % Optical Area Nx Ny Nz in_flag
fstrSubSol = '%d %s %d %12s,%f,%f,%f'; % D(d) M(s) Y(d) t(12s), Az(f), El(f), SubSol(f)
fStrEarthAngles = '%d %s %d %12s,%f,%f'; % D(d) M(s) Y(d) t(12s), Az(f), El(f)
fstrTCo = '%s %s %f %f %f';
fstrDateIn = '%d %s %d %s';
fstrDate = 'dd mmm yyyy HH:MM:SS';
fstrOrder = '%s %s %s %s';

% Header Zeilen
h_Area = Sat_SurfNum + 2;
h_Power = 3;

% % % Eingabe
f_ReloadAllData = 1;            % alle Daten (Winkel, Leistung, etc) einladen (langsam)
f_ReloadMatData = 1;            % Materialdaten einladen (schnell)
f_IgnoreAccessIntervals = 0;    % ignore increased energy load by communication devices
f_UseFixedAlbedo = -1;          % set to < 0 to use correction table
f_UseCelsius = 273.15;           % set to == 0 to use Kelvin


% % % Ausgabe
f_GenerateTFile = 0;
f_ExportPlot = 0;
f_SaveResults = 0;
f_PlotLineWidth = 2;

% Laden von Daten
if (f_ReloadAllData == 1 || f_ReloadMatData)
    Sat_Mat = struct();         % Materialien, Absorption und Emision
    Sat_Cmp = struct();         % Komponenten, Wärmeentwicklung
    Sat_Struct = struct();      % Oberflächen, Verknüpfung mit Mat und Cmp, Flächen und Gewicht
end

% % % Data from "SC Thermal Control Handbook", p. 28 (24h, (C)old and (H)ot
% % % cases)
% Direct sunlight
Sol_Flux = [1322 1414];         % [W/m^2] Cold and hot cases
T_Space = 2.7;                  % [K] Temperature of space
% Albedo
Alb_CorrSubsolar = [0:10:90; 0:0.01:0.05 0.08 0.13 0.2 0.31]';
Alb_OrbitInc_C = [30 60 90; 0.17 0.2 0.2]';
Alb_OrbitInc_H = [30 60 90; 0.19 0.23 0.23]';
% Earth IR
EIR_OrbitInc_C = [30 60 90; 236 226 225]';  % [W/m^2]
EIR_OrbitInc_H = [30 60 90; 257 241 230]';  % [W/m^2]

% Physical Constants
consts;
re = 6371000; % [m] Earth's Radius
r = re/(re+Alt*1e3);

% Zeit
t_PlotScale = 24*3600;
t_DaySec = 24*3600;

% maximaler Subsolarwinkel für Albedo
AlbMaxAngle = 100; % [deg] (Above 100 deg the formula generates complex numbers)

%% Daten einlesen
if (f_ReloadAllData == 1)
    fprintf('Einlesen der Daten:\n');
    
    % % Subsolar Angles
    fprintf(' ... Subsolar Angles ...\n');
    dat_temp = ReadCSV(d_SubSol,fstrSubSol,1);
    dat_SolAng = dat_temp(end-2:end-1); % Azimuth, Elev % Not required
    dat_SubSol = dat_temp(end);
    
    % Check temporal resolution of Subsolar Angles
    subSolStartDate = sprintf(fstrDateIn,dat_temp{1,1}(1),...
        dat_temp{1,2}{1},dat_temp{1,3}(1),dat_temp{1,4}{1});
    subSolStepDate = sprintf(fstrDateIn,dat_temp{1}(2),dat_temp{2}{2},...
        dat_temp{1,3}(2),dat_temp{1,4}{2});
    sEnd = length(dat_temp{1});
    subSolEndDate = sprintf(fstrDateIn,dat_temp{1}(sEnd),...
        dat_temp{2}{sEnd},dat_temp{1,3}(sEnd),dat_temp{1,4}{sEnd});
    subSolStartTime = datevec(subSolStartDate,fstrDate);
    subSolStepTime = datevec(subSolStepDate,fstrDate);
    subSolEndTime = datevec(subSolEndDate,fstrDate);
    subSolDuration = duration(subSolStepTime(4:6)-subSolStartTime(4:6));
    subSolRes = seconds(subSolDuration);
    
    % Get full time data
    dat_Time{4} = dat_temp{1,4};
    dat_Time{3} = dat_temp{1,3};
    dat_Time{2} = dat_temp{1,2};
    dat_Time{1} = dat_temp{1,1};
    
    % Compare against input
    if subSolRes ~= t_Res
        fprintf(['\n !!! WARNING! !!!\n' ' !!! t_Res does not match the '...
            'temporal resolution of the Solar Angles file. '...
            'Using temporal resolution from file. !!!\n\n']);
        t_Res = subSolRes;
    end
    
    % % Local Zenith
    fprintf(' ... Local Zenith ...\n');
    dat_temp = ReadCSV(d_EarthAngles,fStrEarthAngles,1);
    dat_EarthAngles = dat_temp(end-1:end);
    
    % Check temporal resolution of Local Zenith
    zenithStartDate = sprintf(fstrDateIn,dat_temp{1,1}(1),...
        dat_temp{1,2}{1},dat_temp{1,3}(1),dat_temp{1,4}{1});
    zenithStepDate = sprintf(fstrDateIn,dat_temp{1}(2),dat_temp{2}{2},...
        dat_temp{1,3}(2),dat_temp{1,4}{2});
    zEnd = length(dat_temp{1});
    zenithEndDate = sprintf(fstrDateIn,dat_temp{1}(zEnd),...
        dat_temp{2}{zEnd},dat_temp{1,3}(zEnd),dat_temp{1,4}{zEnd});
    zenithStartTime = datevec(zenithStartDate,fstrDate);
    zenithStepTime = datevec(zenithStepDate,fstrDate);
    zenithEndTime = datevec(zenithEndDate,fstrDate);
    zenithDuration = duration(zenithStepTime(4:6)-zenithStartTime(4:6));
    zenithRes = seconds(zenithDuration);
    
    % Compare against Subsolar
    if subSolRes ~= zenithRes
        error(['The temporal resolution of the Earth Angles file does '...
            'not match the Solar Angles file.\n']);
    end
    if any((zenithEndTime - subSolEndTime) ~= 0) || ...
            any((zenithStartTime - subSolStartTime) ~= 0)
        error(['The Solar Angles dates do not match the ones from the '...
            'Earth Angles file.\n']);
    end
    
    % % % % % % % % % Moved to the next loading step % % % % % % % % % % %
    % Flächen und Leistung
    % Reihenfolge extrahieren
    %     fid = fopen(d_AreaS);
    %     dat_temp = textscan(fid, fstrOrder, h_Area-2, 'HeaderLines', 2);
    %     dat_Order = dat_temp(3);
    %     fclose(fid);
    %
    %     % Flächen extrahieren
    %     fprintf(' ... Flächen aus Sonnensicht ...\n');
    %     dat_temp = ReadCSV(d_AreaS,fstrAreas,h_Area);
    %     dat_AreaS = dat_temp(5:end);
    %     dat_Time = dat_temp(1:4);
    %     %fprintf(' ... Flächen aus Erdsicht ...\n');
    %     %dat_temp = ReadCSV(d_AreaE,fstrAreas,h_Area);
    %     %dat_AreaE = dat_temp(5:end);
    %     if (f_PlotGenPower == 1)
    %        fprintf(' ... Erzeugte Leistung ...\n');
    %        dat_temp = ReadCSV(d_Power,fstrPower,h_Power);
    %        dat_Power = dat_temp(5:end);
    %     end
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    
    % % Zugriffszeiten
    if (f_IgnoreAccessIntervals == 0)
        fprintf(' ... Zugriffszeiten ...\n');
        dat_temp = ReadCSV(d_Access,fstrAccess,1);
        % Übersetzen in Sekunden seit Start
        dat_Access = zeros(numel(dat_temp{1}),2); % Start [s] und Ende [s]
        startDate = sprintf(fstrDateIn,dat_Time{1}(1),dat_Time{2}{1},dat_Time{3}(1),dat_Time{4}{1}(1:end-4));
        endDate = sprintf(fstrDateIn,dat_Time{1}(end),dat_Time{2}{size(dat_Time{1},1)},dat_Time{3}(end),dat_Time{4}{size(dat_Time{1},1)}(1:end-4));
        startTime = datevec(startDate,fstrDate);
        for ii = 1:numel(dat_temp{1})
            t = datevec(sprintf(fstrDateOut,dat_temp{2}(ii),dat_temp{3}{ii},dat_temp{4}(ii),dat_temp{5}(ii),dat_temp{6}(ii),dat_temp{7}(ii)),fstrDate);
            dat_Access(ii,1) = etime(t,startTime);
            dat_Access(ii,2) = dat_Access(ii,1) + round(dat_temp{14}(ii));
        end
        ACtot = sum(dat_temp{:,14});
    end

end

if (f_ReloadMatData == 1 || f_ReloadAllData == 1)
    Sat_Mat = struct();     % Materials, Absorption and Emission
    Sat_Struct = struct();
    
    % % Structure: Name, materials, size, etc.
    fprintf(' ... Reading Structure ...\n');
    fid = fopen(d_Str);
    dat_temp = textscan(fid,fstrStr,'CommentStyle','%');
    fclose(fid);
    StructNum = numel(dat_temp{1});
    % Number of simulated parts
    if (f_IncludedParts == 0)
        f_IncludedParts = 1:StructNum;
    end
    
    fstrAreas = ['%d %s %d %12s' repmat(';%f',1,StructNum)];
    h_Area = StructNum + 2;
    
    % SunView
    fprintf(' ... Reading Structure in Sunlight ...\n');
    dat_temp2 = ReadCSV(d_AreaS,fstrAreas,h_Area);
    dat_AreaS = dat_temp2(5:end);
    dat_Time = dat_temp2(1:4);
    
    % Extract order
    fid = fopen(d_AreaS);
    dat_temp2 = textscan(fid, fstrOrder, h_Area-2, 'HeaderLines', 2);
    dat_Order = dat_temp2(3);
    fclose(fid);
    
    fprintf(' ... Creating Structure ...\n');
    for ii = 1:StructNum
        Sat_Struct(ii).name = dat_temp{1}{ii};
        Sat_Struct(ii).optical = dat_temp{2}{ii};	% Cmps have NA surface
        Sat_Struct(ii).bulk = dat_temp{3}{ii};
        Sat_Struct(ii).size = dat_temp{4}(ii);
        Sat_Struct(ii).mass = dat_temp{5}(ii);
        Sat_Struct(ii).cmp = strsplit(dat_temp{6}{ii},',');                % CHECK
        Sat_Struct(ii).internal = dat_temp{7}(ii);
        % Position in area files
        Sat_Struct(ii).AFileIdx = find(strcmp(dat_Order{:},Sat_Struct(ii).name));
        if (isempty(Sat_Struct(ii).AFileIdx))
            fprintf('\t\t+Interne Struktur: %s \n',Sat_Struct(ii).name);
        end
        
        % Read individual surfaces
        if strcmp(Sat_Struct(ii).cmp,'NA')
            temp_path = fullfile(d_StrFolder,[Sat_Struct(ii).name '.txt']);
            if exist(temp_path,'file') ~= 0
                fid = fopen(temp_path);
                dat_temp2 = textscan(fid,fstrStrSurf,'CommentStyle','%');
                fclose(fid);
                tLength = length(dat_temp2{1});
                tempStruct = struct();
                
                for iii = 1:tLength
                    tempStruct(iii).optical = dat_temp2{1}{iii};
                    tempStruct(iii).area = dat_temp2{2}(iii);
                    tempStruct(iii).normalV = ...
                        [dat_temp2{3}(iii), dat_temp2{4}(iii), dat_temp2{5}(iii)];
                    tempStruct(iii).internal = dat_temp2{6}(iii);
                end
                Sat_Struct(ii).surfaces = tempStruct;
            end
        else
            Sat_Struct(ii).surfaces = struct([]);
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % % Komponenten
    fprintf(' ... Komponenten ...\n');
    fid = fopen(d_Cmp);
    dat = textscan(fid, fstrCmp, 'CommentStyle','%');
    fclose(fid);
    for ii = 1:numel(dat{1})
        Sat_Cmp(ii).name = dat{1}{ii};
        Sat_Cmp(ii).mode = dat{2}(ii);
        if (Sat_Cmp(ii).mode == 1)
            Sat_Cmp(ii).loadFrac = ACtot / etime(datevec(endDate),datevec(startDate));
        else
            Sat_Cmp(ii).loadFrac = 0;
        end
        Sat_Cmp(ii).loadActive = dat{3}(ii);
        Sat_Cmp(ii).loadIdle = dat{4}(ii);
    end
    
    % % Materialien
    fprintf(' ... Materialdaten ...\n');
    fid = fopen(d_Mat);
    dat = textscan(fid, fstrMat, 'CommentStyle','%');
    fclose(fid);
    for ii = 1:numel(dat{1})
        Sat_Mat(ii).name = dat{1}{ii};
        Sat_Mat(ii).abs = dat{2}(ii);
        Sat_Mat(ii).emi = dat{3}(ii);
        Sat_Mat(ii).cap = dat{4}(ii);
    end
    
    % % Thermische Kopplung der Oberflächen
    if (f_TCo == 1)
        fprintf(' ... Thermische Kopplung ...\n');
        fid = fopen(d_TCo);
        dat = textscan(fid, fstrTCo, 'CommentStyle','%');
        fclose(fid);
        % 1 Kopplungsmatrix erstellen
        Sat_TCo = cell(numel(Sat_Struct));
        % 2 Standard Struktur reinschreiben
        s = struct();
        s.Mode = 0;
        for ii = 1:numel(Sat_Struct)
            for jj = 1:numel(Sat_Struct)
                Sat_TCo{ii,jj} = s;
            end
        end
        % 3 Mit Werten aus Liste füllen
        for ii = 1:numel(dat{1})
            % Indices finden für beide Oberflächen
            idx1 = find(strcmp({Sat_Struct.name}, dat{1}{ii}));
            idx2 = find(strcmp({Sat_Struct.name}, dat{2}{ii}));
            % Daten eintragen
            s = struct();
            s.Mode = 1;
            s.Length = dat{3}(ii);
            s.Area = dat{4}(ii);
            s.Cond = dat{5}(ii);
            % Speichern
            Sat_TCo{idx1,idx2} = s;
            Sat_TCo{idx2,idx1} = s;
        end
    end
    
    % % Internal Radiation
    fprintf(' ... Internal Radiation ...\n');
    fid = fopen(d_IntRad);
    dat_temp = textscan(fid,'%f %f %f %f %f %f %f %f %f','CommentStyle','%');
    fclose(fid);
    internal_vf = zeros(StructNum);
    for iii = 1:StructNum
        internal_vf(:,iii) = dat_temp{1,iii};
    end
    % Normalization for fully internal components
    for iii = 1:StructNum
        if Sat_Struct(iii).internal == 1
            internal_vf(iii,:) = internal_vf(iii,:)./(sum(internal_vf(iii,:)));
        end
    end
    
    % % extra Lastfälle
    fprintf(' ... Lastfälle ...\n');
    fid = fopen(d_XLo);
    dat_XLo = textscan(fid, fstrXLo, 'CommentStyle','%');
    fclose(fid);
    
    % Full weight
    fprintf(' Full weight:  %.2f kg ...\n',sum([Sat_Struct(:).mass]));
    % from components
    mass = 0;
    for ii = 1:numel(Sat_Struct)
        if (strcmp(Sat_Struct(ii).optical,'NA') == 1)
            mass = mass + Sat_Struct(ii).mass;
        end
    end
    fprintf(' \t Structure:  %.2f kg ...\n',sum([Sat_Struct(:).mass]) - mass);
    fprintf(' \t Components: %.2f kg ...\n',mass);
    
    % Print complete structure
    fprintf('\t Structure:\n');
    for ii = 1:StructNum
        fprintf('\t\t%d: %s\n',ii,Sat_Struct(ii).name);
    end
    
end 

%% Zeitvektor

t_Vec = (1:size(dat_Time{1},1)) * t_Res / t_PlotScale;
N = numel(t_Vec);
% Simulationszeitraum festlegen
t_sim = floor(t_Range*N);
if (t_sim(1) < 1)
    t_sim(1) = 1;
end
if (t_sim(end) > N)
    t_sim(end) = N;
end
ran = t_sim(1):t_Step:t_sim(end);
t_Vec = t_Vec(ran);
% Startzeit
t_Offset = t_Vec(1) * t_PlotScale;

%% Simulation    

fprintf('Simulation ');

% Temperaturen für hot (1) und cold (2) case, jede Fläche, jeder Zeitschritt
T = zeros(2,numel(Sat_Struct),N+1);
% Starttemperatur
T(:,:,ran(1)) = T_Start;
PowerValues = T;

% Oberflächen sortieren
f_IncludedParts = sort(f_IncludedParts);

% Delta t
dT = t_Res * t_Step;

% Correction for the inclination
inc_c = Inc;
if(Inc > 90)
    inc_c = 180 - Inc;
end

% Influence from the inclination (Earth IR)
facInc_C_EIR = EIR_OrbitInc_C(find(EIR_OrbitInc_C(:,1) <= inc_c,1,'last'),2);
facInc_H_EIR = EIR_OrbitInc_H(find(EIR_OrbitInc_H(:,1) <= inc_c,1,'last'),2);

% Influence from the inclination (Albedo)
facInc_C_Alb = Alb_OrbitInc_C(find(Alb_OrbitInc_C(:,1) <= inc_c,1,'last'),2);
facInc_H_Alb = Alb_OrbitInc_H(find(Alb_OrbitInc_H(:,1) <= inc_c,1,'last'),2);


% Time-steps
modDivisor = floor(numel(ran)/10);
for tt = ran
    % Print progess
    if mod(tt,modDivisor) == 0
        fprintf('.');
    end
    
    % Current time (used for internal components)
    %t_Now = (tt - t_Step) * dT;
    
    % Local zenith
    [earthPos(1),earthPos(2),earthPos(3)] = ...
        sph2cart(deg2rad(dat_EarthAngles{1}(tt)), ...
        deg2rad(dat_EarthAngles{2}(tt)),1);
    % Negative of the vector that points towards the Earth
    earthPos = - earthPos;
    
    % Add correction for albedo (from subsolar angle)
    if (f_UseFixedAlbedo < 0)
        a_c = Alb_CorrSubsolar(find(Alb_CorrSubsolar(:,1) <= dat_SubSol{1}(tt),1,'last'),2);
        facInc_C_Alb_tt = facInc_C_Alb + a_c;
        facInc_H_Alb_tt = facInc_H_Alb + a_c;
    else
        facInc_C_Alb_tt = f_UseFixedAlbedo;
        facInc_H_Alb_tt = f_UseFixedAlbedo;
    end
    
    % % % (A) Isolierte Betrachtung
    for ss = f_IncludedParts
        % Received power
        P_C = 0;
        P_H = 0;
        
        % Get materials for ss
        sIdx = find(strcmp({Sat_Mat(:).name}',Sat_Struct(ss).optical));
        bIdx = find(strcmp({Sat_Mat(:).name}',Sat_Struct(ss).bulk));
        if (isempty(sIdx))
            fprintf('Error: Material <%s> not found. Changed to <%s>\n.',Sat_Struct(ss).optical, Sat_Mat(1).name);
            sIdx = 1;
        end
        if (isempty(bIdx))
            fprintf('Error: Material <%s> not found. Changed to <%s>\n.',Sat_Struct(ss).bulk, Sat_Mat(1).name);
            bIdx = 1;
        end
        
        % % % (1) Direct Sunlight (Optimized by Out_AreaSunView)
        if f_Sun == 1
            % Index der Komponente in Flächendatei
            cIdx = Sat_Struct(ss).AFileIdx;
            % Wenn interne Struktur, dann auf -1 setzen
            if (isempty(cIdx))
                cIdx = -1;
            end
            
            if (cIdx > 0)
                A = dat_AreaS{1,cIdx}(tt);
            else
                A = 0;
            end
            
            %             % geringerer Energieeingang für Solarzellen
            %             if (strcmp(Sat_Struct(ss).name(1:end-1),Sat_CellName) == 1)
            %                 A = A * (1 - Sat_CellEff);
            %             end

            ss_abs = Sat_Mat(sIdx).abs;
            ss_emi = Sat_Mat(sIdx).emi;
            
            P_C = P_C + A * ss_abs * Sol_Flux(1);
            P_H = P_H + A * ss_abs * Sol_Flux(2);
        end
        
        % % % (2) Komponenten
        if ss == 2 %%%%%%%%%% EXPERIMENTAL
            P_C = P_C + 2;
            P_H = P_H + 2;
        end
        if (f_Cmp == 1)
            % Durchschnittsverbrauch ermitteln
            P = 0;
            % Alle zugeordneten Komponenten durchgehen
            for cc = Sat_Struct(ss).cmp
                % Nur, wenn Komponenten angegeben sind
                if (strcmp(cc,'NA') == 0)
                    % Index finden
                    comIdx = find(strcmp({Sat_Cmp(:).name}',cc));
                    if (isempty(comIdx))
                        fprintf('Error: Komponente <%s> nicht gefunden.',cc);
                        comIdx = 1;
                    end
                    % Last addieren
                    if (Sat_Cmp(comIdx).mode == 0) % konstanter Fall
                        P = P + Sat_Cmp(comIdx).loadActive;
                    else % nicht konstanter Fall
                        if (f_UseMeanLoads == 1) % gemittelt
                            P = P + (Sat_Cmp(comIdx).loadActive * Sat_Cmp(comIdx).loadFrac + Sat_Cmp(comIdx).loadIdle * (1 - Sat_Cmp(comIdx).loadFrac)); 
                        else
                            switch (Sat_Cmp(comIdx).mode)
                                case 1 % zeitaufgelöst, nach Access
                                    dur = 0;
                                    % Sekunden zählen, die die Komponente eingeschaltet ist
                                    for tt_c = tt*t_Res:t_ResAcc:(tt+1)*t_Res
                                        idx = find(dat_Access(:,1) <= tt_c & dat_Access(:,2) > tt_c);
                                        % add
                                        if (~isempty(idx))
                                            dur = dur + t_ResAcc;
                                        end
                                    end
                                    % normieren auf Zeitskala
                                    dur = dur / (t_Res);
                                    if (dur > 1)
                                        dur = 1;
                                    end
                                    % in Leistung umrechnen
                                    P = P + Sat_Cmp(comIdx).loadActive * dur + Sat_Cmp(comIdx).loadIdle * (1 - dur);
%                                 case 2 % zeitaufgelöst, nach Sentinelwinkel
%                                     if (dat_Sentinel{1}(tt) < ACrit)
%                                         P = P + Sat_Cmp(comIdx).loadActive;
%                                     else
%                                         P = P + Sat_Cmp(comIdx).loadIdle;
%                                     end
                            end
                        end
                    end
                    % % % Extra Last
                    if (f_XLo == 1)
                        xIdx = find(strcmp(dat_XLo{1},cc));
                        % Check ob vorhanden
                        if (~isempty(xIdx))
                            for ii = xIdx'
                                % Check ob in Intervall
                                if (dat_XLo{2}(ii) <= t_Now && t_Now - dat_XLo{2}(ii) < dat_XLo{3}(ii))
                                    P = P + dat_XLo{4}(ii);
                                end
                            end
                        end
                    end
                end
            end
            % Leistungszuwachs
            P_C = P_C + P;
            P_H = P_H + P;
        end
        
        % % % (3) IR Emission (internal only, optimized by _internal_matrix)
        if f_Emi == 1
            % IGNORED CODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % INCORPORATE PROPER IMPLEMENTATION
            % Update area to effective area in case of modified radiator
            %if (strcmp(Sat_Struct(ss).name,Sat_RadName) == 1)
            %    A = A * Sat_RadEffArea;
            %end
            % IGNORED CODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            for pp = f_IncludedParts
                if ss == pp % Object is itself, calculation would be 0, skip
                    continue;
                end
                % Object is another volume, use view factor
                pIdx = find(strcmp({Sat_Mat(:).name}',Sat_Struct(pp).optical));
                pp_emi = Sat_Mat(pIdx).emi;
                pp_A = Sat_Struct(pp).size; % This is the total surf. area
                
                % Internal view factor (source, target)
                % Radiated power
                P_C = P_C + ksb * (T(1,pp,tt)^4 - T(1,ss,tt)^4) * ss_emi * pp_emi * pp_A * internal_vf(pp,ss);
                P_H = P_H + ksb * (T(2,pp,tt)^4 - T(2,ss,tt)^4) * ss_emi * pp_emi * pp_A * internal_vf(pp,ss);
            end
        end
        
        % % % IR Emissions, Earth IR & Albedo
        % Check if the piece is divided into several surfaces
        if numel(Sat_Struct(ss).surfaces) > 0
            for int_ss = 1:numel(Sat_Struct(ss).surfaces) % Cycle through all the surfaces
                % Only take into account external surfaces
                if Sat_Struct(ss).surfaces(int_ss).internal == 0
                    % Area
                    A = Sat_Struct(ss).surfaces(int_ss).area;
                    
                    % Optical Coating
                    sIdx = find(strcmp({Sat_Mat(:).name}',Sat_Struct(ss).surfaces(int_ss).optical));
                    int_ss_emi = Sat_Mat(sIdx).emi;
                    int_ss_abs = Sat_Mat(sIdx).abs;
                    
                    % Normal Vector
                    normalV = Sat_Struct(ss).surfaces(int_ss).normalV;
                    
                    % View Factor
                    rho = rad2deg(atan2(norm(cross(normalV,earthPos)), ...
                        dot(normalV,earthPos)));
                    vF = viewFactor(r,rho);
                    
                    % % % (3) IR Emission (external)
                    if f_Emi == 1
                        P_C = P_C - ksb * (T(1,ss,tt)^4 - T_Space^4) * A * int_ss_emi;
                        P_H = P_H - ksb * (T(2,ss,tt)^4 - T_Space^4) * A * int_ss_emi;
                    end
                    
                    % % % (4) Earth IR
                    if f_EIR == 1
                        P_C = P_C + A * int_ss_emi * facInc_C_EIR .* vF;
                        P_H = P_H + A * int_ss_emi * facInc_H_EIR .* vF;
                    end
                    
                    % % % (5) Albedo
                    if f_Alb == 1 && dat_SubSol{1}(tt) < AlbMaxAngle
                        % Compute Albedo flux
                        AlbedoFlux_C = albedoSolarFlux(Sol_Flux(1).*facInc_C_Alb_tt,dat_SubSol{1}(tt),vF);
                        AlbedoFlux_H = albedoSolarFlux(Sol_Flux(2).*facInc_H_Alb_tt,dat_SubSol{1}(tt),vF);
                        
                        P_C = P_C + A * int_ss_abs * AlbedoFlux_C;
                        P_H = P_H + A * int_ss_abs * AlbedoFlux_H;
                    end
                end
            end
        else % This part is not divided into individual surfaces
            % Would there be data where there it's not divided into
            % surfaces? Add non separated case (?)
        end
        
        % % % (X) Temperature Calculation
        
        % Energy change calculation
        dE_C = P_C * dT;
        dE_H = P_H * dT;
        
        % Temperature change calculation
        dT_C = dE_C / (Sat_Mat(bIdx).cap * Sat_Struct(ss).mass);
        dT_H = dE_H / (Sat_Mat(bIdx).cap * Sat_Struct(ss).mass);
        T(1,ss,tt+t_Step) = T(1,ss,tt) + dT_C;
        T(2,ss,tt+t_Step) = T(2,ss,tt) + dT_H;
        PowerValues(1,ss,tt) = P_C;
        PowerValues(2,ss,tt) = P_H;
        
        % Prevent negative temperatures
        for ii = 1:2
            if (T(ii,ss,tt+t_Step) < T_Space)
                T(ii,ss,tt+t_Step) = T_Space;
            end
        end 
    end
    
    % % % (B) Wärmeleitung, etc. zwischen Einzelkomponenten
    if (f_TCo == 1)
        % alle betrachteten Strukturen
        for ss = f_IncludedParts
            % finde Index, nur ab dort
            idx = find(f_IncludedParts == ss);
            % nur solche Strukturen, die auch simuliert werden
            for pp = f_IncludedParts(idx:end)
                % nicht auf sich selbst anwenden
                if (ss == pp)
                    continue;
                end
                % checken, ob Kopplung vorhanden
                if (Sat_TCo{ss,pp}.Mode == 0)
                    continue;
                end
                % Bulk Material Parameter Index
                bIdx_s = find(strcmp({Sat_Mat(:).name}',Sat_Struct(ss).bulk));
                bIdx_p = find(strcmp({Sat_Mat(:).name}',Sat_Struct(pp).bulk));
                % Parameter extrahieren
                m_s = Sat_Struct(ss).mass;
                m_p = Sat_Struct(pp).mass;
                hc_s = Sat_Mat(bIdx_s).cap;
                hc_p = Sat_Mat(bIdx_p).cap;
                C_s = m_s * hc_s;
                C_p = m_p * hc_p;
                dt = t_Res * t_Step;
                % Zeitskala abschätzen, auf der Temperaturausgleich stattfindet
                R = Sat_TCo{ss,pp}.Length / (Sat_TCo{ss,pp}.Area * Sat_TCo{ss,pp}.Cond);
                tau = R * mean([C_s C_p]);
                t_Ratio = dt / tau;
                % heißer und kalter Fall
                for cc = 1:2
                    % schneller Ausgleich -> beide Strukturen haben dieselbe Temperatur
                    if (t_Ratio > t_IntLim(2))
                        T_prime = (C_s * T(cc,ss,tt+t_Step) + C_p * T(cc,pp,tt+t_Step))/(C_s + C_p);
                        T(cc,ss,tt+t_Step) = T_prime;
                        T(cc,pp,tt+t_Step) = T_prime;
                    elseif (t_Ratio < t_IntLim(1)) % langsamer Ausgleich -> lineare Näherung
                        deltaT = T(cc,pp,tt+t_Step) - T(cc,ss,tt+t_Step);
                        T(cc,ss,tt+t_Step) = T(cc,ss,tt+t_Step) + dt/C_s*deltaT/R;
                        T(cc,pp,tt+t_Step) = T(cc,pp,tt+t_Step) - dt/C_p*deltaT/R;
                    else % ungefähr gleiche Zeitskala, runterbrechen
                        dt2 = t_IntLim(1) * dt;
                        for ii = dt2:dt2:dt
                            deltaT = T(cc,pp,tt+t_Step) - T(cc,ss,tt+t_Step);
                            T(cc,ss,tt+t_Step) = T(cc,ss,tt+t_Step) + dt2/C_s*deltaT/R;
                            T(cc,pp,tt+t_Step) = T(cc,pp,tt+t_Step) - dt2/C_p*deltaT/R;
                        end                        
                    end
                    % Verhältnis von dt zu tau angeben
                    if (f_Verbose == 1 && tt == ran(1) && cc == 1)
                        fprintf('\nWärmeleitung: <%s>, <%s> (dt/tau = %f).',Sat_Struct(ss).name, Sat_Struct(pp).name, t_Res/tau);
                    end
                end
            end
            % Negative Temperaturen verhindern
            for ii = 1:2
                if (T(ii,ss,tt+t_Step) < T_Space)
                    T(ii,ss,tt+t_Step) = T_Space;
                end
            end 
        end
    end
    
end
fprintf('\n ... fertich.\n');

    
%% Plotten der Ergebnisse
    
fprintf('Ergebnisse zeichnen ...');
f_DrawParts = f_IncludedParts;

figure(f_FigNum);
if (f_DrawOnTop == 1)
    hold on;
end

tempOffset = 0;
if (f_UseCelsius ~= 0)
    tempOffset = f_UseCelsius;
end    

if (f_DrawCaseIndex == 4)
    leg = cell(numel(f_DrawParts) * 2,1);
else
    leg = cell(numel(f_DrawParts),1);
end
c = 1;
for ii = f_IncludedParts(f_DrawParts)
    % cold case
    switch (f_DrawCaseIndex)
        case 1 % mean case
            y = squeeze(mean(T(:,ii,ran(1) + 1:end),1));
            y = y(y ~= 0);
            plot(t_Vec, y - tempOffset,'LineWidth',f_PlotLineWidth);
            if (ii == f_IncludedParts(f_DrawParts(1)))
                hold on;
            end
            % Legende
            leg{c} = [Sat_Struct(ii).name ' (mean)'];
            c = c + 1;
        case 3 % cold case
            y = squeeze(T(:,ii,ran(1) + 1:end));
            plot(t_Vec, y(1,(y(1,:) ~= 0)) - tempOffset,'LineWidth',f_PlotLineWidth);
            if (ii == f_IncludedParts(f_DrawParts(1)))
                hold on;
            end
            % Legende
            leg{c} = [Sat_Struct(ii).name ' (cold)'];
            c = c + 1;
        case 2 % hot case
            y = squeeze(T(:,ii,ran(1) + 1:end));
            plot(t_Vec, y(2,(y(2,:) ~= 0)) - tempOffset,'LineWidth',f_PlotLineWidth);
            if (ii == f_IncludedParts(f_DrawParts(1)))
                hold on;
            end
            % Legende
            leg{c} = [Sat_Struct(ii).name ' (hot)'];
            c = c + 1; 
        case 4 % both cases
            y = squeeze(T(:,ii,ran(1) + 1:end));
            plot(t_Vec, y(1,(y(1,:) ~= 0)) - tempOffset,'LineWidth',f_PlotLineWidth);
            if (ii == f_IncludedParts(f_DrawParts(1)))
                hold on;
            end
            plot(t_Vec, y(2,(y(2,:) ~= 0)) - tempOffset,'LineWidth',f_PlotLineWidth);
            % Legende
            leg{c} = [Sat_Struct(ii).name ' (cold)'];
            c = c + 1;
            leg{c} = [Sat_Struct(ii).name ' (hot)'];
            c = c + 1;
    end
end
% xlim([t_Vec(1) t_Vec(end)]);
% zeichne Temperaturegrenzen
for tt = 2:numel(f_TemperatureLimits)
    x = xlim;
    y = f_TemperatureLimits;
    p = patch([x(1) x(2) x(2) x(1)],[y(tt-1) y(tt-1) y(tt) y(tt)] - tempOffset,f_TemperatureLimitsC(tt-1));
    set(p,'FaceAlpha',0.05,'EdgeAlpha',0);
    %plot(xlim,[f_TemperatureLimits(tt) f_TemperatureLimits(tt)],[f_TemperatureLimitsC(tt) '--']);
end
hold off;
% Titel und Achsen
title(sprintf('inc = %.0f deg, alt = %.0f km, raan = %.0f deg, temp res = %.0f s',Inc,Alt,RAA,t_Res));
xlabel 'time [d]'
if (f_UseCelsius == 0)
    ylabel 'temperature [K]'
else
    ylabel 'temperature [°C]'
end
legend(leg);
set(gcf,'Color','w');
set(gca,'FontName','Calibri','FontSize',18);

% exportiere Plot
if (f_ExportPlot == 1)
    fprintf('Exportiere ...\n');
    figure(f_FigNum);
    name = sprintf('thermal_a%.0fkm_i%.0fdeg_r%03.0fdeg_t%.0fs.png',Alt,Inc,RAA,t_Res * t_Step);
    fprintf(' ... %s\n',name);
    export_fig(name,'-r150');
end

% Plot der generierten Leistung
if (f_PlotGenPower == 1)
    figure(f_FigNum + 1);
    if (f_DrawOnTop == 1)
        hold on;
    end
    % plot solar panel power generation
    relPower = dat_Power{1}';
    relPower = relPower(ran);
%     relPower = movingmean(relPower,t_DaySec/(t_Res * t_Step),2,[]);
%     relPower = smooth(relPower,10);
    plot(t_Vec, relPower);
    xlim([t_Vec(1) t_Vec(end)]);
    xlabel 'time [d]'
    ylabel 'daily mean solar power [W]'
    hold off;
    
    % exportiere Plot
    if (f_ExportPlot == 1)
        fprintf('Exportiere ...\n');
        figure(f_FigNum + 1);
        name = sprintf('power_a%.0fkm_i%.0fdeg_r%03.0fdeg_t%.0fs.png',Alt,Inc,RAA,t_Res * t_Step);
        fprintf(' ... %s\n',name);
        export_fig(name,'-r150');
    end
    
    % längste Zeit ohne Sonne
    zeroPower = relPower == 0;
    % Indices mit 0
    idx0 = find(zeroPower);
    % Differenz nicht-konsekutiver Indices
    sum0 = [0 cumsum(diff(idx0)~=1)];
    % Indices von längster Sequenz
    if (numel(sum0) > 1)
        seq0 = idx0(sum0==mode(sum0));
        % Zeitspanne
        noPower = t_Res * numel(seq0);
        fprintf('\n\tLängste Zeit ohne Sonnenlicht: %.1f min,\n\tMissionzeit: %.3f -> %.3f d.\n',noPower/60,t_Vec(seq0(1)),t_Vec(seq0(end)));
    else
        fprintf('\n\tImmer Sonnenlicht vorhanden.\n');        
    end
end
fprintf(' fertich.\n');
    
%% Exportiere Daten

if (f_GenerateTFile == 1)
    fprintf('Schreibe Daten ');
    % Kopfzeile schreiben
    fid = fopen(d_TEx,'W');
    fprintf(fid, '%%# Temperaturdaten erzeugt am %s\n',date);
    fprintf(fid, '%%# Inc = %.2f deg, Alt = %.2f km, RAAN = %.2f deg, dt = %.0f s\n',Inc, Alt, RAA, t_Res * t_Step);
    % Komponenten listen
    fprintf(fid, '%%# %d: Zeit\n',0);
    count = 1;
    for ss = f_IncludedParts
        fprintf(fid, '%% %d: %s\n',count,Sat_Struct(ss).name);
        count = count + 1;
    end
    fprintf(fid, '%%# %% %% DATEN %% %% %%\n');
    fprintf(fid, '%%# kalter Fall [K]: warmer Fall [K]: Durchschnitt [K]\n');
    % Daten schreiben
    count = 0;
    for tt = ran
        % Fortschritt
        if (mod(count,floor(numel(ran)/10)) == 0)
            fprintf('.');
        end
        count = count + 1;
        % Daten
        fprintf(fid, fstrDateIn, dat_Time{1}(tt),dat_Time{2}{tt},dat_Time{3}(tt),dat_Time{4}{tt});
        for ss = f_IncludedParts
            % kalter Fall:warmer Fall:Durchschnitt
            fprintf(fid, ';%.1f:%.1f:%.1f',T(1,ss,tt),T(2,ss,tt),mean(T(:,ss,tt)));
        end
        fprintf(fid,'\n');
    end
    fclose(fid);
    fprintf(' fertich.\n');
end

%% Speichere Daten
if (f_SaveResults == 1)
    pfad = [d_Par 'Daten' fmt_Base '.mat'];
    save(pfad);
end

%% Elapsed Time
if measureTime == 1
    toc
end