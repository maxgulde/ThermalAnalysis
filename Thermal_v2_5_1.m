% Nanosatellite Thermal Analysis Script v2.5.1
%
% Last updated: 2018-07-04

% Assumptions
% - PENDING

% To-do
% - TO-DO

% To optimize
% - EarthAngles into LocalZenith (vector components instead of angles)

% Changes
% 2.5.1
% - Added separated surfaces for each part
% - Optimized surface calculation by avoiding loop repetition
% 2.5
% - Cleaner code
% - English version, translated from German
% 2.4
% - Changed albedo to fitted function of Lockheed Martin's Data (includes
%   view factor)
% - Added view factors for internal radiation
% 2.3
% - kompatibel gemacht zun Sidos Berechnungen f�r Validierung via Thermal
%   Desktop
% - Sol_Flux (hot case), fixer Albedowert, Celsius als Einheit, keine
%   Zugriffszeiten, 
% - Plotunterschiede zwischen hot und cold
% - Fix: Albedokorrektur Subsolarwinkel
% - Fix: Maximaler Albedowinkel eingef�hrt
% - Fix: Zeitkonstante korrigiert f�r thermal coupling

%% Settings
measureTime = 1;
if measureTime == 1
    tic
end

% % % Orbit
Inc = 98;   % [deg]
Alt = 700;  % [km]
RAA = 10;    % [deg]

% % % Satellite
satName = 'ERNST';          % Name of the satellite
T_Start = 273.15+52.7;              % [K] Initial temperature of the satellite
% IGNORED CODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sat_RadEffArea = 1.36;        % effektive Fl�che des Radiators, Pyramide
% Sat_RadName = 'Radiator';
% Sat_CellEff = 0.34;             % Effizienz Solarzellen
% Sat_CellName = 'SolarCells';
% IGNORED CODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % % Simulation
t_Res = 60;             % [s] Temporal resolution                          % CHECK
t_Range = [0 1];    % Simulated time [start end], full file = [0 1]
t_Step = 1;             % Step width (positive int)                        %%%%% MIGHT NOT WORK PROPERLY

% % % Options (activate / deactivate effects)
f_Sun = 1;           % Direct sunlight
f_Cmp = 0;          % Component-generated heat
f_Alb = 1;           % Albedo
f_EIR = 1;           % Earth IR
f_Emi = 1;           % IR emissions (satellite)
f_TCo = 0;          % Thermal coupling
f_IncludedParts = 0;    % Indices of simulated parts, 0 = all

% % % Input
f_ReloadAllData = 1;         % Load all data (not slow in this version) % CHECK
f_ReloadMatData = 1;         % Load material data only (fast)           % CHECK
f_UseAccessIntervals = 0;    % Calculate energy load by comms          % DISABLED
f_UseFixedAlbedo = -1;          % Set to < 0 to use correction table
f_UseCelsius = 1;            % Use Celsius at display                   % CHECK

% % % Output
f_PlotLineWidth = 2;

% % % Presentation
f_FigNum = 1;           % Figure number
f_DrawOnTop = 0;    % Keep previous plot
f_PlotGenPower = 0; % Plot power generated by satellite
f_DrawCaseIndex = 1;    % 1 = mean, 2 = hot, 3 = cold, 4 = both
f_Verbose = 0;      % Print verbose to the console
f_TemperatureLimits = [253 273 293 308];    % Temperature limits regions
f_TemperatureLimitsC = 'bgr';               % Colors for the regions
t_PlotScale = 24*3600;  % [s] Temporal scale of the plot                   

% % % Paths
d_Bas = sprintf('Data_%s_i%i_a%i_r%i_t%i',satName,Inc,Alt,RAA,t_Res); % Base path
d_StrFolder = fullfile(d_Bas,'Structure');

% % % Data Files
d_SubSol = fullfile(d_Bas,'SunAngles.csv');
d_EarthAngles = fullfile(d_Bas,'EarthAngles.csv');
d_AreaS = fullfile(d_Bas,'Out_AreaSunView.txt');
d_Mat = fullfile(d_Bas,'_materials.txt');        % Materials file
d_Str = fullfile(d_StrFolder,'_structure.txt');      % Structure file
d_TCo = fullfile(d_Bas,'_conduction.txt');    % Thermal conductivity
d_IntRad = fullfile(d_Bas,'_internal_matrix.txt'); % Internal viewfactors

% % % Format Strings
fstrMat = '%s %f %f %f'; % Name Abs Emi HCap
fstrStr = '%s %s %s %f %f %s'; % Name Optical Bulk TotArea Weight Component
fstrStrSurf = '%s %s %f %f %f %f %d'; % SurfName Optical Area Nx Ny Nz in_flag
fstrSubSol = '%d %s %d %12s,%f,%f,%f'; % D(d) M(s) Y(d) t(12s), Az(f), El(f), SubSol(f)
fStrEarthAngles = '%d %s %d %12s,%f,%f'; % D(d) M(s) Y(d) t(12s), Az(f), El(f)
fstrTCo = '%s %s %f %f %f';
fstrDateIn = '%d %s %d %s';
fstrDate = 'dd mmm yyyy HH:MM:SS';
fstrOrder = '%s %s %s %s';
%fstrDateOut = '%.0f %s %.0f %02.0f:%02.0f:%02.0f';                        % DISABLED
%fstrAccess = '%f,%f %s %f %f:%f:%f,%f %s %f %f:%f:%f,%f';                 % DISABLED
% IGNORED CODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fstrPower = '%d %s %d %12s;%f';
% fstrCmp = '%s %d %f %f';
% fstrXLo = '%s %f %f %f';
% IGNORED CODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

% Maximal Albedo Angle
AlbMaxAngle = 100; % [deg] cos(0.9*theta) <-- Up to 100 deg not in shadow

%% Data reading

if f_ReloadAllData == 1
    fprintf('Reading data:\n');
    
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
    
%     % % Access Intervals
%     if f_UseAccessIntervals
%         fprintf(' ... Access Intervals ...\n');
%         dat_temp = ReadCSV(d_Access,fstrAccess,1);
%         % Translation in seconds since start
%         dat_Access = zeros(numel(dat_temp{1}),2); % Start [s] und Ende [s]
%         startDate = sprintf(fstrDateIn,dat_Time{1}(1),dat_Time{2}{1},dat_Time{3}(1),dat_Time{4}{1}(1:end-4));
%         endDate = sprintf(fstrDateIn,dat_Time{1}(end),dat_Time{2}{size(dat_Time{1},1)},dat_Time{3}(end),dat_Time{4}{size(dat_Time{1},1)}(1:end-4));
%         startTime = datevec(startDate,fstrDate);
%         for ii = 1:numel(dat_temp{1})
%             t = datevec(sprintf(fstrDateOut,dat_temp{2}(ii),dat_temp{3}{ii},dat_temp{4}(ii),dat_temp{5}(ii),dat_temp{6}(ii),dat_temp{7}(ii)),fstrDate);
%             dat_Access(ii,1) = etime(t,startTime);
%             dat_Access(ii,2) = dat_Access(ii,1) + round(dat_temp{14}(ii));
%         end
%         ACtot = sum(dat_temp{:,14});
%     end
end

if f_ReloadAllData || f_ReloadMatData
    Sat_Mat = struct();     % Materials, Absorption and Emission
    % IGNORED CODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Sat_Cmp = struct();         % Components                              % CHECK
    % IGNORED CODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Sat_Struct = struct();                                                 % PROPER COMMENT NEEDED Oberfl�chen, Verkn�pfung mit Mat und Cmp, Fl�chen und Gewicht
    
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
                    tempStruct(iii).name = dat_temp2{1}{iii};
                    tempStruct(iii).optical = dat_temp2{2}{iii};
                    tempStruct(iii).area = dat_temp2{3}(iii);
                    tempStruct(iii).normalV = ...
                        [dat_temp2{4}(iii), dat_temp2{5}(iii), dat_temp2{6}(iii)];
                    tempStruct(iii).internal = dat_temp2{7}(iii);
                end
                Sat_Struct.surfaces = tempStruct;
            end
        else
            Sat_Struct.surfaces = struct([]);
        end
    end
    
    % % Components
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CHECK / TO-DO
    
    % % Materials
    fprintf(' ... Material data ...\n');
    fid = fopen(d_Mat);
    dat_temp = textscan(fid,fstrMat,'CommentStyle','%');
    fclose(fid);
    for ii = 1:numel(dat_temp{1})
        Sat_Mat(ii).name = dat_temp{1}{ii};
        Sat_Mat(ii).abs = dat_temp{2}(ii);
        Sat_Mat(ii).emi = dat_temp{3}(ii);
        Sat_Mat(ii).cap = dat_temp{4}(ii);
    end
    
    % % Thermal Coupling of the surfaces
    if f_TCo == 1
        fprintf(' ... Thermal Coupling ...\n');
        fid = fopen(d_TCo);
        dat_temp = textscan(fid,fstrTCo,'CommentStyle','%');
        fclose(fid);
        % 1 Create coupling matrix
        Sat_TCo = cell(numel(Sat_Struct));
        % 2 Fill with standard struct
        s = struct();
        s.Mode = 0;
        for ii = 1:numel(Sat_Struct)
            for jj = 1:numel(Sat_Struct)
                Sat_TCo{ii,jj} = s;
            end
        end
        % 3 Fill with values from the list
        for ii = 1:numel(dat_temp{1})
            % Find indices for both surfaces
            idx1 = find(strcmp({Sat_Struct.name}, dat_temp{1}{ii}));
            idx2 = find(strcmp({Sat_Struct.name}, dat_temp{2}{ii}));
            % Input data
            s = struct();
            s.Mode = 1;
            s.Length = dat_temp{3}(ii);
            s.Area = dat_temp{4}(ii);
            s.Cond = dat_temp{5}(ii);
            % Save
            Sat_TCo{idx1,idx2} = s;
            Sat_TCo{idx2,idx1} = s;
        end
    end
    
    % % Internal Radiation                                                 % CHECK AND PROPERLY INCORPORATE
    fprintf(' ... Internal Radiation TESTING ...\n');
    fid = fopen(d_IntRad);
    dat_temp = textscan(fid,'%f %f %f %f %f %f %f %f %f','CommentStyle','%');
    fclose(fid);
    internal_vf = zeros(9);
    for iii = 1:9
        internal_vf(:,iii) = dat_temp{1,iii};
    end
    
    % IGNORED CODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % extra Lastf�lle
    %fprintf(' ... Lastf�lle ...\n');
    %fid = fopen(d_XLo);
    %dat_XLo = textscan(fid, fstrXLo, 'CommentStyle','%');
    %fclose(fid);
    % IGNORED CODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
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

%% Time vector

t_Vec = (1:size(dat_Time{1},1)) * t_Res / t_PlotScale;
N = numel(t_Vec);

% Only take into account the desired simulated time
t_sim = floor(t_Range*N);
if (t_sim(1) < 1)
    %t_sim(1) = t_Step;                                                    % CHECK IF THIS WOULD BE BETTER
    t_sim(1) = 1;
end
if (t_sim(end) > N)
    t_sim(end) = N;
end

% Simulated time index vector, t_Vec trimming
t_sim_index = t_sim(1):t_Step:t_sim(end);
t_Vec = t_Vec(t_sim_index);

% IGNORED CODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% % CHCECK
% % Startzeit
% t_Offset = t_Vec(1) * t_PlotScale;
% IGNORED CODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Simulation

fprintf('Simulation ');

% Temperature for cold (1) and hot (2) cases, each part, each time-step
T = zeros(2,numel(Sat_Struct),N+1);

% Initial temperature
T(:,:,1) = T_Start;

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
modDivisor = floor(numel(t_sim_index)/10);
for tt = t_sim_index
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
    
    % % % (A) Isolated Consideration
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
            % Index der Komponente in Fl�chendatei
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
            
            %             % geringerer Energieeingang f�r Solarzellen
            %             if (strcmp(Sat_Struct(ss).name(1:end-1),Sat_CellName) == 1)
            %                 A = A * (1 - Sat_CellEff);
            %             end

            ss_abs = Sat_Mat(sIdx).abs;
            
            P_C = P_C + A * ss_abs * Sol_Flux(1);
            P_H = P_H + A * ss_abs * Sol_Flux(2);
        end
        
        % % % (2) Components (Surface independent)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                 % TO-DO
        
        % % % (3) IR Emission (internal only, optimized by _internal_matrix)
        if f_Emi == 1
            % IGNORED CODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % INCORPORATE PROPER IMPLEMENTATION
            % Update area to effective area in case of modified radiator
            %if (strcmp(Sat_Struct(ss).name,Sat_RadName) == 1)
            %    A = A * Sat_RadEffArea;
            %end
            % IGNORED CODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            for pp = f_IncludedParts(find(f_IncludedParts == ss):end)
                if ss == pp % Object is itself, calculation would be 0, skip
                    continue;
                end
                % Object is another volume, use view factor
                pIdx = find(strcmp({Sat_Mat(:).name}',Sat_Struct(pp).optical));
                pp_emi = Sat_Mat(pIdx).emi;
                pp_A = Sat_Struct(pp).size; % This is the total surf. area
                
                % Internal view factor (source, target)
                % Radiated power                                           % CHECK PROPER AREA VALUES (SHOULD BE FULL VOLUME, ADD EXTERNAL SURFACE?)
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
        
        % Prevent negative temperatures
        for ii = 1:2
            if (T(ii,ss,tt+t_Step) < T_Space)
                T(ii,ss,tt+t_Step) = T_Space;
            end
        end 
    end
    
    % % % (B) Conduction, etc. between individual components
    %%%%%%%%%%%%%%%%%%%%%%%%% TEMP CODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
                % Zeitskala absch�tzen, auf der Temperaturausgleich stattfindet
                R = Sat_TCo{ss,pp}.Length / (Sat_TCo{ss,pp}.Area * Sat_TCo{ss,pp}.Cond);
                tau = R * mean([C_s C_p]);
                t_Ratio = dt / tau;
                % hei�er und kalter Fall
                for cc = 1:2
                    % schneller Ausgleich -> beide Strukturen haben dieselbe Temperatur
                    if (t_Ratio > t_IntLim(2))
                        T_prime = (C_s * T(cc,ss,tt+t_Step) + C_p * T(cc,pp,tt+t_Step))/(C_s + C_p);
                        T(cc,ss,tt+t_Step) = T_prime;
                        T(cc,pp,tt+t_Step) = T_prime;
                    elseif (t_Ratio < t_IntLim(1)) % langsamer Ausgleich -> lineare N�herung
                        deltaT = T(cc,pp,tt+t_Step) - T(cc,ss,tt+t_Step);
                        T(cc,ss,tt+t_Step) = T(cc,ss,tt+t_Step) + dt/C_s*deltaT/R;
                        T(cc,pp,tt+t_Step) = T(cc,pp,tt+t_Step) - dt/C_p*deltaT/R;
                    else % ungef�hr gleiche Zeitskala, runterbrechen
                        dt2 = t_IntLim(1) * dt;
                        for ii = dt2:dt2:dt
                            deltaT = T(cc,pp,tt+t_Step) - T(cc,ss,tt+t_Step);
                            T(cc,ss,tt+t_Step) = T(cc,ss,tt+t_Step) + dt2/C_s*deltaT/R;
                            T(cc,pp,tt+t_Step) = T(cc,pp,tt+t_Step) - dt2/C_p*deltaT/R;
                        end                        
                    end
                    % Verh�ltnis von dt zu tau angeben
                    if (f_Verbose == 1 && tt == ran(1) && cc == 1)
                        fprintf('\nW�rmeleitung: <%s>, <%s> (dt/tau = %f).',Sat_Struct(ss).name, Sat_Struct(pp).name, t_Res/tau);
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
    %%%%%%%%%%%%%%%%%%%%%%%%% TEMP CODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
fprintf('\n ... done.\n');

%% Plots of the results

fprintf('Making plots ...');
f_DrawParts = 1;

% Set figure
figure(f_FigNum);
if f_DrawOnTop == 1
    hold on;
end

% Set temperature offset
if f_UseCelsius == 1
    tempOffset = 273.15;
else
    tempOffset = 0;
end

% Draw according to case
c = 1;
for ii = f_IncludedParts(f_DrawParts)
    switch f_DrawCaseIndex
        case 1 % Mean case
            leg = cell(numel(f_DrawParts),1);
            
            y = squeeze(mean(T(:,ii,t_sim_index(1) + 1:end),1));
            y = y(y ~= 0);
            plot(t_Vec, y - tempOffset,'LineWidth',f_PlotLineWidth);
            if (ii == f_IncludedParts(f_DrawParts(1)))
                hold on;
            end
            % Legend
            leg{c} = [Sat_Struct(ii).name ' (mean)'];
            c = c + 1;
            
        case 2 % Hot case
            leg = cell(numel(f_DrawParts),1);
            
            y = squeeze(T(:,ii,t_sim_index(1) + 1:end));
            plot(t_Vec, y(2,(y(2,:) ~= 0)) - tempOffset,'LineWidth',f_PlotLineWidth);
            if (ii == f_IncludedParts(f_DrawParts(1)))
                hold on;
            end
            % Legend
            leg{c} = [Sat_Struct(ii).name ' (hot)'];
            c = c + 1;
            
        case 3 % Cold case
            leg = cell(numel(f_DrawParts),1);
            
            y = squeeze(T(:,ii,t_sim_index(1) + 1:end));
            plot(t_Vec, y(1,(y(1,:) ~= 0)) - tempOffset,'LineWidth',f_PlotLineWidth);
            if (ii == f_IncludedParts(f_DrawParts(1)))
                hold on;
            end
            % Legend
            leg{c} = [Sat_Struct(ii).name ' (cold)'];
            c = c + 1;
            
        case 4 % Both cases
            leg = cell(numel(f_DrawParts)*2,1);
            
            y = squeeze(T(:,ii,t_sim_index(1) + 1:end));
            plot(t_Vec, y(1,(y(1,:) ~= 0)) - tempOffset,'LineWidth',f_PlotLineWidth);
            if (ii == f_IncludedParts(f_DrawParts(1)))
                hold on;
            end
            plot(t_Vec, y(2,(y(2,:) ~= 0)) - tempOffset,'LineWidth',f_PlotLineWidth);
            % Legends
            leg{c} = [Sat_Struct(ii).name ' (cold)'];
            c = c + 1;
            leg{c} = [Sat_Struct(ii).name ' (hot)'];
            c = c + 1;
            
        otherwise
            fprintf('\nError making plot, draw case mismatch.\n');
    end
end
% Temperature limits
for tt = 2:numel(f_TemperatureLimits)
    x = xlim;
    y = f_TemperatureLimits;
    p = patch([x(1) x(2) x(2) x(1)],[y(tt-1) y(tt-1) y(tt) y(tt)] - tempOffset,f_TemperatureLimitsC(tt-1));
    set(p,'FaceAlpha',0.05,'EdgeAlpha',0);
end
hold off;

% Title and axis
title(sprintf('inc = %.0f deg, alt = %.0f km, raan = %.0f deg, temp res = %.0f s',Inc,Alt,RAA,t_Res));
xlabel 'time [d]'
if f_UseCelsius == 1
    ylabel 'temperature [�C]'
else
    ylabel 'temperature [K]'
end
legend(leg);
set(gcf,'Color','w');
set(gca,'FontName','Calibri','FontSize',18);

% IGNORED CODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % exportiere Plot
% if (f_ExportPlot == 1)
%     fprintf('Exportiere ...\n');
%     figure(f_FigNum);
%     name = sprintf('thermal_a%.0fkm_i%.0fdeg_r%03.0fdeg_t%.0fs.png',Alt,Inc,RAA,t_Res * t_Step);
%     fprintf(' ... %s\n',name);
%     export_fig(name,'-r150');
% end
% 
% % Plot der generierten Leistung
% if (f_PlotGenPower == 1)
%     figure(f_FigNum + 1);
%     if (f_DrawOnTop == 1)
%         hold on;
%     end
%     % plot solar panel power generation
%     relPower = dat_Power{1}';
%     relPower = relPower(ran);
% %     relPower = movingmean(relPower,t_DaySec/(t_Res * t_Step),2,[]);
% %     relPower = smooth(relPower,10);
%     plot(t_Vec, relPower);
%     xlim([t_Vec(1) t_Vec(end)]);
%     xlabel 'time [d]'
%     ylabel 'daily mean solar power [W]'
%     hold off;
%     
%     % exportiere Plot
%     if (f_ExportPlot == 1)
%         fprintf('Exportiere ...\n');
%         figure(f_FigNum + 1);
%         name = sprintf('power_a%.0fkm_i%.0fdeg_r%03.0fdeg_t%.0fs.png',Alt,Inc,RAA,t_Res * t_Step);
%         fprintf(' ... %s\n',name);
%         export_fig(name,'-r150');
%     end
%     
%     % l�ngste Zeit ohne Sonne
%     zeroPower = relPower == 0;
%     % Indices mit 0
%     idx0 = find(zeroPower);
%     % Differenz nicht-konsekutiver Indices
%     sum0 = [0 cumsum(diff(idx0)~=1)];
%     % Indices von l�ngster Sequenz
%     if (numel(sum0) > 1)
%         seq0 = idx0(sum0==mode(sum0));
%         % Zeitspanne
%         noPower = t_Res * numel(seq0);
%         fprintf('\n\tL�ngste Zeit ohne Sonnenlicht: %.1f min,\n\tMissionzeit: %.3f -> %.3f d.\n',noPower/60,t_Vec(seq0(1)),t_Vec(seq0(end)));
%     else
%         fprintf('\n\tImmer Sonnenlicht vorhanden.\n');        
%     end
% end
% IGNORED CODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(' done.\n');

% IGNORED CODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Exportiere Daten
% 
% if (f_GenerateTFile == 1)
%     fprintf('Schreibe Daten ');
%     % Kopfzeile schreiben
%     fid = fopen(d_TEx,'W');
%     fprintf(fid, '%%# Temperaturdaten erzeugt am %s\n',date);
%     fprintf(fid, '%%# Inc = %.2f deg, Alt = %.2f km, RAAN = %.2f deg, dt = %.0f s\n',Inc, Alt, RAA, t_Res * t_Step);
%     % Komponenten listen
%     fprintf(fid, '%%# %d: Zeit\n',0);
%     count = 1;
%     for ss = f_IncludedParts
%         fprintf(fid, '%% %d: %s\n',count,Sat_Struct(ss).name);
%         count = count + 1;
%     end
%     fprintf(fid, '%%# %% %% DATEN %% %% %%\n');
%     fprintf(fid, '%%# kalter Fall [K]: warmer Fall [K]: Durchschnitt [K]\n');
%     % Daten schreiben
%     count = 0;
%     for tt = ran
%         % Fortschritt
%         if (mod(count,floor(numel(ran)/10)) == 0)
%             fprintf('.');
%         end
%         count = count + 1;
%         % Daten
%         fprintf(fid, fstrDateIn, dat_Time{1}(tt),dat_Time{2}{tt},dat_Time{3}(tt),dat_Time{4}{tt});
%         for ss = f_IncludedParts
%             % kalter Fall:warmer Fall:Durchschnitt
%             fprintf(fid, ';%.1f:%.1f:%.1f',T(1,ss,tt),T(2,ss,tt),mean(T(:,ss,tt)));
%         end
%         fprintf(fid,'\n');
%     end
%     fclose(fid);
%     fprintf(' fertich.\n');
% end
% 
% %% Speichere Daten
% if (f_SaveResults == 1)
%     pfad = [d_Par 'Daten' fmt_Base '.mat'];
%     save(pfad);
% end
% IGNORED CODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if measureTime == 1
    toc
end