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
% 2.3
% - kompatibel gemacht zun Sidos Berechnungen für Validierung via Thermal Desktop
% - Sol_Flux (hot case), fixer Albedowert, Celsius als Einheit, keine Zugriffszeiten, 
% - Plotunterschiede zwischen hot und cold
% - Fix: Albedokorrektur Subsolarwinkel
% - Fix: Maximaler Albedowinkel eingeführt
% - Fix: Zeitkonstante korrigiert für thermal coupling

% clc;
% clear;

%% Einstellungen

% % % Orbit
Inc = 97;
Alt = 500;
RAA = 0;

% % % Satellit
satName = 'ERNST';          % Name des Satelliten
T_Start = 293;                  % Starttemperatur des Satelliten
Sat_RadEffArea = 1.36;          % effektive Fläche des Radiators, Pyramide
Sat_RadName = 'Radiator';
Sat_CellEff = 0.34;             % Effizienz Solarzellen
Sat_CellName = 'SolarCells';
Sat_SurfNum = 21;               % Anzahl der Oberflächen wie in Simulation berechnet

% % % Simulation
t_Res = 120;                    % [s] zeitliche Auflösung
t_ResAcc = t_Res/6;             % [s] zeitliche Auflösung Zugriff, nur wenn f_UseMeanLoads == 0
t_Range = [0 2/365] + 0.0;      % simulierte Zeit [start ende], [0 1] voll
t_Step = 1;                     % Schrittweite
t_IntLim = [1/60 5];            % Grenzen der Zeitkonstanten für die Simulation der thermischen Kopplung

% % % Paths %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Testing simple cube sat
d_Bas = 'Data\';
d_Dat = [d_Bas 'ERNST_i97_a500_r00_t120_stiff\'];
d_Par = 'Data\';
d_Cmp = [d_Par '_Komponenten.txt'];
d_Mat = [d_Par 'new_materials.txt'];
d_Sur = [d_Par 'new_structure.txt'];
d_TCo = [d_Par 'new_Waermeleitung.txt'];
d_XLo = [d_Par '_XLoads.txt'];
d_TEx = [d_Par 'Temperatur'];
% d_Suff = ' - short';
d_Suff = '';

% % % Optionen zum Zu- und Abschalten bestimmter Effekt
f_Sun = 1;      % Sonneneinstrahlung
f_Cmp = 0;      % Abwärme Komponenten
f_Alb = 1;      % Albedo
f_EIR = 1;      % Erde IR
f_Emi = 1;      % Emission Oberflächen
f_TCo = 0;      % Thermische Kopplung
f_XLo = 0;      % extra Loads, die nur kurzzeitig anfallen
f_IncludedParts = 0;    % Indices simulierte Strukturteile, 0 = alle
% f_DrawParts = [1:10 11];

% % % Eingabe
f_ReloadAllData = 1;            % alle Daten (Winkel, Leistung, etc) einladen (langsam)
f_ReloadMatData = 1;            % Materialdaten einladen (schnell)
f_IgnoreAccessIntervals = 0;    % ignore increased energy load by communication devices
f_UseFixedAlbedo = -1;          % set to < 0 to use correction table
f_UseCelsius = 273.15;           % set to == 0 to use Kelvin

% % % Darstellung
f_FigNum = 2;           % Nummer der figure
f_DrawOnTop = 0;        % drüberzeichnen
f_PlotGenPower = 0;
f_DrawCaseIndex = 1;    % 1: mean, 2: hot, 3: cold, 4: both
f_Verbose = 0;
f_UseMeanLoads = 0;
f_TemperatureLimits = [253 273 293 308];
f_TemperatureLimitsC = 'bgr';

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

% % % Daten aus "SC Thermal Control Handbook", p. 28 (24 h, (C)old und (H)ot case)
% Direktes Sonnenlicht
Sol_Flux = [1322 1414];     % [W/m^2] cold and hot case
T_Space = 2.7;              % [K] Temperatur Weltraum

% Albedo
Alb_KorrSubsolar = [0:10:90; 0:0.01:0.05 0.08 0.13 0.2 0.31]';
Alb_OrbitInc_C = [30 60 90; 0.17 0.2 0.2]';
Alb_OrbitInc_H = [30 60 90; 0.19 0.23 0.23]';
% Erde IR
EIR_OrbitInc_C = [30 60 90; 236 226 225]';  % [W/m^2]
EIR_OrbitInc_H = [30 60 90; 257 241 230]';

% Daten
fmt_Base = sprintf('_i%.0f_a%.0f_r%02.0f_t%03.0f',Inc,Alt,RAA,t_Res);
d_TEx = [d_TEx fmt_Base '.txt'];
d_AreaS = sprintf('%sOut_AreaSunView%s%s.txt',d_Dat,fmt_Base,d_Suff);
d_AreaE = sprintf('%sOut_AreaEarthView%s%s.txt',d_Dat,fmt_Base,d_Suff);
d_Power = sprintf('%sOut_Power%s%s.txt',d_Dat,fmt_Base,d_Suff);
d_Access = sprintf('%s%s%s Access.csv',d_Dat,satName,fmt_Base(1:end-5));
d_SubSol = sprintf('%sSunAngles.csv',d_Dat);
d_EarthAngles = sprintf('%sEarthAngles.csv',d_Dat);

% Formatstrings
fstrAreas = ['%d %s %d %12s' repmat(';%f',1,Sat_SurfNum)];
fstrPower = '%d %s %d %12s;%f';
fstrAccess = '%f,%f %s %f %f:%f:%f,%f %s %f %f:%f:%f,%f';
fstrDateIn = '%.0f %s %.0f %s';
fstrDateOut = '%.0f %s %.0f %02.0f:%02.0f:%02.0f';
fstrDate = 'dd mmm yyyy HH:MM:SS';
fstrCmp = '%s %d %f %f';
fstrMat = '%s %f %f %f';
fstrStr = '%s %s %s %f %f %s %f %f';
fstrSubSol = '%d %s %d %12s,%f,%f,%f'; % day (d) month (s) year (d) time (12s), azimuthal (f), elevation (f), subsolar (f)
fStrEarthAngles = '%d %s %d %12s,%f,%f'; % day (d) month (s) year (d) time (12s), azimuthal (f), elevation (f)
fstrOrder = '%s %s %s %s';
fstrTCo = '%s %s %f %f %f';
fstrXLo = '%s %f %f %f';
% Header Zeilen
h_Area = Sat_SurfNum + 2;
h_Power = 3;

% Zeit
t_PlotScale = 24*3600;
t_DaySec = 24*3600;

% Physikalische Konstanten
consts;
clearConstants = true;
re = 6371000;
r = (re)/(re+700*1e3); % Testing
%r = (RE)/(RE+Alt*1e3);

% maximaler Subsolarwinkel für Albedo
AlbMaxAngle = 100; % [deg] (Above 100 deg the formula generates complex numbers)

%% Daten einlesen
if (f_ReloadAllData == 1)
    fprintf('Einlesen der Daten:\n');
    
    % % Flächen und Leistung
    % Reihenfolge extrahieren
    fid = fopen(d_AreaS);
    dat_temp = textscan(fid, fstrOrder, h_Area-2, 'HeaderLines', 2);
    dat_Order = dat_temp(3);
    fclose(fid);
    
    % Flächen extrahieren
    fprintf(' ... Flächen aus Sonnensicht ...\n');
    dat_temp = ReadCSV(d_AreaS,fstrAreas,h_Area);
    dat_AreaS = dat_temp(5:end);
    dat_Time = dat_temp(1:4);
    fprintf(' ... Flächen aus Erdsicht ...\n');
    dat_temp = ReadCSV(d_AreaE,fstrAreas,h_Area);
    dat_AreaE = dat_temp(5:end);
    if (f_PlotGenPower == 1)
        fprintf(' ... Erzeugte Leistung ...\n');
        dat_temp = ReadCSV(d_Power,fstrPower,h_Power);
        dat_Power = dat_temp(5:end);
    end
    
    % % Subsolarwinkel
    fprintf(' ... Subsolarwinkel ...\n');
    dat_temp = ReadCSV(d_SubSol,fstrSubSol,1);
    dat_SolAng = dat_temp(end-2:end-1); % Azimuth, Elev
    dat_SubSol = dat_temp(end);
    
    % % Local Zenith
    fprintf(' ... Local Zenith ...\n');
    dat_temp = ReadCSV(d_EarthAngles,fStrEarthAngles,1);
    dat_EarthAngles = dat_temp(end-1:end);
    
    
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
    else
        
    end

end

if (f_ReloadMatData == 1 || f_ReloadAllData == 1)
    
    % % Flächen: Name, Material, Größe, Komponenten
    fprintf(' ... Struktur ...\n');
    fid = fopen(d_Sur);
    dat = textscan(fid, fstrStr, 'CommentStyle','%');
    fclose(fid);
    StructNum = numel(dat{1});
    % Anzahl der simulierten Teile setzen, falls noch nicht geschehen
    if (f_IncludedParts == 0)
        f_IncludedParts = 1:StructNum;
    end
    for ii = 1:StructNum
        Sat_Struct(ii).name = dat{1}{ii};
        Sat_Struct(ii).surf = dat{2}{ii};   % components have NA surface
        Sat_Struct(ii).bulk = dat{3}{ii};
        Sat_Struct(ii).size = dat{4}(ii);
        Sat_Struct(ii).mass = dat{5}(ii);
        Sat_Struct(ii).cmp = strsplit(dat{6}{ii},',');
        Sat_Struct(ii).azimuth = dat{7}(ii);
        Sat_Struct(ii).elevation = dat{8}(ii);
        % Position in Flächendatei
        Sat_Struct(ii).AFileIdx = find(strcmp(dat_Order{:},Sat_Struct(ii).name));
        if (isempty(Sat_Struct(ii).AFileIdx))
            fprintf('\t\t+Interne Struktur: %s \n',Sat_Struct(ii).name);
        end
    end

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
        %         fprintf(' ... Thermische Kopplung ...\n');
        %         fid = fopen(d_TCo);
        %         dat = textscan(fid, fstrTCo, 'CommentStyle','%');
        %         fclose(fid);
        %         % 1 Kopplungsmatrix erstellen
        %         Sat_TCo = cell(numel(Sat_Struct));
        %         % 2 Standard Struktur reinschreiben
        %         s = struct();
        %         s.Mode = 0;
        %         for ii = 1:numel(Sat_Struct)
        %             for jj = 1:numel(Sat_Struct)
        %                 Sat_TCo{ii,jj} = s;
        %             end
        %         end
        %         % 3 Mit Werten aus Liste füllen
        %         for ii = 1:numel(dat{1})
        %             % Indices finden für beide Oberflächen
        %             idx1 = find(strcmp({Sat_Struct.name}, dat{1}{ii}));
        %             idx2 = find(strcmp({Sat_Struct.name}, dat{2}{ii}));
        %             % Daten eintragen
        %             s = struct();
        %             s.Mode = 1;
        %             s.Length = dat{3}(ii);
        %             s.Area = dat{4}(ii);
        %             s.Cond = dat{5}(ii);
        %             % Speichern
        %             Sat_TCo{idx1,idx2} = s;
        %             Sat_TCo{idx2,idx1} = s;
        %         end
    end
    
    fprintf(' ... Internal Radiation TESTING ...\n');
    fid = fopen([d_Par 'matrix.txt']);
    dat = textscan(fid, '%f %f %f %f %f %f %f %f %f', 'CommentStyle','%');
    fclose(fid);
    internal_vf = zeros(9);
    for iii = 1:9
        internal_vf(:,iii) = dat{1,iii};
    end
    
    % % extra Lastfälle
    fprintf(' ... Lastfälle ...\n');
    fid = fopen(d_XLo);
    dat_XLo = textscan(fid, fstrXLo, 'CommentStyle','%');
    fclose(fid);
    
    % Gesamtgewicht
    fprintf(' Gesamtgewicht:  %.2f kg ...\n',sum([Sat_Struct(:).mass]));
    % von Komponenten
    mass = 0;
    for ii = 1:numel(Sat_Struct)
        if (strcmp(Sat_Struct(ii).surf,'NA') == 1)
            mass = mass + Sat_Struct(ii).mass;
        end
    end
    fprintf(' \t Struktur:    %.2f kg ...\n',sum([Sat_Struct(:).mass]) - mass);
    fprintf(' \t Komponenten: %.2f kg ...\n',mass);
    
    % einmal komplette Struktur anzeigen
    fprintf('\tStruktur:\n');
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

% Oberflächen sortieren
f_IncludedParts = sort(f_IncludedParts);

% Zeitschritte
count = 1;
for tt = ran
    % Fortschritt
    if (mod(count,floor(numel(ran)/10)) == 0)
        fprintf('.');
    end
    % aktuelle Zeit
    t_Now = (count - 1) * t_Res * t_Step;
    count = count + 1;
    
    % Local zenith
    [localZenith(1),localZenith(2),localZenith(3)] = ...
        sph2cart(deg2rad(dat_EarthAngles{1}(tt)),deg2rad(dat_EarthAngles{2}(tt)),1);
    localZenith = - localZenith; % Negative of the vector that points towards the Earth
    
    % Sun vector
    [sunVector(1),sunVector(2),sunVector(3)] = ...
        sph2cart(deg2rad(dat_SolAng{1}(tt)),deg2rad(dat_SolAng{2}(tt)),1);
    
    % % % (A) Isolierte Betrachtung
    for ss = f_IncludedParts
        
        % Leistungszuwachs
        dP_C = 0;
        dP_H = 0;
        
        % View Factor for Albedo and Earth IR
        [normalV(1),normalV(2),normalV(3)] = ...
            sph2cart(deg2rad(Sat_Struct(ss).azimuth),deg2rad(Sat_Struct(ss).elevation),1);
        rho = rad2deg(atan2(norm(cross(normalV,localZenith)), dot(normalV,localZenith)));
        vF = viewFactor(r,rho);
        
        % Index der Komponente in Flächendatei
        cIdx = Sat_Struct(ss).AFileIdx;
        % Wenn interne Struktur, dann auf -1 setzen
        if (isempty(cIdx))
            cIdx = -1;
        end
        
        % Material finden
        sIdx = find(strcmp({Sat_Mat(:).name}',Sat_Struct(ss).surf));
        bIdx = find(strcmp({Sat_Mat(:).name}',Sat_Struct(ss).bulk));
        if (isempty(sIdx))
            fprintf('Error: Material <%s> nicht gefunden. Wechsle zu <%s>\n.',Sat_Struct(ss).surf, Sat_Mat(1).name);
            sIdx = 1;
        end
        if (isempty(bIdx))
            fprintf('Error: Material <%s> nicht gefunden. Wechsle zu <%s>\n.',Sat_Struct(ss).bulk, Sat_Mat(1).name);
            bIdx = 1;
        end
        
        % % % (1) Direktes Sonnenlicht
        %         if (f_Sun == 1)
        %             % Fläche
        %             if (cIdx >= 0)
        %                 A = dat_AreaS{1,cIdx}(tt);
        %             else
        %                 A = 0;
        %             end
        %             % geringerer Energieeingang für Solarzellen
        %             if (strcmp(Sat_Struct(ss).name(1:end-1),Sat_CellName) == 1)
        %                 A = A * (1 - Sat_CellEff);
        %             end
        %             % cold case
        %             % aufgenommene Leistung
        %             P_C = A * Sat_Mat(sIdx).abs * Sol_Flux(1);
        %             % wenn nach außen gerichtete Fläche null, dann auch aufgenommene Leistung 0
        %             P_C = P_C * (Sat_Struct(ss).size > 0);
        %             % hot case
        %             P_H = A * Sat_Mat(sIdx).abs * Sol_Flux(2);
        %             P_H = P_H * (Sat_Struct(ss).size > 0);
        %             % Zuwachs
        %             dP_C = dP_C + P_C;
        %             dP_H = dP_H + P_H;
        %         end
        
        % % % (1) Direktes Sonnenlicht (v2)
        if (f_Sun == 1)
            % Fläche
            A = Sat_Struct(ss).size;
            ss_abs = Sat_Mat(sIdx).abs;
            
            % geringerer Energieeingang für Solarzellen
            if (strcmp(Sat_Struct(ss).name(1:end-1),Sat_CellName) == 1)
                A = A * (1 - Sat_CellEff);
            end
            
            surf_sun_angle = rad2deg(atan2(norm(cross(normalV,sunVector)), dot(normalV,sunVector)));
            
            if surf_sun_angle < 90 || surf_sun_angle > -90
                dP_C = dP_C + A * ss_abs * Sol_Flux(1) * cosd(surf_sun_angle);
                dP_H = dP_H + A * ss_abs * Sol_Flux(2) * cosd(surf_sun_angle);
            end
        end
        
        % % % (2) Komponenten
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
            dP_C = dP_C + P;
            dP_H = dP_H + P;
        end
        
        %         % % % (3) IR Abstrahlung des Satelliten
        %         if (f_Emi == 1)
        %             % Fläche * Emissivität
        % %             if (cIdx >= 0)
        %                 A = Sat_Struct(ss).size; %%%%%%%%%%%%%% MAKE SURE TO IGNORE INTERNAL COMPONENTS!!!
        % %             else
        % %                 A = 0;
        % %             end
        %             % Radiator hat strukturierte Oberfläche
        %             if (strcmp(Sat_Struct(ss).name,Sat_RadName) == 1)
        %                 A = A * Sat_RadEffArea;
        %             end
        %             % Leistungsänderung
        %             dP_C = dP_C - ksb * (T(1,ss,tt)^4 - T_Space^4) * A * Sat_Mat(sIdx).emi;
        %             dP_H = dP_H - ksb * (T(2,ss,tt)^4 - T_Space^4) * A * Sat_Mat(sIdx).emi;
        %         end
        
        % % % (3) IR Emission of the Satellite (v2, including internal radiation)
        if (f_Emi == 1)
            % Get the area of the radiating (and absorbing) surface
            A = Sat_Struct(ss).size; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NEED TO UPDATE, MATRIX REPRESENTS VF OF TOTAL OBJECT
            ss_emi = Sat_Mat(sIdx).emi;
            
            % Update area to effective area in case of modified radiator
            if (strcmp(Sat_Struct(ss).name,Sat_RadName) == 1)
                A = A * Sat_RadEffArea;
            end
            
            for pp = f_IncludedParts
                if ss == 1 || ss == 2 || ss == 3
                    continue;
                end
                if ss == pp % Object is itself, radiate power
                    % Change in incoming power
                    dP_C = dP_C - ksb * T(1,ss,tt)^4 * A * ss_emi;
                    dP_H = dP_H - ksb * T(2,ss,tt)^4 * A * ss_emi;
                else % Object is another surface, use view factor
                    pIdx = find(strcmp({Sat_Mat(:).name}',Sat_Struct(pp).surf));
                    pp_emi = Sat_Mat(pIdx).emi;
                    A_src = Sat_Struct(pp).size; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NEED TO UPDATE, MATRIX REPRESENTS VF OF TOTAL OBJECT
                    
                    % Internal view factor (emitter, receiver)
                    % Change in incoming power
                    dP_C = dP_C + ksb * T(1,ss,tt)^4 * ss_emi * pp_emi * A_src * internal_vf(pp,ss);
                    dP_H = dP_H + ksb * T(2,ss,tt)^4 * ss_emi * pp_emi * A_src * internal_vf(pp,ss);
                end
            end
            
            % If external surface, radiate to the outside
            if ~isnan(Sat_Struct(ss).azimuth) % External surfaces have a non NaN azimuth angle
                dP_C = dP_C - ksb * (T(1,ss,tt)^4 - T_Space^4) * A * ss_emi;
                dP_H = dP_H - ksb * (T(2,ss,tt)^4 - T_Space^4) * A * ss_emi;
            end
        end
        
        % % % (4) Erde IR
        if (f_EIR == 1)
            % Fläche
            A = Sat_Struct(ss).size;
            % Korrektur um Inklination
            inc_c = Inc;
            if(Inc > 90)
                inc_c = 180 - Inc;
            end
            facInc_C = EIR_OrbitInc_C(find(EIR_OrbitInc_C(:,1) <= inc_c,1,'last'),2);
            facInc_H = EIR_OrbitInc_H(find(EIR_OrbitInc_H(:,1) <= inc_c,1,'last'),2);
            % Leistungsänderung (emi gibt den Wert der Emission/Absorption im IR an)
            if (isnan(Sat_Struct(ss).azimuth) || isnan(Sat_Struct(ss).elevation)) % Internal component, does not receive albedo
                P_C = 0;
                P_H = 0;
            else % No internal component, receives albedo
                P_C = A * Sat_Mat(sIdx).emi * facInc_C .* vF; %%%%%%%%%%%%%%%%%%%%%%%%%%% CHECK %%%%%%%%%%%%%%%%%%%%
                P_H = A * Sat_Mat(sIdx).emi * facInc_H .* vF; %%%%%%%%%%%%%%%%%%%%%%%%%%% CHECK %%%%%%%%%%%%%%%%%%%%
            end
            % Energieänderung
            dP_C = dP_C + P_C;
            dP_H = dP_H + P_H;
        end
        
        % % % (5) Albedo (nur wenn Subsolarwinkel kleiner als maximaler Winkel)
        if (f_Alb == 1 && dat_SubSol{1}(tt) < AlbMaxAngle)
            % Fläche
            A = Sat_Struct(ss).size;
            % Korrektur um Inklination für SSO
            inc_c = Inc;
            if(Inc > 90)
                inc_c = 180 - Inc;
            end
            if (f_UseFixedAlbedo < 0)
                %Einfluss anhand der Inklination
                facInc_C = Alb_OrbitInc_C(find(Alb_OrbitInc_C(:,1) <= inc_c,1,'last'),2);
                facInc_H = Alb_OrbitInc_H(find(Alb_OrbitInc_H(:,1) <= inc_c,1,'last'),2);
                %Korrektur anhand des Subsolarwinkels
                alpha = dat_SubSol{1}(tt);
                a_c = Alb_KorrSubsolar(find(Alb_KorrSubsolar(:,1) <= alpha,1,'last'),2);
                facInc_C = facInc_C + a_c;
                facInc_H = facInc_H + a_c;
            else
                facInc_C = f_UseFixedAlbedo;
                facInc_H = f_UseFixedAlbedo;
            end
            % Compute Albedo flux
            AlbedoFlux_C = albedoSolarFlux(Sol_Flux(1).*facInc_C,dat_SubSol{1}(tt),vF);
            AlbedoFlux_H = albedoSolarFlux(Sol_Flux(2).*facInc_H,dat_SubSol{1}(tt),vF);
            % Leistungsänderung
            if (isnan(Sat_Struct(ss).azimuth) || isnan(Sat_Struct(ss).elevation)) % Internal component, does not receive albedo
                P_C = 0;
                P_H = 0;
            else % No internal component, receives albedo
                P_C = A * Sat_Mat(sIdx).abs * AlbedoFlux_C;
                P_H = A * Sat_Mat(sIdx).abs * AlbedoFlux_H;
            end
            % Leistungszuwachs
            dP_C = dP_C + P_C;
            dP_H = dP_H + P_H;
        end
        
        
        % % % (X) Temperaturberechnung
    
        % Energie erhalten durch Skalierung Leistung mit Zeitschritt
        E_C = dP_C * (t_Res * t_Step);
        E_H = dP_H * (t_Res * t_Step);

        % Umrechnen in Temperaturänderung
        dT_C = E_C / (Sat_Mat(bIdx).cap * Sat_Struct(ss).mass);
        dT_H = E_H / (Sat_Mat(bIdx).cap * Sat_Struct(ss).mass);
        T(1,ss,tt+t_Step) = T(1,ss,tt) + dT_C;
        T(2,ss,tt+t_Step) = T(2,ss,tt) + dT_H;
        % Negative Temperaturen verhindern
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
%             sIdx = find(strcmp({Sat_Mat(:).name}',Sat_Struct(ss).surf));
%             targ_emi = Sat_Mat(sIdx).emi;
%             
%             dP_C = 0;
%             dP_H = 0;
%             for pp = f_IncludedParts
%                 % Remove self rad, add incoming rad (only taking into accound internal sides)
%                 if ss == pp
%                     dP_C = dP_C - ksb * T(1,pp,tt)^4 * Sat_Struct(pp).size * targ_emi;
%                     dP_H = dP_H - ksb * T(2,pp,tt)^4 * Sat_Struct(pp).size * targ_emi;
%                 else
%                     pIdx = find(strcmp({Sat_Mat(:).name}',Sat_Struct(pp).surf));
%                     A_emi = Sat_Struct(pp).size;
%                     source_emi = Sat_Mat(pIdx).emi;
%                     dP_C = dP_C + target_view_factor(pp,ss) * ksb * T(1,pp,tt)^4 * A_emi * source_emi * targ_emi;
%                     dP_H = dP_H + target_view_factor(pp,ss) * ksb * T(2,pp,tt)^4 * A_emi * source_emi * targ_emi;
%                 end
%             end
%             E_C = dP_C * (t_Res * t_Step);
%             E_H = dP_H * (t_Res * t_Step);
%             
%             % Umrechnen in Temperaturänderung
%             dT_C = E_C / (Sat_Mat(sIdx).cap * Sat_Struct(ss).mass);
%             dT_H = E_H / (Sat_Mat(sIdx).cap * Sat_Struct(ss).mass);
%             T(1,ss,tt+t_Step) = T(1,ss,tt+t_Step) + dT_C;
%             T(2,ss,tt+t_Step) = T(2,ss,tt+t_Step) + dT_H;
%             % Negative Temperaturen verhindern
%             for ii = 1:2
%                 if (T(ii,ss,tt+t_Step) < T_Space)
%                     T(ii,ss,tt+t_Step) = T_Space;
%                 end
%             end
            %             % finde Index, nur ab dort
            %             idx = find(f_IncludedParts == ss);
            %             % nur solche Strukturen, die auch simuliert werden
            %             for pp = f_IncludedParts(idx:end)
            %                 % nicht auf sich selbst anwenden
            %                 if (ss == pp)
            %                     continue;
            %                 end
            %
            %                 % checken, ob Kopplung vorhanden
            %                 if (Sat_TCo{ss,pp}.Mode == 0)
            %                     continue;
            %                 end
            %                 % Bulk Material Parameter Index
            %                 bIdx_s = find(strcmp({Sat_Mat(:).name}',Sat_Struct(ss).bulk));
            %                 bIdx_p = find(strcmp({Sat_Mat(:).name}',Sat_Struct(pp).bulk));
            %                 % Parameter extrahieren
            %                 m_s = Sat_Struct(ss).mass;
            %                 m_p = Sat_Struct(pp).mass;
            %                 hc_s = Sat_Mat(bIdx_s).cap;
            %                 hc_p = Sat_Mat(bIdx_p).cap;
            %                 C_s = m_s * hc_s;
            %                 C_p = m_p * hc_p;
            %                 dt = t_Res * t_Step;
            %                 % Zeitskala abschätzen, auf der Temperaturausgleich stattfindet
            %                 R = Sat_TCo{ss,pp}.Length / (Sat_TCo{ss,pp}.Area * Sat_TCo{ss,pp}.Cond);
            %                 tau = R * mean([C_s C_p]);
            %                 t_Ratio = dt / tau;
            %                 % heißer und kalter Fall
            %                 for cc = 1:2
            %                     % schneller Ausgleich -> beide Strukturen haben dieselbe Temperatur
            %                     if (t_Ratio > t_IntLim(2))
            %                         T_prime = (C_s * T(cc,ss,tt+t_Step) + C_p * T(cc,pp,tt+t_Step))/(C_s + C_p);
            %                         T(cc,ss,tt+t_Step) = T_prime;
            %                         T(cc,pp,tt+t_Step) = T_prime;
            %                     elseif (t_Ratio < t_IntLim(1)) % langsamer Ausgleich -> lineare Näherung
            %                         deltaT = T(cc,pp,tt+t_Step) - T(cc,ss,tt+t_Step);
            %                         T(cc,ss,tt+t_Step) = T(cc,ss,tt+t_Step) + dt/C_s*deltaT/R;
            %                         T(cc,pp,tt+t_Step) = T(cc,pp,tt+t_Step) - dt/C_p*deltaT/R;
            %                     else % ungefähr gleiche Zeitskala, runterbrechen
            %                         dt2 = t_IntLim(1) * dt;
            %                         for ii = dt2:dt2:dt
            %                             deltaT = T(cc,pp,tt+t_Step) - T(cc,ss,tt+t_Step);
            %                             T(cc,ss,tt+t_Step) = T(cc,ss,tt+t_Step) + dt2/C_s*deltaT/R;
            %                             T(cc,pp,tt+t_Step) = T(cc,pp,tt+t_Step) - dt2/C_p*deltaT/R;
            %                         end
            %                     end
            %                     % Verhältnis von dt zu tau angeben
            %                     if (f_Verbose == 1 && tt == ran(1) && cc == 1)
            %                         fprintf('\nWärmeleitung: <%s>, <%s> (dt/tau = %f).',Sat_Struct(ss).name, Sat_Struct(pp).name, t_Res/tau);
            %                     end
            %                 end
            %             end
        end
    end
    
end
fprintf('\n ... fertich.\n');

    
%% Plotten der Ergebnisse
    
fprintf('Ergebnisse zeichnen ...');
f_DrawParts = 1:7;

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

%% Clear workspace clutter
if clearConstants
    clearConsts
end