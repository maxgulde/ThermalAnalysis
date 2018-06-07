% thermale Analyse ERNST 2.3-te N�herung
% Autor: Max Gulde
% Stand: 2017-09-20

% Annahmen
% - Energietransport zwischen internen Komponenten �ber Strahlung vernachl�ssigbar
%   - konservativer Ansatz, Strahlungstransport entspannt zus�tzlich
% - Perfekte W�rmeleitung innerhalb eines Bauteils
%   - Zeitskala f�r Tempertatur�nderungen langsam

% todo
% - ViewFactors einbetten f�r Aussenfl�chen in ThermoSim

% �nderungen
% 2.4
% - Change #3
% 2.3
% - kompatibel gemacht zun Sidos Berechnungen f�r Validierung via Thermal Desktop
% - Sol_Flux (hot case), fixer Albedowert, Celsius als Einheit, keine Zugriffszeiten, 
% - Plotunterschiede zwischen hot und cold
% - Fix: Albedokorrektur Subsolarwinkel
% - Fix: Maximaler Albedowinkel eingef�hrt
% - Fix: Zeitkonstante korrigiert f�r thermal coupling

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
Sat_RadEffArea = 1.36;          % effektive Fl�che des Radiators, Pyramide
Sat_RadName = 'Radiator';
Sat_CellEff = 0.34;             % Effizienz Solarzellen
Sat_CellName = 'SolarCells';
Sat_SurfNum = 21;               % Anzahl der Oberfl�chen wie in Simulation berechnet

% % % Simulation
t_Res = 120;                    % [s] zeitliche Aufl�sung
t_ResAcc = t_Res/6;             % [s] zeitliche Aufl�sung Zugriff, nur wenn f_UseMeanLoads == 0
t_Range = [0 0.5/365] + 0.0;      % simulierte Zeit [start ende], [0 1] voll
t_Step = 1;                     % Schrittweite
t_IntLim = [1/60 5];            % Grenzen der Zeitkonstanten f�r die Simulation der thermischen Kopplung

% % % Pfade
d_Bas = 'Data\';
d_Dat = [d_Bas 'ERNST_i97_a500_r00_t120_stiff\'];
d_Par = 'Data\';
d_Cmp = [d_Par '_Komponenten.txt'];
d_Mat = [d_Par '_Materialien.txt'];
d_Sur = [d_Par '_Struktur.txt'];
d_TCo = [d_Par '_Waermeleitung.txt'];
d_XLo = [d_Par '_XLoads.txt'];
d_TEx = [d_Par 'Temperatur'];
% d_Suff = ' - short';
d_Suff = '';

% % % Optionen zum Zu- und Abschalten bestimmter Effekt
f_Sun = 1;      % Sonneneinstrahlung
f_Cmp = 1;      % Abw�rme Komponenten
f_Alb = 1;      % Albedo
f_EIR = 1;      % Erde IR
f_Emi = 1;      % Emission Oberfl�chen
f_TCo = 1;      % Thermische Kopplung
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
f_DrawOnTop = 0;        % dr�berzeichnen
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
    Sat_Cmp = struct();         % Komponenten, W�rmeentwicklung
    Sat_Struct = struct();      % Oberfl�chen, Verkn�pfung mit Mat und Cmp, Fl�chen und Gewicht
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
d_SubSol = sprintf('%s%s%s SunAngles.csv',d_Dat,satName,fmt_Base);

% Data and formatStrings(NEW and FIXED for Albedo)
d_SubSol_fix = sprintf('%sSunAngles.csv',d_Dat);
fStr_SubSol = '%d %s %d %12s,%f,%f,%f'; % day (d) month (s) year (d) time (12s), azimuthal (f), elevation (f), subsolar (f)
d_EarthAngles = sprintf('%sEarthAngles.csv',d_Dat);
fStr_EarthAngles = '%d %s %d %12s,%f,%f'; % day (d) month (s) year (d) time (12s), azimuthal (f), elevation (f)

% Formatstrings
fstrAreas = ['%d %s %d %12s' repmat(';%f',1,Sat_SurfNum)];
fstrPower = '%d %s %d %12s;%f';
fstrAccess = '%f,%f %s %f %f:%f:%f,%f %s %f %f:%f:%f,%f';
fstrDateIn = '%.0f %s %.0f %s';
fstrDateOut = '%.0f %s %.0f %02.0f:%02.0f:%02.0f';
fstrDate = 'dd mmm yyyy HH:MM:SS';
fstrCmp = '%s %d %f %f';
fstrMat = '%s %f %f %f';
fstrStr = '%s %s %s %f %f %s';
fstrSubSol = '%d %s %d %12s,%f,%f,%f,%f';
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
r = (RE)/(RE+Alt*1e3);
%%%%%% CONSTANT RHO TEST
rho = 140;

% maximaler Subsolarwinkel f�r Albedo
AlbMaxAngle = 100; % [deg] NEW (albedo formula depends on cos(0.9*theta)^1.5, complex numbers above 100 deg
%acosd(RE/au) + acosd(RE/(RE+Alt*1e3)); % OLD

%% Daten einlesen
if (f_ReloadAllData == 1)
    fprintf('Einlesen der Daten:\n');
    
    % % Fl�chen und Leistung
    % Reihenfolge extrahieren
    fid = fopen(d_AreaS);
    dat_temp = textscan(fid, fstrOrder, h_Area-2, 'HeaderLines', 2);
    dat_Order = dat_temp(3);
    fclose(fid);
    
    % Fl�chen extrahieren
    fprintf(' ... Fl�chen aus Sonnensicht ...\n');
    dat_temp = ReadCSV(d_AreaS,fstrAreas,h_Area);
    dat_AreaS = dat_temp(5:end);
    dat_Time = dat_temp(1:4);
    fprintf(' ... Fl�chen aus Erdsicht ...\n');
    dat_temp = ReadCSV(d_AreaE,fstrAreas,h_Area);
    dat_AreaE = dat_temp(5:end);
    if (f_PlotGenPower == 1)
        fprintf(' ... Erzeugte Leistung ...\n');
        dat_temp = ReadCSV(d_Power,fstrPower,h_Power);
        dat_Power = dat_temp(5:end);
    end
    
    % % Subsolarwinkel
    fprintf(' ... Subsolarwinkel ...\n');
    dat_temp = ReadCSV(d_SubSol,fstrSubSol,1); % OLD
    dat_SubSol = dat_temp(end); % OLD
    % FIXED, NEW CODE
    dat_temp = ReadCSV(d_SubSol_fix,fStr_SubSol,1); % NEW
    dat_SubSol_fix = dat_temp(end); % NEW
    dat_temp = ReadCSV(d_EarthAngles,fStr_EarthAngles,1); % NEW
    dat_EarthAngles = dat_temp(end-1:end); % NEW
    
    
    % % Zugriffszeiten
    if (f_IgnoreAccessIntervals == 0)
        fprintf(' ... Zugriffszeiten ...\n');
        dat_temp = ReadCSV(d_Access,fstrAccess,1);
        % �bersetzen in Sekunden seit Start
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
    
    % % Fl�chen: Name, Material, Gr��e, Komponenten
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
        % Position in Fl�chendatei
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
    
    % % Thermische Kopplung der Oberfl�chen
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
        % 3 Mit Werten aus Liste f�llen
        for ii = 1:numel(dat{1})
            % Indices finden f�r beide Oberfl�chen
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
    
    % % extra Lastf�lle
    fprintf(' ... Lastf�lle ...\n');
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

% Temperaturen f�r hot (1) und cold (2) case, jede Fl�che, jeder Zeitschritt
T = zeros(2,numel(Sat_Struct),N+1);
% Starttemperatur
T(:,:,ran(1)) = T_Start;

% Oberfl�chen sortieren
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
    % % % (A) Isolierte Betrachtung
    for ss = f_IncludedParts
        
        % Leistungszuwachs
        dP_C = 0;
        dP_H = 0;
        
        % Index der Komponente in Fl�chendatei
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
        if (f_Sun == 1)
            % Fl�che
            if (cIdx >= 0)
                A = dat_AreaS{1,cIdx}(tt);
            else
                A = 0;
            end
            % geringerer Energieeingang f�r Solarzellen
            if (strcmp(Sat_Struct(ss).name(1:end-1),Sat_CellName) == 1)
                A = A * (1 - Sat_CellEff);
            end
            % cold case
            % aufgenommene Leistung
            P_C = A * Sat_Mat(sIdx).abs * Sol_Flux(1);
            % wenn nach au�en gerichtete Fl�che null, dann auch aufgenommene Leistung 0
            P_C = P_C * (Sat_Struct(ss).size > 0);
            % hot case
            P_H = A * Sat_Mat(sIdx).abs * Sol_Flux(2);
            P_H = P_H * (Sat_Struct(ss).size > 0);
            % Zuwachs
            dP_C = dP_C + P_C;
            dP_H = dP_H + P_H;
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
                                case 1 % zeitaufgel�st, nach Access
                                    dur = 0;
                                    % Sekunden z�hlen, die die Komponente eingeschaltet ist
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
%                                 case 2 % zeitaufgel�st, nach Sentinelwinkel
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
        
        % % % (3) IR Abstrahlung des Satelliten
        if (f_Emi == 1)
            % Fl�che * Emissivit�t
            if (cIdx >= 0)
                A = Sat_Struct(ss).size;
            else
                A = 0;
            end
            % Radiator hat strukturierte Oberfl�che
            if (strcmp(Sat_Struct(ss).name,Sat_RadName) == 1)
                A = A * Sat_RadEffArea;              
            end            
            % Leistungs�nderung
            dP_C = dP_C - ksb * (T(1,ss,tt)^4 - T_Space^4) * A * Sat_Mat(sIdx).emi;
            dP_H = dP_H - ksb * (T(2,ss,tt)^4 - T_Space^4) * A * Sat_Mat(sIdx).emi;
        end
        
        % % % (4) Erde IR
        if (f_EIR == 1)
            % Fl�che aus Sicht der Erde
            if (cIdx >= 0)
                A = dat_AreaE{1,cIdx}(tt);
            else
                A = 0;
            end
            % Korrektur um Inklination
            inc_c = Inc;
            if(Inc > 90)
                inc_c = 180 - Inc;
            end
            facInc_C = EIR_OrbitInc_C(find(EIR_OrbitInc_C(:,1) <= inc_c,1,'last'),2);
            facInc_H = EIR_OrbitInc_H(find(EIR_OrbitInc_H(:,1) <= inc_c,1,'last'),2);
            % Leistungs�nderung (emi gibt den Wert der Emission/Absorption im IR an)
            P_C = A * Sat_Mat(sIdx).emi * facInc_C;
            P_H = A * Sat_Mat(sIdx).emi * facInc_H;
            % wenn nach au�en gerichtete Fl�che null, dann auch aufgenommene Leistung 0
            P_C = P_C * (Sat_Struct(ss).size > 0);
            P_H = P_H * (Sat_Struct(ss).size > 0);
            % Energie�nderung
            dP_C = dP_C + P_C;
            dP_H = dP_H + P_H;
        end
        
        % % % (5) Albedo (nur wenn Subsolarwinkel kleiner als maximaler Winkel)
        if (f_Alb == 1 && dat_SubSol_fix{1}(tt) < AlbMaxAngle)
            % Fl�che aus Sicht der Erde
            if (cIdx >= 0)
                A = dat_AreaE{1,cIdx}(tt);
            else
                A = 0;
            end
            % Korrektur um Inklination f�r SSO
            inc_c = Inc;
            if(Inc > 90)
                inc_c = 180 - Inc;
            end
            if (f_UseFixedAlbedo < 0)
                % Einfluss anhand der Inklination
                facInc_C = Alb_OrbitInc_C(find(Alb_OrbitInc_C(:,1) <= inc_c,1,'last'),2);
                facInc_H = Alb_OrbitInc_H(find(Alb_OrbitInc_H(:,1) <= inc_c,1,'last'),2);
                % Korrektur anhand des Subsolarwinkels
                alpha = dat_SubSol_fix{1}(tt);
                a_c = Alb_KorrSubsolar(find(Alb_KorrSubsolar(:,1) <= alpha,1,'last'),2);
                facInc_C = facInc_C + a_c;
                facInc_H = facInc_H + a_c;
            else
                facInc_C = f_UseFixedAlbedo;
                facInc_H = f_UseFixedAlbedo;
            end
            % Leistungs�nderung
            P_C = A * Sat_Mat(sIdx).abs * albedoSolarFlux(Sol_Flux(1),dat_SubSol_fix{1}(tt),r,rho); % NEW
            P_H = A * Sat_Mat(sIdx).abs * albedoSolarFlux(Sol_Flux(2),dat_SubSol_fix{1}(tt),r,rho); % NEW
            %P_C = A * Sat_Mat(sIdx).abs * Sol_Flux(1) * facInc_C; % OLD
            %P_H = A * Sat_Mat(sIdx).abs * Sol_Flux(2) * facInc_H; % OLD
            % wenn nach au�en gerichtete Fl�che null, dann auch aufgenommene Leistung 0
            P_C = P_C * (Sat_Struct(ss).size > 0);
            P_H = P_H * (Sat_Struct(ss).size > 0);
            % Leistungszuwachs
            dP_C = dP_C + P_C;
            dP_H = dP_H + P_H;
        end
        
        
        % % % (X) Temperaturberechnung
    
        % Energie erhalten durch Skalierung Leistung mit Zeitschritt
        E_C = dP_C * (t_Res * t_Step);
        E_H = dP_H * (t_Res * t_Step);

        % Umrechnen in Temperatur�nderung
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
    
    % % % (B) W�rmeleitung, etc. zwischen Einzelkomponenten
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
    
end
fprintf('\n ... fertich.\n');

    
%% Plotten der Ergebnisse
    
fprintf('Ergebnisse zeichnen ...');
f_DrawParts = [1 4 14 15 22];

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
    ylabel 'temperature [�C]'
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
    
    % l�ngste Zeit ohne Sonne
    zeroPower = relPower == 0;
    % Indices mit 0
    idx0 = find(zeroPower);
    % Differenz nicht-konsekutiver Indices
    sum0 = [0 cumsum(diff(idx0)~=1)];
    % Indices von l�ngster Sequenz
    if (numel(sum0) > 1)
        seq0 = idx0(sum0==mode(sum0));
        % Zeitspanne
        noPower = t_Res * numel(seq0);
        fprintf('\n\tL�ngste Zeit ohne Sonnenlicht: %.1f min,\n\tMissionzeit: %.3f -> %.3f d.\n',noPower/60,t_Vec(seq0(1)),t_Vec(seq0(end)));
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