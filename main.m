close all

%% Setup
% Forecasted temperature
temp_for = readmatrix('temp.csv'); 

% Actual temperature
temp_act = readmatrix('temp.csv');

% Forecasted wind
wind_for = readmatrix('wind.csv');

% Actual wind
wind_act = readmatrix('wind.csv');

% Ambient air temperature [°C]
Ta_for = temp_for(566995:566995+8759,3); 
Ta_act = temp_act(566995:566995+8759,3); 

% Wind velocity [m/s]
Vw_for = wind_for(577612:577612+8759,5); 
Vw_act = wind_act(577612:577612+8759,5);

% Wind direction. Angle between the wind direction and the conductor axis [deg]. Look at page 11 if issues occur.
phi_for = mod(90+wind_for(577612:577612+8759,3),360); 
phi_act = mod(90+wind_act(577612:577612+8759,3),360); 

% Outer conductor diameter [m]
D0 = 27.0002*10^-3;

% Projected area of conductor [m^2/m]
Area = 27.0002*10^-3; 

% Conductor surface temperature [°C]
Ts = 60; 

% Low average conductor temperature for which ac resistance is specified [°C]
Tlow = 25; 

% High average conductor temperature for which ac resistance is specified [°C]
Thigh = 100; 

% Lower temperature resistance [ohm/m]
R_Tlow = 0.0736325*10^-3 ; 

% Higher temperature resistance [ohm/m]
R_Thigh = 0.088359*10^-3; 

% AC resistance of conductor at temperature Ts [ohm/m]
R_Ts = ((R_Thigh-R_Tlow)/(Thigh-Tlow))*(Ts-Tlow)+R_Tlow; 
    
% Azimuth of line [deg]
Zl = 90; 

% Degrees of latitude [deg]
Lat = 65; 

% Elevation of conductor above sea level [m]
He = 40; 

% Solar absorptivity (.23 to .91) []
alpha = 0.8; 

% Emissivity (.23 to .91) []
epsilon = 0.8;

% Voltage [kV]
U = 400*10^3;

% Power factor []
PF = 0.95;

%% Calculate forecasted capacity

I_for = zeros(length(Ta_for),1);
qc_for = zeros(length(Ta_for),1);
qr_for = zeros(length(Ta_for),1);
qs_for = zeros(length(Ta_for),1);

for i=1:length(Ta_for)
    % Convection heat loss
    qc_for(i) = qcCalculation(Ts,Ta_for(i),D0,He,Vw_for(i),phi_for(i));
       
    % Radiated heat loss
    qr_for(i) = qrCalculation(Ts,Ta_for(i),D0,epsilon);

    % Day of the year
    N = floor((i - 1) / 24) + 1; 
    
    % Hour of the day
    hour = mod(i - 1, 24); 

    % Hour angle relative to noon [deg]. 15*(Time-12), at 11AM, Time = 11 and the Hour angle= –15 deg
    omega = 15 * (hour - 12); 

    % Solar heat gain
    qs_for(i) = qsCalculation(N,Lat,omega,Zl,Area,He,alpha);
    
    % Ampacity
    I_for(i) = sqrt((qc_for(i)+qr_for(i)-qs_for(i))/R_Ts); % [A]
end

power_for = sqrt(3)*U*I_for*PF*10^-6;

%% Calculate actual capacity

I_act = zeros(length(Ta_act),1);
qc_act = zeros(length(Ta_act),1);
qr_act = zeros(length(Ta_act),1);
qs_act = zeros(length(Ta_act),1);

for i=1:length(Ta_act)
    % Convection heat loss
    qc_act(i) = qcCalculation(Ts,Ta_act(i),D0,He,Vw_act(i),phi_act(i));
       
    % Radiated heat loss
    qr_act(i) = qrCalculation(Ts,Ta_act(i),D0,epsilon);

    % Day of the year
    N = floor((i - 1) / 24) + 1; 
    
    % Hour of the day
    hour = mod(i - 1, 24); 

    % Hour angle relative to noon [deg]. 15*(Time-12), at 11AM, Time = 11 and the Hour angle= –15 deg
    omega = 15 * (hour - 12); 

    % Solar heat gain
    qs_act(i) = qsCalculation(N,Lat,omega,Zl,Area,He,alpha);
    
    % Ampacity
    I_act(i) = sqrt((qc_act(i)+qr_act(i)-qs_act(i))/R_Ts); % [A]
end

power_act = sqrt(3)*U*I_act*PF*10^-6;

%% Calculate difference
power_diff = power_for-power_act;

%% Plots

% Transmission capacity
figure
plot(1:length(Ta_for),power_for)
ylim([0 3000])
title('Transmission capacity during 2022')
xlabel('Time [h]')
ylabel('Capacity [MW]')

% Ambient temperature
figure
plot(1:length(Ta_for),Ta_for)
title('Ambient temperature during 2022')
xlabel('Time [h]')
ylabel('Temperature [C]')

% Wind speed
figure
plot(1:length(Ta_for),Vw_for)
title('Wind speed during 2022')
xlabel('Time [h]')
ylabel('Wind speed [m/s]')

% Wind direction
figure
plot(1:length(Ta_for),phi_for)
title('Wind direction during 2022')
xlabel('Time [h]')
ylabel('Wind direction [deg]')

% Convection heat loss
figure
plot(1:length(Ta_for),qc_for)
title('Convection heat loss during 2022')
xlabel('Time [h]')
ylabel('qc [W/m]')

% Radiative heat loss
figure
plot(1:length(Ta_for),qr_for)
title('Radiative heat loss during 2022')
xlabel('Time [h]')
ylabel('qr [W/m]')

% Solar heat gain
figure
plot(1:length(Ta_for),qs_for)
title('Solar heat gain during 2022')
xlabel('Time [h]')
ylabel('qs [W/m]')