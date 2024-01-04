close all

%% Setup
% Forecast data
data_for = readmatrix('24hr_forecast_November1st_January5th_Boden.csv');

% Actual data
data_act = readmatrix('historical_November1st_January5th_Boden.csv');

% Forecast temperature
temp_for = data_for(:,2); 

% Actual temperature
temp_act = data_act(:,3);

% Forecast wind speed
windspeed_for = data_for(:,3);

% Actual wind speed
windspeed_act = data_act(:,5);

% Forecast wind direction
winddir_for = data_for(:,4);

% Actual wind direction
winddir_act = data_act(:,4);

% Ambient air temperature [°C]
Ta_for = temp_for; 
Ta_act = temp_act; 

% Wind velocity [m/s]
Vw_for = round(windspeed_for); 
Vw_act = round(windspeed_act);

% Wind direction. Angle between the wind direction and the conductor axis [deg]. Look at page 11 if issues occur.
phi_for = round(mod(90+winddir_for,360),-1); 
phi_act = round(mod(90+winddir_act,360),-1); 

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
He = 1000; 

% Solar absorptivity (.23 to .91) []
alpha = 0.8; 

% Emissivity (.23 to .91) []
epsilon = 0.8;

% Voltage [kV]
U = 400*10^3;

% Power factor []
PF = 0.95;

%% Calculate Forecast capacity

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
Ta_diff = Ta_for-Ta_act;
speed_diff = Vw_for-Vw_act;
dir_diff = mod(phi_for-phi_act,360);
qc_diff = qc_for-qc_act;
qr_diff = qr_for-qr_act;
qs_diff = qs_for-qs_act;



%% Plots

% Transmission capacity
figure

subplot(1,2,1)
plot(1:length(Ta_for),power_for)
title('Forecast transmission capacity')
xlabel('Time [h]')
ylabel('Capacity [MW]')
ylim([0 2200]); 

subplot(1,2,2)
plot(1:length(Ta_act),power_act)
title('Actual transmission capacity')
xlabel('Time [h]')
ylabel('Capacity [MW]')
ylim([0 2200]);

% Difference in transmission capacity
figure
plot(1:length(Ta_for),power_diff)
title('Difference in transmission capacity')
xlabel('Time [h]')
ylabel('Capacity [MW]')

% Ambient temperature
figure

subplot(1,2,1)
plot(1:length(Ta_for),Ta_for)
title('Forecast ambient temperature')
xlabel('Time [h]')
ylabel('Temperature [C]')

subplot(1,2,2)
plot(1:length(Ta_act),Ta_act)
title('Actual ambient temperature')
xlabel('Time [h]')
ylabel('Temperature [C]')

% Ambient temperature difference
figure

plot(1:length(Ta_for),Ta_diff)
title('Ambient temperature difference')
xlabel('Time [h]')
ylabel('Temperature [C]')

% Wind speed
figure
subplot(1,2,1)

plot(1:length(Ta_for),Vw_for)
title('Forecast wind speed')
xlabel('Time [h]')
ylabel('Wind speed [m/s]')
ylim([0 14])

subplot(1,2,2)
plot(1:length(Ta_act),Vw_act)
title('Actual wind speed')
xlabel('Time [h]')
ylabel('Wind speed [m/s]')
ylim([0 14])

% Wind speed difference
figure

plot(1:length(Ta_for),speed_diff)
title('Wind speed difference')
xlabel('Time [h]')
ylabel('Wind speed [m/s]')

% Wind direction
figure

subplot(1,2,1)
scatter(1:length(Ta_for),phi_for,'filled')
title('Forecast wind direction')
xlabel('Time [h]')
ylabel('Wind direction [deg]')

subplot(1,2,2)
scatter(1:length(Ta_act),phi_act,'filled')
title('Actual wind direction')
xlabel('Time [h]')
ylabel('Wind direction [deg]')

% Wind direction difference
figure

scatter(1:length(Ta_for),dir_diff,'filled')
title('Wind direction difference')
xlabel('Time [h]')
ylabel('Wind direction [deg]')

% Convection heat loss
figure

subplot(1,2,1)
plot(1:length(Ta_for),qc_for)
title('Forecast convection heat loss')
xlabel('Time [h]')
ylabel('qc [W/m]')
ylim([0 800])

subplot(1,2,2)
plot(1:length(Ta_act),qc_act)
title('Actual convection heat loss')
xlabel('Time [h]')
ylabel('qc [W/m]')
ylim([0 800])

% Convection heat loss difference
figure

plot(1:length(Ta_for),qc_diff)
title('Convection heat loss difference')
xlabel('Time [h]')
ylabel('qc [W/m]')

% Radiative heat loss
figure

subplot(1,2,1)
plot(1:length(Ta_for),qr_for)
title('Forecast radiative heat loss')
xlabel('Time [h]')
ylabel('qr [W/m]')

subplot(1,2,2)
plot(1:length(Ta_act),qr_act)
title('Actual radiative heat loss')
xlabel('Time [h]')
ylabel('qr [W/m]')

% Radiative heat loss difference
figure

plot(1:length(Ta_for),qr_diff)
title('Radiative heat loss difference')
xlabel('Time [h]')
ylabel('qr [W/m]')

% Solar heat gain
figure

subplot(1,2,1)
plot(1:length(Ta_for),qs_for)
title('Forecast solar heat gain')
xlabel('Time [h]')
ylabel('qs [W/m]')

subplot(1,2,2)
plot(1:length(Ta_act),qs_act)
title('Actual solar heat gain')
xlabel('Time [h]')
ylabel('qs [W/m]')

% Solar heat gain difference
figure

plot(1:length(Ta_for),qs_diff)
title('Solar heat gain difference')
xlabel('Time [h]')
ylabel('qs [W/m]')