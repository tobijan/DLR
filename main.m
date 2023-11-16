close all

%% Setup
tempdata = readmatrix('temp.csv'); 
winddata = readmatrix('wind.csv');
solardata = readmatrix('solar.csv');

D0 = 27.0002*10^-3; % Outer conductor diameter [m]
% Dc = 1; % Conductor core diameter [m]
Area = 27.0002*10^-3; % Projected area of conductor [m^2/m]

Ts = 80; % Conductor surface temperature [°C]

% 566995:566995+8760
Ta = tempdata(566995:566995+8759,3); % Ambient air temperature [°C]
Tlow = 25; % Low average conductor temperature for which ac resistance is specified [°C]
Thigh = 100; % High average conductor temperature for which ac resistance is specified [°C]
R_Tlow = 0.0736325*10^-3 ; % Lower temperature resistance [ohm/m]
R_Thigh = 0.088359*10^-3; % Higher temperature resistance [ohm/m]

% 341881:341881+8759
irradiation = solardata(341881:341881+8759,3);
Zl = 90; % Azimuth of line [deg]
Lat = 65; % Degrees of latitude [deg]
He = 40; % Elevation of conductor above sea level [m]

% 577612:577612+8759
Vw = winddata(577612:577612+8759,5); % Wind velocity [m/s]
phi = mod(90+winddata(577612:577612+8759,3),360); % Angle between the wind direction and the conductor axis [deg]. Look at page 11 if issues occur.

alpha = 0.8; % Solar absorptivity (.23 to .91) []
epsilon = 0.8; % Emissivity (.23 to .91) []

I = zeros(length(Ta),1);
qc = zeros(length(Ta),1);
qr = zeros(length(Ta),1);
qs = zeros(length(Ta),1);
for i=1:length(Ta)
    N = floor((i - 1) / 24) + 1; % Day of the year
    hour = mod(i - 1, 24); % Hour of the day
    omega = 15 * (hour - 12); % Hour angle relative to noon [deg]. 15*(Time-12), at 11AM, Time = 11 and the Hour angle= –15 deg
    
    % Natural convection heat loss
    Tfilm = (Ts+Ta(i))/2; % Average temperature of boundary layer [°C]
    rhof = (1.293-1.525*10^-4*He+6.379*10^-9*He^2)/(1+0.00367*Tfilm); % Air density [kg/m3]
    qcn(i) = 3.645*rhof^0.5*D0^0.75*(Ts-Ta(i))^1.25; % [W/m]
    
    % Forced convection
    if Vw(i) == 0 % No wind
        qcf(i) = 0;
    else
        muf = (1.458*10^-6*(Tfilm+273)^1.5)/(Tfilm+383.4); % Absolute (dynamic) viscosity of air [kg/m-s]
        kf = 2.424*10^-2+7.477*10^-5*Tfilm-4.407*10^-9*Tfilm^2; % Thermal conductivity of air at temperature Tfilm [W/(m-°C)
        Kangle(i) = 1.194-cosd(phi(i))+0.194*cosd(2*phi(i))+0.368*sind(2*phi(i));
        NRe(i) = (D0*rhof*Vw(i))/muf; % Reynolds number []. Transition between qc1 and qc2 at 1000, or compare and take highest?
        if NRe(i)<1000
            qcf(i) = Kangle(i)*(1.01+1.35*NRe(i)^0.52)*kf*(Ts-Ta(i)); % [W/m]
        else
            qcf(i) = Kangle(i)*0.754*NRe(i)^0.6*kf*(Ts-Ta(i)); % [W/m]
        end
    end
    
    qcoptions = [qcn(i) qcf(i)];

    qc(i) = max(qcoptions);
    
    % Radiated heat loss
    qr(i) = 17.8*D0*epsilon*(((Ts+273)/100)^4-((Ta(i)+273)/100)^4); % [W/m]
    
    % Solar heat gain
    delta = 23.46*sind((284+N)/365*360); % Solar declination (–23.45 to +23.45) [deg]
    Hc = asind(cosd(Lat)*cosd(delta)*cosd(omega)+sind(Lat)*sind(delta)); % Altitude of sun [deg]
    if Hc > 1 % Sun is above the horizon
        chi = sind(omega)/(sind(Lat)*cosd(omega)-cosd(Lat)*tand(delta)); % Solar azimuth variable []
        if -180 <= omega && omega < 0
            if chi >= 0
                C=0; % Solar azimuth constant [deg]
            else
                C=180;
            end
        else
            if chi >= 0
                C=180;
            else
                C=360;
            end
        end
        
        Zc = C + atand(chi); % Azimuth of sun [deg]
        
        A = -42.2391;
        B = 63.8044;
        C = -1.9220;
        D = 3.46921e-2;
        E = -3.61118e-4;
        F = 1.94318e-6;
        G = -4.07608e-9;
        Qs = A + B*Hc + C*Hc^2 + D*Hc^3 + E*Hc^4 + F*Hc^5 + G*Hc^6;
        
        a = 1;
        b = 1.148*10^-4;
        c = -1.108*10^-8;
        Ksolar = a + b*He + c*He^2; % Solar altitude correction factor []
        Qse = Ksolar*Qs; 
        
        theta = acosd(cosd(Hc)*cosd(Zc-Zl)); % Effective angle of incidence of the sun rays [deg]
        qs(i) = alpha*Qse*sind(theta)*Area; % Solar heat gain [W/m]
        %qs(i) = irradiation(i)*Area*alpha;
    else
        % Sun is down, set solar heat gain to zero
        qs(i) = 0;
        %qs(i) = irradiation(i)*Area*alpha;
    end
    
    % Conductor resistance
    R_Ts = ((R_Thigh-R_Tlow)/(Thigh-Tlow))*(Ts-Tlow)+R_Tlow; % AC resistance of conductor at temperature Tavg [ohm/m]
    
    % Ampacity
    I(i) = sqrt((qc(i)+qr(i)-qs(i))/R_Ts); % [A]
end

U = 400*10^3;
PF = 0.95;

power = sqrt(3)*U*I*PF*10^-6;

figure
plot(1:length(Ta),power)
ylim([0 3000])
title('Transmission capacity during 2022')
xlabel('Time [h]')
ylabel('Capacity [MW]')

figure
plot(1:length(Ta),Ta)
title('Ambient temperature during 2022')
xlabel('Time [h]')
ylabel('Temperature [C]')


figure
plot(1:length(Ta),Vw)
title('Wind speed during 2022')
xlabel('Time [h]')
ylabel('Wind speed [m/s]')

figure
plot(1:length(Ta),phi)
title('Wind direction during 2022')
xlabel('Time [h]')
ylabel('Wind direction [deg]')

figure
plot(1:length(Ta),irradiation)
title('Solar irradiation during 2022')
xlabel('Time [h]')
ylabel('Irradiation [W/m^2]')

figure
plot(1:length(Ta),qc)
title('Convection heat loss during 2022')
xlabel('Time [h]')
ylabel('qc [W/m]')

figure
plot(1:length(Ta),qr)
title('Radiative heat loss during 2022')
xlabel('Time [h]')
ylabel('qr [W/m]')

figure
plot(1:length(Ta),qs)
title('Solar heat gain during 2022')
xlabel('Time [h]')
ylabel('qs [W/m]')


