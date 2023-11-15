close all

%% Setup
tempdata = readmatrix('temp.csv'); 
winddata = readmatrix('wind.csv');

D0 = 27.0002*10^-3; % Outer conductor diameter [m]
% Dc = 1; % Conductor core diameter [m]
Area = 27.0002*10^-3; % Projected area of conductor [m^2/m]

Ts = 60; % Conductor surface temperature [°C]


Ta = tempdata(566995:566995+8760,3); % Ambient air temperature [°C]
Tlow = 25; % Low average conductor temperature for which ac resistance is specified [°C]
Thigh = 100; % High average conductor temperature for which ac resistance is specified [°C]
R_Tlow = 0.0736325*10^-3 ; % Lower temperature resistance [ohm/m]
R_Thigh = 0.088359*10^-3; % Higher temperature resistance [ohm/m]

Zl = 90; % Azimuth of line [deg]
Lat = 65; % Degrees of latitude [deg]
He = 40; % Elevation of conductor above sea level [m]

Vw = winddata(566995:566995+8760,5); % Wind velocity [m/s]
phi = winddata(566995:566995+8760,3); % Angle between the wind direction and the conductor axis [deg]. Look at page 11 if issues occur.

alpha = 0.5; % Solar absorptivity (.23 to .91) []
epsilon = 0.5; % Emissivity (.23 to .91) []

I = zeros(length(Ta),1);
for i=1:length(Ta)
    N = floor((i - 1) / 24) + 1; % Day of the year
    hour = mod(i - 1, 24); % Hour of the year
    omega = 15 * (hour - 12); % Hour angle relative to noon [deg]. 15*(Time-12), at 11AM, Time = 11 and the Hour angle= –15 deg
    
    % Natural convection heat loss
    Tfilm = (Ts+Ta(i))/2; % Average temperature of boundary layer [°C]
    rhof = (1.293-1.525*10^-4*He+6.379*10^-9*He^2)/(1+0.00367*Tfilm); % Air density [kg/m3]
    qcn = 3.645*rhof^0.5*D0^0.75*(Ts-Ta(i))^1.25; % [W/m]
    
    % Forced convection
    muf = (1.458*10^-6*(Tfilm+273)^1.5)/(Tfilm+383.4); % Absolute (dynamic) viscosity of air [kg/m-s]
    kf = 2.424*10^-2+7.477*10^-5*Tfilm-4.407*10^-9*Tfilm^2; % Thermal conductivity of air at temperature Tfilm [W/(m-°C)
    Kangle = 1.194-cos(phi(i))+0.194*cos(2*phi(i))+0.368*sin(2*phi(i));
    NRe = (D0*rhof*Vw(i))/muf; % Reynolds number []. Transition between qc1 and qc2 at 1000, or compare and take highest?
    if NRe<1000
        qcf = Kangle*(1.01+1.35*NRe^0.52)*kf*(Ts-Ta(i)); % [W/m]
    else
        qcf = Kangle*0.754*NRe^0.6*kf*(Ts-Ta(i)); % [W/m]
    end
    
    qcoptions = [qcn qcf];

    qc = max(qcoptions);
    
    % Radiated heat loss
    qr = 17.8*D0*epsilon*(((Ts+273)/100)^4-((Ta(i)+273)/100)^4); % [W/m]
    
    % Solar heat gain
    delta = 23.46*sin((284+N)/365*360); % Solar declination (–23.45 to +23.45) [deg]
    Hc = asin(cos(Lat)*cos(delta)*cos(omega)+sin(Lat)*sin(delta)); % Altitude of sun [deg]
    chi = sin(omega)/(sin(Lat)*cos(omega)-cos(Lat)*tan(delta)); % Solar azimuth variable []
    
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
    
    Zc = C + atan(chi); % Azimuth of sun [deg]
    
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
    
    theta = acos(cos(Hc)*cos(Zc-Zl)); % Effective angle of incidence of the sun rays [deg]
    qs = alpha*Qse*sin(theta)*Area; % Solar heat gain [W/m]
    
    % Conductor resistance
    R_Ts = ((R_Thigh-R_Tlow)/(Thigh-Tlow))*(Ts-Tlow)+R_Tlow; % AC resistance of conductor at temperature Tavg [ohm/m]
    
    % Ampacity
    I(i) = sqrt((qc+qr-qs)/R_Ts); % [A]
end

U = 400*10^3;
PF = 0.95;

power = sqrt(3)*U*I*PF;

