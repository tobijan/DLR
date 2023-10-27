close all

%% Setup
D0 = 1; % Outer conductor diameter [m]
Dc = 1; % Conductor core diameter [m]
Area = 1; % Projected area of conductor [m^2/m]

Ts = 1; % Conductor surface temperature [°C]
Ta = 1; % Ambient air temperature [°C]
Tlow = 1; % Low average conductor temperature for which ac resistance is specified [°C]
Thigh = 2; % High average conductor temperature for which ac resistance is specified [°C]
R_Tlow = 1; % Lower temperature resistance [ohm/m]
R_Thigh = 2; % Higher temperature resistance [ohm/m]
tao = 1; % Thermal time constant of the conductor [min]. Usually between 5 and 20

Zl = 1; % Azimuth of line [deg]
Lat = 1; % Degrees of latitude [deg]
He = 1; % Elevation of conductor above sea level [m]

Vw = 1; % Wind velocity [m/s]
phi = 1; % Angle between the wind direction and the conductor axis [deg]. Look at page 11 if issues occur.

alpha = 1; % Solar absorptivity (.23 to .91) []
epsilon = 1; % Emissivity (.23 to .91) []

N = 1; % Day of the year
omega = 1; % Hour angle relative to noon [deg]. 15*(Time-12), at 11AM, Time = 11 and the Hour angle= –15 deg



%% Natural convection heat loss
Tfilm = (Ts+Ta)/2; % Average temperature of boundry layer [°C]
rhof = (1.293-1.525*10^-4*He+6.379*10^-9*He^2)/(1+0.00367*Tfilm); % Air density [kg/m3]
qcn = 3.645*rhof^0.5*D0^0.75*(Ts-Ta)^1.25; % [W/m]

%% Forced convection
muf = (1.458*10^-6*(Tfilm+273)^1.5)/(Tfilm+383.4); % Absolute (dynamic) viscosity of air [kg/m-s]
kf = 2.424*10^-2+7.477*10^-5*Tfilm-4.407*10^-9*Tfilm^2; % Thermal conductivity of air at temperature Tfilm [W/(m-°C)
Kangle = 1.194-cos(phi)+0.194*cos(2*phi)+0.368*sin(2*phi);
NRe = (D0*rhof*Vw)/muf; % Reynolds number []. Transition between qc1 and qc2 at 1000, or compare and take highest?
qc1 = Kangle*(1.01+1.35*NRe^0.52)*kf*(Ts-Ta); % [W/m]
qc2 = Kangle*0.754*NRe^0.6*kf*(Ts-Ta); % [W/m]

qc = 1; % Choice depending on largest one 

%% Radiated heat loss
qr = 17.8*D0*epsilon*(((Ts+273)/100)^4-((Ta+273)/100)^4); % [W/m]

%% Solar heat gain
delta = 23.46*sin((284+N)/365*360); % Solar declination (–23.45 to +23.45) [deg]
Hc = asin(cos(Lat)*cos(delta)*cos(omega)+sin(Lat)*sin(delta)); % Altitude of sun [deg]
xeta = sin(omega)/(sin(Lat)*cos(omega)-cos(Lat)*tan(delta)); % Solar azimuth variable []

if -180 <= omega && omega < 0
    if xeta >= 0
        C=0; % Solar azimuth constant [deg]
    else
        C=180;
    end
else
    if xeta >= 0
        C=180;
    else
        C=360;
    end
end

Zc = C + atan(xeta); % Azimuth of sun [deg]

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

theta = acos(cos(Hc)*cos(Zc-Zl)); % Effective angle of incidence of the sun s rays [deg]
qs = alpha*Qse*sin(theta)*Area; % Solar heat gain [W/m]

%% Conductor resistance
Tavg = 1; % Average temperature of aluminum strand layers [°C]
R_Tavg = ((R_Thigh-R_Tlow)/(Thigh-Tlow))*(Tavg-Tlow)+R_Tlow; % AC resistance of conductor at temperature Tavg [ohm/m]

%% Ampacity
I = sqrt((qc+qr-qs)/R_Tavg); % [A]

