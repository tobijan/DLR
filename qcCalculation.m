function qc = qcCalculation(Ts, Ta, D0, He, Vw, phi)
% Average temperature of boundary layer [°C]
Tfilm = (Ts+Ta)/2; 

% Air density [kg/m3]
rhof = (1.293-1.525*10^-4*He+6.379*10^-9*He^2)/(1+0.00367*Tfilm); 

% Natural convection heat loss [W/m]
qcn = 3.645*rhof^0.5*D0^0.75*(Ts-Ta)^1.25; 

if Vw == 0 % No wind
    qcf = 0;
else
    % Absolute (dynamic) viscosity of air [kg/m-s]
    muf = (1.458*10^-6*(Tfilm+273)^1.5)/(Tfilm+383.4);
    
    % Thermal conductivity of air at temperature Tfilm [W/(m-°C)
    kf = 2.424*10^-2+7.477*10^-5*Tfilm-4.407*10^-9*Tfilm^2; 
    
    Kangle = 1.194-cosd(phi)+0.194*cosd(2*phi)+0.368*sind(2*phi);
    
    % Reynolds number []. Transition qcf equation between qc1 and qc2 at 1000
    NRe = (D0*rhof*Vw)/muf; 

    % Forced convection
    if NRe<1000
        qcf = Kangle*(1.01+1.35*NRe^0.52)*kf*(Ts-Ta); % [W/m]
    else
        qcf = Kangle*0.754*NRe^0.6*kf*(Ts-Ta); % [W/m]
    end
end
% Choose the highest of the natural and forced convection
qcoptions = [qcn qcf];
qc = max(qcoptions);
end
