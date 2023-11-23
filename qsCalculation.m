function qs = qsCalculation(N,Lat,omega,Zl,Area,He,alpha)
% Solar declination (–23.45 to +23.45) [deg]
delta = 23.46*sind((284+N)/365*360); 

% Altitude of sun [deg]
Hc = asind(cosd(Lat)*cosd(delta)*cosd(omega)+sind(Lat)*sind(delta)); 

if Hc > 1 % Sun is above the horizon
    % Solar azimuth variable []
    chi = sind(omega)/(sind(Lat)*cosd(omega)-cosd(Lat)*tand(delta)); 

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
    % Azimuth of sun [deg]
    Zc = C + atand(chi);
    
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
    % Solar altitude correction factor []
    Ksolar = a + b*He + c*He^2; 
    Qse = Ksolar*Qs; 
    
    theta = acosd(cosd(Hc)*cosd(Zc-Zl)); % Effective angle of incidence of the sun rays [deg]
    qs = alpha*Qse*sind(theta)*Area; % Solar heat gain [W/m]
    
else
    % Sun is down, set solar heat gain to zero
    qs = 0;
end
end