function qr = qrCalculation(Ts, Ta, D0,epsilon)
% Radiated heat loss [W/m]
qr = 17.8*D0*epsilon*(((Ts+273)/100)^4-((Ta+273)/100)^4);
end