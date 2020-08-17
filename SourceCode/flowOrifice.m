function [q_m, q_V, epsilon, Re_D, A, C, error, i]=flowOrifice(dP, rho, Kappa, p1, mu, errorInDischargeCoeff,sensitivity)
%************************Dagfinn 2018*******************************
%ISO 5167 - 2
%[q_m, q_V, epsilon, Re_D, A, C, error, i]=flowOrifice(350, 16.563, 1.267, 20.5, 0.00001213);
%[q_m, q_V, epsilon, Re_D, A, C, error, i]=flowOrifice(15474.2, 23.67187, 1.27456, 30, 0.000010414);
%[q_m] = kg/s
%[q_V] = m^3/s
%[dP] = mBar, Pressure differaence on across orifice
%[rho] = kg/m3
%[Kappa] = -  Isentropic Exponent
% [P1] = Pa , Pressure upstream orifice
% [P2] = Pa , Pressure downstream orifice
% [mu] = Pa*s, Dynamic viscosity oppstrøms orifice
%*******************************************************************

p1 = (p1+1.01325)*100000; %Convert from barg to Pa
dP = dP*100; %Convert from mBar to Pa

d = 0.165924; %Inner diameter flow orifice
D = 0.273418; %Pipe Diameter
beta=d/D;
p2 = p1-dP; %Pressure downstream flow orifice

C = 0.6;%initial value of C


epsilon = 1-(0.351+0.256*beta^4+0.93*beta^8)*(1-(p2/p1)^(1/Kappa));
q_m = (C*epsilon*pi*d^2*sqrt(2*dP*rho))/(sqrt(1-beta^4)*4);

error=3;%initial value of error

i=0;
while error>0.00000000000001
    i=i+1;
    Re_D = (4*q_m)/(pi*mu*D);
    A = (19000*beta/Re_D)^0.8;
    C=0.5961+0.0261*beta^2-0.216*beta^8+0.000521*(10^6*beta/Re_D)^0.7+(0.0188+0.0063*A)*beta^3.5*(10^6/Re_D)^0.3;
    q_m_Old=q_m;
    q_m = (C*epsilon*pi*d^2*sqrt(2*dP*rho))/(sqrt(1-beta^4)*4);  
    error=abs(q_m_Old-q_m)/q_m_Old;
end

if errorInDischargeCoeff %use uncertainty in the discharge coeff
    if ~sensitivity
        sigmaC = (1.667*beta - 0.5)*C/100;
        C = normrnd(C,sigmaC);
        q_m = (C*epsilon*pi*d^2*sqrt(2*dP*rho))/(sqrt(1-beta^4)*4);
    end
end

q_V=q_m/rho;
