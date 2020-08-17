function [Compressor] = CompressorPerformanceEquationNode(FP,IS,Constants)
%% Calculate compressor performance parameters

%polytropic exponent
Compressor.PVTsim.n = log((FP.Out.PVTsim.Pressure+Constants.atm)./(FP.In.PVTsim.Pressure+Constants.atm))./log((FP.Out.PVTsim.Result.mixDensity)./(FP.In.PVTsim.Result.mixDensity));

% Polytropic head
Compressor.PVTsim.PolyHead = 100*(Compressor.PVTsim.n./(Compressor.PVTsim.n-1)).*((FP.Out.PVTsim.Pressure+Constants.atm)./FP.Out.PVTsim.Result.mixDensity-(FP.In.PVTsim.Pressure+Constants.atm)./FP.In.PVTsim.Result.mixDensity);

% Actual head thermodynamic
Compressor.PVTsim.ActHead_thermo = FP.Out.PVTsim.Result.mixMassEnthalpy-FP.In.PVTsim.Result.mixMassEnthalpy;

% Actual head mechanical
Compressor.PVTsim.ActHead_Mech = (IS.PowerTorque./FP.In.PVTsim.Result.totalMassFlow);

% Power themodynamic
Compressor.PVTsim.Power_thermo = Compressor.PVTsim.ActHead_thermo.*FP.In.PVTsim.Result.totalMassFlow;

% Power mechanical
Compressor.PVTsim.Power_Mech = IS.PowerTorque;

% Efficiency themodynamic
Compressor.PVTsim.PolyEff_thermo = Compressor.PVTsim.PolyHead./Compressor.PVTsim.ActHead_thermo;

% Efficiency mechanical
Compressor.PVTsim.PolyEff_Mech = Compressor.PVTsim.PolyHead./(IS.PowerTorque./FP.In.PVTsim.Result.totalMassFlow);

% Impeller Tip Speed
Compressor.PVTsim.ImpTipSpeed = (pi*IS.ImpellerDiam.*IS.Speed)./60;

% compressor actual flow
Compressor.PVTsim.ActualFlow = 3600*FP.In.PVTsim.Result.totalMassFlow./FP.In.PVTsim.Result.mixDensity;

% phi dry
Compressor.PVTsim.FlowCoeffDry = 4*(FP.In.PVTsim.Result.gasMassFlow./FP.In.PVTsim.Result.gasDensity)./...
    (pi*IS.ImpellerDiam.^2.*Compressor.PVTsim.ImpTipSpeed);

% phi wet
Compressor.PVTsim.FlowCoeffWet = 4*(FP.In.PVTsim.Result.totalMassFlow./FP.In.PVTsim.Result.mixDensity)./...
    (pi*IS.ImpellerDiam.^2.*Compressor.PVTsim.ImpTipSpeed);

% Ma dry (mach number dry)
u = Compressor.PVTsim.ImpTipSpeed;
a_gas = FP.In.PVTsim.Result.gasSoundSpeed;
Compressor.PVTsim.MachNumberDry = u./a_gas;

% Ma wet (mach number wet based woods Hundseid multphase speed of sound) 
u = Compressor.PVTsim.ImpTipSpeed;
alpha_gas = FP.In.PVTsim.Result.gasVolumeFraction;
alpha_oil = FP.In.PVTsim.Result.oilVolumeFraction;
alpha_water = FP.In.PVTsim.Result.waterVolumeFraction;
rhoLiq = (FP.In.PVTsim.Result.oilVolumeFraction.*FP.In.PVTsim.Result.oilDensity + ...
    FP.In.PVTsim.Result.waterVolumeFraction.*FP.In.PVTsim.Result.waterDensity)./...
    (FP.In.PVTsim.Result.oilVolumeFraction + FP.In.PVTsim.Result.waterVolumeFraction);
rhoGas = FP.In.PVTsim.Result.gasDensity;
OneOverdelta = rhoLiq./rhoGas;
konst = (1 - alpha_gas)./alpha_gas;
a_woods = (1+konst)./(alpha_gas.*(1+(konst.*OneOverdelta)));
Compressor.PVTsim.MachNumberWet = u./a_woods;

% Machine/tip Reynolds number Dry
b = 0.016; %Impeller Tip width
visc_gas = FP.In.PVTsim.Result.gasViscosity;
Compressor.PVTsim.Re_M_Dry = (rhoGas.*u.*b)./visc_gas;

% Machine/tip Reynolds number Wet
rho_oil = FP.In.PVTsim.Result.oilDensity;
rho_water = FP.In.PVTsim.Result.waterDensity;
beta_gas = FP.In.PVTsim.Result.gasWtFraction;
beta_oil = FP.In.PVTsim.Result.oilWtFraction;
beta_water = FP.In.PVTsim.Result.waterWtFraction;
rho_m = alpha_gas.*rhoGas + alpha_oil.*rho_oil + alpha_water.*rho_water;
visc_oil = FP.In.PVTsim.Result.oilViscosity;
visc_water = FP.In.PVTsim.Result.waterViscosity;
visc_m_MacAdams = (beta_gas./visc_gas + beta_oil./visc_oil + beta_water./visc_water).^-1; % source Effective property models for homogeneous two-phase flows
visc_m_Bingham = (alpha_gas./visc_gas + alpha_oil./visc_oil + alpha_water./visc_water).^-1; %source TwoPhase Flow
visc_m_Cicchitti = beta_gas.*visc_gas + beta_oil.*visc_oil + beta_water.*visc_water; 
visc_m_My = alpha_gas.*visc_gas + alpha_oil.*visc_oil + alpha_water.*visc_water;
Compressor.PVTsim.Re_M_wet_MacAdams = (rho_m.*u.*b)./visc_m_MacAdams;
Compressor.PVTsim.Re_M_wet_Bingham = (rho_m.*u.*b)./visc_m_Bingham;
Compressor.PVTsim.Re_M_wet_Cicchitti = (rho_m.*u.*b)./visc_m_Cicchitti;
Compressor.PVTsim.Re_M_wet_My = (rho_m.*u.*b)./visc_m_My;

%Specific volume ratio machine
Compressor.PVTsim.SpecVolRatio = FP.Out.PVTsim.Result.mixDensity./FP.In.PVTsim.Result.mixDensity;

%Density Ratio
Compressor.PVTsim.Dr = rhoLiq./rhoGas;

%LMF
Compressor.PVTsim.LMF = 1-FP.In.PVTsim.Result.gasWtFraction;

%GVF
Compressor.PVTsim.GVF = FP.In.PVTsim.Result.gasVolumeFraction;

end

