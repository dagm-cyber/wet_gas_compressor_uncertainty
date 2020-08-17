function [Compressor] = CompressorPerformanceStatistics(Compressor)
%% Calculater compressor performance statistics

%polytropic exponent
Compressor.PVTsimSTD.n = std(Compressor.PVTsim.n);
Compressor.PVTsimSTD_Rel.n = Compressor.PVTsimSTD.n./mean(Compressor.PVTsim.n);
Compressor.PVTsimSTD_Perc.n = 100*Compressor.PVTsimSTD_Rel.n;
[Compressor.PVTsimLower.n, Compressor.PVTsimUpper.n] = coverage(Compressor.PVTsim.n,95);

% Polytropic head
Compressor.PVTsimSTD.PolyHead = std(Compressor.PVTsim.PolyHead);
Compressor.PVTsimSTD_Rel.PolyHead = Compressor.PVTsimSTD.PolyHead./mean(Compressor.PVTsim.PolyHead);
Compressor.PVTsimSTD_Perc.PolyHead = 100*Compressor.PVTsimSTD_Rel.PolyHead;
[Compressor.PVTsimLower.PolyHead, Compressor.PVTsimUpper.PolyHead] = coverage(Compressor.PVTsim.PolyHead,95);
% Actual head thermodynamic

Compressor.PVTsimSTD.ActHead_thermo = std(Compressor.PVTsim.ActHead_thermo);
Compressor.PVTsimSTD_Rel.ActHead_thermo = Compressor.PVTsimSTD.ActHead_thermo./mean(Compressor.PVTsim.ActHead_thermo);
Compressor.PVTsimSTD_Perc.ActHead_thermo = 100*Compressor.PVTsimSTD_Rel.ActHead_thermo;
[Compressor.PVTsimLower.ActHead_thermo, Compressor.PVTsimUpper.ActHead_thermo] = coverage(Compressor.PVTsim.ActHead_thermo,95);
% Actual head mechanical

Compressor.PVTsimSTD.ActHead_Mech = std(Compressor.PVTsim.ActHead_Mech);
Compressor.PVTsimSTD_Rel.ActHead_Mech = Compressor.PVTsimSTD.ActHead_Mech./mean(Compressor.PVTsim.ActHead_Mech);
Compressor.PVTsimSTD_Perc.ActHead_Mech = 100*Compressor.PVTsimSTD_Rel.ActHead_Mech;
[Compressor.PVTsimLower.ActHead_Mech, Compressor.PVTsimUpper.ActHead_Mech] = coverage(Compressor.PVTsim.ActHead_Mech,95);

% Power themodynamic
Compressor.PVTsimSTD.Power_thermo = std(Compressor.PVTsim.Power_thermo);
Compressor.PVTsimSTD_Rel.Power_thermo = Compressor.PVTsimSTD.Power_thermo./mean(Compressor.PVTsim.Power_thermo);
Compressor.PVTsimSTD_Perc.Power_thermo = 100*Compressor.PVTsimSTD_Rel.Power_thermo;
[Compressor.PVTsimLower.Power_thermo, Compressor.PVTsimUpper.Power_thermo] = coverage(Compressor.PVTsim.Power_thermo,95);

% Power mechanical
Compressor.PVTsimSTD.Power_Mech = std(Compressor.PVTsim.Power_Mech);
Compressor.PVTsimSTD_Rel.Power_Mech = Compressor.PVTsimSTD.Power_Mech./mean(Compressor.PVTsim.Power_Mech);
Compressor.PVTsimSTD_Perc.Power_Mech = 100*Compressor.PVTsimSTD_Rel.Power_Mech;
[Compressor.PVTsimLower.Power_Mech, Compressor.PVTsimUpper.Power_Mech] = coverage(Compressor.PVTsim.Power_Mech,95);

% Efficiency themodynamic
Compressor.PVTsimSTD.PolyEff_thermo = std(Compressor.PVTsim.PolyEff_thermo);
Compressor.PVTsimSTD_Rel.PolyEff_thermo = Compressor.PVTsimSTD.PolyEff_thermo./mean(Compressor.PVTsim.PolyEff_thermo);
Compressor.PVTsimSTD_Perc.PolyEff_thermo = 100*Compressor.PVTsimSTD_Rel.PolyEff_thermo;
[Compressor.PVTsimLower.PolyEff_thermo, Compressor.PVTsimUpper.PolyEff_thermo] = coverage(Compressor.PVTsim.PolyEff_thermo,95);

% Efficiency mechanical
Compressor.PVTsimSTD.PolyEff_Mech = std(Compressor.PVTsim.PolyEff_Mech);
Compressor.PVTsimSTD_Rel.PolyEff_Mech = Compressor.PVTsimSTD.PolyEff_Mech./mean(Compressor.PVTsim.PolyEff_Mech);
Compressor.PVTsimSTD_Perc.PolyEff_Mech = 100*Compressor.PVTsimSTD_Rel.PolyEff_Mech;
[Compressor.PVTsimLower.PolyEff_Mech, Compressor.PVTsimUpper.PolyEff_Mech] = coverage(Compressor.PVTsim.PolyEff_Mech,95);

% Impeller Tip Speed
Compressor.PVTsimSTD.ImpTipSpeed = std(Compressor.PVTsim.ImpTipSpeed);
Compressor.PVTsimSTD_Rel.ImpTipSpeed = Compressor.PVTsimSTD.ImpTipSpeed./mean(Compressor.PVTsim.ImpTipSpeed);
Compressor.PVTsimSTD_Perc.ImpTipSpeed = 100*Compressor.PVTsimSTD_Rel.ImpTipSpeed;
[Compressor.PVTsimLower.ImpTipSpeed, Compressor.PVTsimUpper.ImpTipSpeed] = coverage(Compressor.PVTsim.ImpTipSpeed,95);

% Actual flow
Compressor.PVTsimSTD.ActualFlow = std(Compressor.PVTsim.ActualFlow);
Compressor.PVTsimSTD_Rel.ActualFlow = Compressor.PVTsimSTD.ActualFlow./mean(Compressor.PVTsim.ActualFlow);
Compressor.PVTsimSTD_Perc.ActualFlow = 100*Compressor.PVTsimSTD_Rel.ActualFlow;
[Compressor.PVTsimLower.ActualFlow, Compressor.PVTsimUpper.ActualFlow] = coverage(Compressor.PVTsim.ActualFlow,95);

% phi dry
Compressor.PVTsimSTD.FlowCoeffDry = std(Compressor.PVTsim.FlowCoeffDry);
Compressor.PVTsimSTD_Rel.FlowCoeffDry = Compressor.PVTsimSTD.FlowCoeffDry./mean(Compressor.PVTsim.FlowCoeffDry);
Compressor.PVTsimSTD_Perc.FlowCoeffDry = 100*Compressor.PVTsimSTD_Rel.FlowCoeffDry;
[Compressor.PVTsimLower.FlowCoeffDry, Compressor.PVTsimUpper.FlowCoeffDry] = coverage(Compressor.PVTsim.FlowCoeffDry,95);

% phi wet
Compressor.PVTsimSTD.FlowCoeffWet = std(Compressor.PVTsim.FlowCoeffWet);
Compressor.PVTsimSTD_Rel.FlowCoeffWet = Compressor.PVTsimSTD.FlowCoeffWet./mean(Compressor.PVTsim.FlowCoeffWet);
Compressor.PVTsimSTD_Perc.FlowCoeffWet = 100*Compressor.PVTsimSTD_Rel.FlowCoeffWet;
[Compressor.PVTsimLower.FlowCoeffWet, Compressor.PVTsimUpper.FlowCoeffWet] = coverage(Compressor.PVTsim.FlowCoeffWet,95);

% Ma dry (mach number dry)
Compressor.PVTsimSTD.MachNumberDry = std(Compressor.PVTsim.MachNumberDry);
Compressor.PVTsimSTD_Rel.MachNumberDry = Compressor.PVTsimSTD.MachNumberDry./mean(Compressor.PVTsim.MachNumberDry);
Compressor.PVTsimSTD_Perc.MachNumberDry = 100*Compressor.PVTsimSTD_Rel.MachNumberDry;
[Compressor.PVTsimLower.MachNumberDry, Compressor.PVTsimUpper.MachNumberDry] = coverage(Compressor.PVTsim.MachNumberDry,95);

% Ma wet (mach number wet based woods Hundseid multphase speed of sound) 
Compressor.PVTsimSTD.MachNumberWet = std(Compressor.PVTsim.MachNumberWet);
Compressor.PVTsimSTD_Rel.MachNumberWet = Compressor.PVTsimSTD.MachNumberWet./mean(Compressor.PVTsim.MachNumberWet);
Compressor.PVTsimSTD_Perc.MachNumberWet = 100*Compressor.PVTsimSTD_Rel.MachNumberWet;
[Compressor.PVTsimLower.MachNumberWet, Compressor.PVTsimUpper.MachNumberWet] = coverage(Compressor.PVTsim.MachNumberWet,95);

% Machine/tip Reynolds number Dry
Compressor.PVTsimSTD.Re_M_Dry = std(Compressor.PVTsim.Re_M_Dry);
Compressor.PVTsimSTD_Rel.Re_M_Dry = Compressor.PVTsimSTD.Re_M_Dry./mean(Compressor.PVTsim.Re_M_Dry);
Compressor.PVTsimSTD_Perc.Re_M_Dry = 100*Compressor.PVTsimSTD_Rel.Re_M_Dry;
[Compressor.PVTsimLower.Re_M_Dry, Compressor.PVTsimUpper.Re_M_Dry] = coverage(Compressor.PVTsim.Re_M_Dry,95);

% Machine/tip Reynolds number Wet
Compressor.PVTsimSTD.Re_M_wet_MacAdams = std(Compressor.PVTsim.Re_M_wet_MacAdams);
Compressor.PVTsimSTD_Rel.Re_M_wet_MacAdams = Compressor.PVTsimSTD.Re_M_wet_MacAdams./mean(Compressor.PVTsim.Re_M_wet_MacAdams);
Compressor.PVTsimSTD_Perc.Re_M_wet_MacAdams = 100*Compressor.PVTsimSTD_Rel.Re_M_wet_MacAdams;

Compressor.PVTsimSTD.Re_M_wet_Bingham = std(Compressor.PVTsim.Re_M_wet_Bingham);
Compressor.PVTsimSTD_Rel.Re_M_wet_Bingham = Compressor.PVTsimSTD.Re_M_wet_Bingham./mean(Compressor.PVTsim.Re_M_wet_Bingham);
Compressor.PVTsimSTD_Perc.Re_M_wet_Bingham = 100*Compressor.PVTsimSTD_Rel.Re_M_wet_Bingham;

Compressor.PVTsimSTD.Re_M_wet_Cicchitti = std(Compressor.PVTsim.Re_M_wet_Cicchitti);
Compressor.PVTsimSTD_Rel.Re_M_wet_Cicchitti = Compressor.PVTsimSTD.Re_M_wet_Cicchitti./mean(Compressor.PVTsim.Re_M_wet_Cicchitti);
Compressor.PVTsimSTD_Perc.Re_M_wet_Cicchitti = 100*Compressor.PVTsimSTD_Rel.Re_M_wet_Cicchitti;

Compressor.PVTsimSTD.Re_M_wet_My = std(Compressor.PVTsim.Re_M_wet_My);
Compressor.PVTsimSTD_Rel.Re_M_wet_My = Compressor.PVTsimSTD.Re_M_wet_My./mean(Compressor.PVTsim.Re_M_wet_My);
Compressor.PVTsimSTD_Perc.Re_M_wet_My = 100*Compressor.PVTsimSTD_Rel.Re_M_wet_My;

[Compressor.PVTsimLower.Re_M_wet_MacAdams, Compressor.PVTsimUpper.Re_M_wet_MacAdams] = coverage(Compressor.PVTsim.Re_M_wet_MacAdams,95);
[Compressor.PVTsimLower.Re_M_wet_Bingham, Compressor.PVTsimUpper.Re_M_wet_Bingham] = coverage(Compressor.PVTsim.Re_M_wet_Bingham,95);
[Compressor.PVTsimLower.Re_M_wet_Cicchitti, Compressor.PVTsimUpper.Re_M_wet_Cicchitti] = coverage(Compressor.PVTsim.Re_M_wet_Cicchitti,95);
[Compressor.PVTsimLower.Re_M_wet_My, Compressor.PVTsimUpper.Re_M_wet_My] = coverage(Compressor.PVTsim.Re_M_wet_My,95);

%Specific volume ratio machine
Compressor.PVTsimSTD.SpecVolRatio = std(Compressor.PVTsim.SpecVolRatio);
Compressor.PVTsimSTD_Rel.SpecVolRatio = Compressor.PVTsimSTD.SpecVolRatio./mean(Compressor.PVTsim.SpecVolRatio);
Compressor.PVTsimSTD_Perc.SpecVolRatio = 100*Compressor.PVTsimSTD_Rel.SpecVolRatio;
[Compressor.PVTsimLower.SpecVolRatio, Compressor.PVTsimUpper.SpecVolRatio] = coverage(Compressor.PVTsim.SpecVolRatio,95);

%Density Ratio
Compressor.PVTsimSTD.Dr = std(Compressor.PVTsim.Dr);
Compressor.PVTsimSTD_Rel.Dr = Compressor.PVTsimSTD.Dr./mean(Compressor.PVTsim.Dr);
Compressor.PVTsimSTD_Perc.Dr = 100*Compressor.PVTsimSTD_Rel.Dr;
[Compressor.PVTsimLower.Dr, Compressor.PVTsimUpper.Dr] = coverage(Compressor.PVTsim.Dr,95);

%LMF
Compressor.PVTsimSTD.LMF = std(Compressor.PVTsim.LMF);
Compressor.PVTsimSTD_Rel.LMF = Compressor.PVTsimSTD.LMF./mean(Compressor.PVTsim.LMF);
Compressor.PVTsimSTD_Perc.LMF = 100*Compressor.PVTsimSTD_Rel.LMF;
[Compressor.PVTsimLower.LMF, Compressor.PVTsimUpper.LMF] = coverage(Compressor.PVTsim.LMF,95);

%GVF
Compressor.PVTsimSTD.GVF = std(Compressor.PVTsim.GVF);
Compressor.PVTsimSTD_Rel.GVF = Compressor.PVTsimSTD.GVF./mean(Compressor.PVTsim.GVF);
Compressor.PVTsimSTD_Perc.GVF = 100*Compressor.PVTsimSTD_Rel.GVF;
[Compressor.PVTsimLower.GVF, Compressor.PVTsimUpper.GVF] = coverage(Compressor.PVTsim.GVF,95);
end

%% Find upper and lower covering limits for Y, this simplified version is ok for large N
function [Lower, Upper] = coverage(Y,Coverage)
% find upper and lower covering limits for Y
% GUM supplement 1 chapter 7.2.2
    Y = sort(Y);
    N = length(Y);
    coverage = Coverage/100; %convert coverage from percent to frac
    idLower = round(((1-coverage)/2)*N);
    idUpper = round(((1+coverage)/2)*N);
    if idLower == 0
        idLower = 1;
    end
    Lower = Y(idLower);
    Upper = Y(idUpper);
end


