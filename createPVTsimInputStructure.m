function [S,Constants] = createPVTsimInputStructure(IM,Config)
 
%% Ranges for compositional data
GasCompRange = Config.GasCompRange;
OilCompRange = Config.OilCompRange;
LastFluidIndex = Config.LastFluidIndex;


%% Create Structure based on input Matrix
S.Gas_X = IM(GasCompRange,:); 
S.Gas_M = IM(GasCompRange(end)+1:GasCompRange(end)+5,:);
S.Gas_D = IM(GasCompRange(end)+6:GasCompRange(end)+10,:);
S.Oil_X = IM(OilCompRange,:);
S.Oil_M = IM(OilCompRange(end)+1:OilCompRange(end)+5,:);
S.Oil_D = IM(OilCompRange(end)+6:OilCompRange(end)+10,:);

S.P_Sep = IM(LastFluidIndex+1,:);
S.T_Sep = IM(LastFluidIndex+2,:);

S.P_Orifice = IM(LastFluidIndex+3,:);
S.T_Orifice = IM(LastFluidIndex+4,:);
S.dP_Orifice = IM(LastFluidIndex+5,:);

S.P_In = IM(LastFluidIndex+6,:);
S.T_In = IM(LastFluidIndex+7,:);

S.P_Out = IM(LastFluidIndex+8,:);
S.T_Out = IM(LastFluidIndex+9,:);

S.m_Gas = IM(LastFluidIndex+10,:);
S.m_Oil = IM(LastFluidIndex+11,:);
S.m_Water = IM(LastFluidIndex+12,:);

S.m_Oil_Cori = IM(LastFluidIndex+13,:);
S.m_Water_Cori = IM(LastFluidIndex+14,:);

S.PowerTorque = IM(LastFluidIndex+15,:);%1400; %kW

S.Speed = IM(LastFluidIndex+16,:);
S.ImpellerDiam = IM(LastFluidIndex+17,:);

%%Gerg index
S.GERG_idx = 2:14;
S.GERG_C10pussIdx = 15:23;



%% Constants
Constants.D_FP3 = 0.2164995;
Constants.D_FP4 = 0.1749755;
Constants.GearRatio = 185/28;
Constants.atm = 1.01325; %athmospheric pressure bara 
Constants.R = 8.314459848;  %Gas Constant 
Constants.GuardArea = pi*1.7^2;% m^2