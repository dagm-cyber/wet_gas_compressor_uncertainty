function [FP, LoopCalcData, IS, Constants, totCompNames] = PVTsimLoopCalc_IsolateF(IM,Config,ProgressText,sensitivity)
% Loop calculation for the monte Carlo simulation
% sensitivity: true for sensitivity false for MCM...true turne off perturbation of discharge coeff C in orifice calculation

%%
n = length(IM(1,:));
f = waitbar(0,ProgressText);
waitbar(0,f,ProgressText)

%% Create Input Structure
[IS,Constants] = createPVTsimInputStructure(IM, Config); 

%% Connect to PVTsim
PVTsimConnect;

%% Preallocate structures for Flash point data
    FP = struct; 
    FP.Sep.PVTsim.Result = struct; 
    FP.Orifice.PVTsim.Result = struct;
    FP.In.PVTsim.Result = struct;
    FP.Out.PVTsim.Result = struct;
    FP.Orifice.PVTsim.Result.GERG =struct;
    FP.Orifice.PVTsim.Result.GERG_Mod =struct;
    
%% Preallocate momory
LoopCalcIterationTime = zeros(1,n);
AccululativeIterationTime = zeros(1,n);
        
%% Create C10+ mixed fluid matrix

% Molecular weights of defined substances
M_const = [18.0153408050000;28.0135231020000;44.0097999570000;16.0428791050000;30.0698146820000;44.0967674260000;58.1236991880000;...
        58.1237068180000;72.1506423950000;72.1506347700000];

    GasX = [zeros(1,length(IS.Gas_X(1,:))); IS.Gas_X]/100; %Add water to gas composition vector
    GasX = GasX./sum(GasX); %Normalize
    GasM = [repmat(M_const,1,length(IS.Gas_M(1,:))); IS.Gas_M]; %Generate Mw array, by concatenating Mw of defined substance and Mw for C6,C7,... given in input matrix
    GasM_avarage = sum(GasX.*GasM); %Calculate average Mw
    Gas_n = IS.m_Gas./GasM_avarage; %Calculate number of moles gas
    
    OilX = [zeros(1,length(IS.Oil_X(1,:))); IS.Oil_X]./100; %Add water to composition vector
    if any(OilX(:)) %No normalization if no oil composition i present
        OilX = OilX./sum(OilX); %Normalize
    end
    OilM = [repmat(M_const,1,length(IS.Oil_M(1,:))); IS.Oil_M]; %%Generate Mw array, by concatenating Mw of defined substance and Mw for C6,C7,... given in input matrix
    OilM_avarage = sum(OilX.*OilM); %Calculate average Mw
    Oil_n = IS.m_Oil./OilM_avarage;%Calculate number of moles oil
    if ~any(OilX(:)) %write oil moles to zero if no oil is present
        Oil_n = zeros(1,length(Oil_n)); 
    end
    WaterX = [repmat(100,1,length(IS.Gas_X(1,:))); zeros(length(GasX(:,1))-1,length(GasX(1,:)))]./100; %Add nonwater components to water composition
    WaterX = WaterX./sum(WaterX); %Normalize
    WaterM = [repmat(M_const,1,length(IS.Gas_M(1,:)));IS.Oil_M];%dummy water Mw array, only Mw_water needed
    WaterM_avarage = sum(WaterX.*WaterM); %Calculate average Mw
    Water_n = IS.m_Water./WaterM_avarage; %Calculate number of moles water
    
    TotComp = Gas_n.*GasX + Oil_n.*OilX + Water_n.*WaterX; %Calculate new tot comp C10+
    TotComp = TotComp./sum(TotComp); %Normalize


%% Create characterized mixed fluid matrix

%Find C10 and distrubute C10+ according to distribution in PVTsim fluid characterization 
for i = 1:Flash.FluidInfo.NumberofComponents %find C10
    if isequal(Flash.FluidInfo.ComponentName(i-1),'C9')
        C10_Idx = i+1;
    end
end

j=0; %create C10+ fraction vector
for i = C10_Idx:Flash.FluidInfo.NumberofComponents
    j=j+1;
    C10_Plus_Frac(j) = Flash.FluidInfo.MixtureMolePercent(i-1); %#ok<AGROW>
end
C10_Plus_NormFrac = C10_Plus_Frac./sum(C10_Plus_Frac); %Normalize

TotCompChar = [TotComp(1:end-1,:); C10_Plus_NormFrac'*TotComp(end,:)]; %Calculate new tot comp  Characterized


%%  LoopCalc
for i= 1:length(IM(1,:))
    
    waitbar(i/n,f,ProgressText)    
    tic    


%% Separator flash

    FP.Sep.PVTsim.Pressure(i) = IS.P_Sep(i);
    FP.Sep.PVTsim.Temperature(i) = IS.T_Sep(i);
    FP_Sep_FlashOutput = PVTsimFlashCalc(TotCompChar(:,i),FP.Sep.PVTsim.Pressure(i),FP.Sep.PVTsim.Temperature(i)...
        ,Flash,FlashInput,FlashOutput); %Perform separator flash of total composition
    [FP.Sep.PVTsim.Result] =  PVTsimFlashCalcProp(FP.Sep.PVTsim.Result,FP_Sep_FlashOutput,i); %Collect properties data from the flash into a structure
    [FP.Sep.PVTsim.Result] =  PVTsimFlashCalcComp(FP.Sep.PVTsim.Result,FP_Sep_FlashOutput,Ncomp,i);%Collect phase compositions from the flash into a structure

 %% Orifice FP
    FP.Orifice.PVTsim.Temperature(i) = IS.T_Orifice(i);
    FP.Orifice.PVTsim.Pressure(i) = IS.P_Orifice(i);

    FP_Orifice_FlashOutput = PVTsimFlashCalc(FP.Sep.PVTsim.Result.gasComposition{i},FP.Orifice.PVTsim.Pressure(i)...
        ,FP.Orifice.PVTsim.Temperature(i),Flash,FlashInput,FlashOutput); %Perform orifice flash of the gas composition from the seperator flash
    [FP.Orifice.PVTsim.Result] =  PVTsimFlashCalcProp(FP.Orifice.PVTsim.Result,FP_Orifice_FlashOutput,i);%Collect properties data from the flash into a structure
    [FP.Orifice.PVTsim.Result] =  PVTsimFlashCalcComp(FP.Orifice.PVTsim.Result,FP_Orifice_FlashOutput,Ncomp,i);%Collect phase compositions from the flash into a structure

%% flow orifice

%determine when error in discharge coefficient is going to be used
    if Config.UseMeanInputAsFirstIteration && i==1
        errorInDischargeCoeff = false;
    elseif Config.UseErrorInOrificeDischargeCoeff
        errorInDischargeCoeff = true;
    elseif ~Config.UseErrorInOrificeDischargeCoeff
        errorInDischargeCoeff = false;
    end
    
%calculate flow rate and put data into structure    
    [m_dot,~,~,~,~,C] = flowOrifice(IS.dP_Orifice(i), FP.Orifice.PVTsim.Result.gasDensity(i), FP.Orifice.PVTsim.Result.gasKappa(i),...
        FP.Orifice.PVTsim.Pressure(i),FP.Orifice.PVTsim.Result.gasViscosity(i),errorInDischargeCoeff,sensitivity); %ISO 5167 - 2 calculation
    FP.Orifice.PVTsim.Result.totalMassFlow(i) = m_dot; 
    FP.Orifice.PVTsim.Result.DischargeCoeff(i) = C;
    FP.Orifice.PVTsim.Result.gasMassFlow(i) =FP.Orifice.PVTsim.Result.gasWtFraction(i)*FP.Orifice.PVTsim.Result.totalMassFlow(i); %Calculate gas mass flow
    FP.Orifice.PVTsim.Result.oilMassFlow(i) =FP.Orifice.PVTsim.Result.oilWtFraction(i)*FP.Orifice.PVTsim.Result.totalMassFlow(i); %Calculate oil mass flow
    FP.Orifice.PVTsim.Result.waterMassFlow(i) =FP.Orifice.PVTsim.Result.waterWtFraction(i)*FP.Orifice.PVTsim.Result.totalMassFlow(i); %Calculate water mass flow
    clear m_dot q_dot
    
%% GERG orifice calculation

%calculate flow rate based on GERG and put data into structure
%NB!! kappa in PVTsim is C_p/C_v and not isentropic constant

    [FP.Orifice.PVTsim.Result]=PVTsimGERG(FP.Orifice.PVTsim.Result,FP.Orifice.PVTsim.Result.mixComposition{i},...
        FP.Orifice.PVTsim.Pressure(i),FP.Orifice.PVTsim.Temperature(i),FlashGERG,FlashInputGERG,... 
        FlashOutputGERG,NcompGERG,IS.GERG_idx,IS.GERG_C10pussIdx,i); %Perform orifice Gerg flash collect properties and phase compositions
    
    [m_dot,~,~,~,~,C] = flowOrifice(IS.dP_Orifice(i), FP.Orifice.PVTsim.Result.GERG.gasDensity(i), FP.Orifice.PVTsim.Result.GERG.gasKappa(i),...
        FP.Orifice.PVTsim.Pressure(i),FP.Orifice.PVTsim.Result.GERG.gasViscosity(i),errorInDischargeCoeff,sensitivity);%ISO 5167 - 2 calculation
    FP.Orifice.PVTsim.Result.GERG.totalMassFlow(i) = m_dot;
    FP.Orifice.PVTsim.Result.GERG.DischargeCoeff(i) = C;
    clear m_dot q_dot
    
%% GERG Modified based on Z from the flash and the true Molecular Weigth
    FP.Orifice.PVTsim.Result.GERG_Mod.gasDensity(i) = FP.Orifice.PVTsim.Result.GERG.gasDensity(i)*FP.Orifice.PVTsim.Result.gasMolarMass(i)/FP.Orifice.PVTsim.Result.GERG.gasMolarMass(i);
    [m_dot,~,~,~,~,C] = flowOrifice(IS.dP_Orifice(i), FP.Orifice.PVTsim.Result.GERG_Mod.gasDensity(i), FP.Orifice.PVTsim.Result.GERG.gasKappa(i),...
        FP.Orifice.PVTsim.Pressure(i),FP.Orifice.PVTsim.Result.GERG.gasViscosity(i),errorInDischargeCoeff,sensitivity);
    FP.Orifice.PVTsim.Result.GERG_Mod.totalMassFlow(i) = m_dot;
    FP.Orifice.PVTsim.Result.GERG_Mod.DischargeCoeff(i) = C;
    clear m_dot q_dot

%% Mix Fluids upstream compressor
    comps = [FP.Sep.PVTsim.Result.gasComposition{i} FP.Sep.PVTsim.Result.oilComposition{i} FP.Sep.PVTsim.Result.waterComposition{i}]; %compositions of streams to be mixed   
    massflows = [FP.Orifice.PVTsim.Result.totalMassFlow(i) IS.m_Oil_Cori(i) IS.m_Water_Cori(i)]; %mass flows of streams to be mixed
    
    if sum(any(comps)) == 1
        newComp = comps(:,any(comps)); %If only one phase has composition no mix mass need to be done
    else
        newComp = PVTsimMixMassFlow(comps,massflows,Ncomp,Fluid);
    end

 %% Compressor inlet FP
    FP.In.PVTsim.Pressure(i) = IS.P_In(i);
    FP.In.PVTsim.Temperature(i) = IS.T_In(i);
    
    FP_In_FlashOutput = PVTsimFlashCalc(newComp, FP.In.PVTsim.Pressure(i),FP.In.PVTsim.Temperature(i)...
        ,Flash,FlashInput,FlashOutput); %Perform compressor inlet flash
    [FP.In.PVTsim.Result] =  PVTsimFlashCalcProp(FP.In.PVTsim.Result,FP_In_FlashOutput,i); %Collect properties data from the flash into a structure
    [FP.In.PVTsim.Result] =  PVTsimFlashCalcComp(FP.In.PVTsim.Result,FP_In_FlashOutput,Ncomp,i); %Collect phase compositions from the flash into a structure
    FP.In.PVTsim.Result.totalMassFlow(i) = sum(massflows);
    FP.In.PVTsim.Result.gasMassFlow(i) = FP.In.PVTsim.Result.gasWtFraction(i)*FP.In.PVTsim.Result.totalMassFlow(i); %Calculate gas mass flow
    FP.In.PVTsim.Result.oilMassFlow(i) = FP.In.PVTsim.Result.oilWtFraction(i)*FP.In.PVTsim.Result.totalMassFlow(i); %Calculate oil mass flow
    FP.In.PVTsim.Result.waterMassFlow(i) = FP.In.PVTsim.Result.waterWtFraction(i)*FP.In.PVTsim.Result.totalMassFlow(i); %Calculate water mass flow

 %% Compressor discharge FP
    FP.Out.PVTsim.Pressure(i) = IS.P_Out(i);
    FP.Out.PVTsim.Temperature(i) = IS.T_Out(i);
    
    FP_Out_FlashOutput = PVTsimFlashCalc(newComp, FP.Out.PVTsim.Pressure(i),FP.Out.PVTsim.Temperature(i)...
        ,Flash,FlashInput,FlashOutput); %Perform compressor discharge flash
    [FP.Out.PVTsim.Result] =  PVTsimFlashCalcProp(FP.Out.PVTsim.Result,FP_Out_FlashOutput,i); %Collect properties data from the flash into a structure
    [FP.Out.PVTsim.Result] =  PVTsimFlashCalcComp(FP.Out.PVTsim.Result,FP_Out_FlashOutput,Ncomp,i); %Collect phase compositions from the flash into a structure
    FP.Out.PVTsim.Result.totalMassFlow(i) = FP.In.PVTsim.Result.totalMassFlow(i);
    FP.Out.PVTsim.Result.gasMassFlow(i) =FP.Out.PVTsim.Result.gasWtFraction(i)*FP.Out.PVTsim.Result.totalMassFlow(i); %Calculate gas mass flow
    FP.Out.PVTsim.Result.oilMassFlow(i) =FP.Out.PVTsim.Result.oilWtFraction(i)*FP.Out.PVTsim.Result.totalMassFlow(i); %Calculate oil mass flow
    FP.Out.PVTsim.Result.waterMassFlow(i) =FP.Out.PVTsim.Result.waterWtFraction(i)*FP.Out.PVTsim.Result.totalMassFlow(i); %Calculate water mass flow


    clear comp massflows
    
    LoopCalcIterationTime(i) = toc;
    AccululativeIterationTime(i) = sum(LoopCalcIterationTime(1:i));

%     if ismember(100*i/n,0:10:100)
%         loopCalcPercentComplete = 100*i/n;
%     end
end
LoopCalcData.LoopCalcIterationTime = LoopCalcIterationTime;
LoopCalcData.AccululativeIterationTime = AccululativeIterationTime;

%% find component names of totatal composition
for k = 1:Flash.FluidInfo.NumberofComponents %loop over components
    totCompNames{k} = Flash.FluidInfo.ComponentName(k-1);
end
    close(f)
end