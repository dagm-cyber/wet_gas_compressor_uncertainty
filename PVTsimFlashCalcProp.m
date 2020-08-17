function [Output] =  PVTsimFlashCalcProp(Output,FlashOutput,i)
%% Extract mixtrure and phase properties for the PVTsim flash output objeect

R = 8.31446261815324;

hasPhaseTypeGas = false;
hasPhaseTypeOil = false;
hasPhaseTypeAqueous = false;

if ~isfield(Output,'outputUnits')
    %% Output units
    
    VariableMix = {'mixMolarMass';'mixMassEnthalpy';'mixMassEntropy';'mixMolarVolume';'mixDensity';'mixZfactor';'mixEnthalpy';'mixEntropy';'mixCp';'mixCv';'mixKappa';'mixViscosity';'mixThermalConductivity'};
    UnitMix = {'g/mol';'kJ/kg';'kJ/kgK';'m3/mole';'kg/m3';'-';'J/mol';'J/mol K';'J/mol K';'J/mol K';'-';'Pa*s';'W/mK'};
    Variablegas = {'gasFractionc';'gasMolarVolume';'gasVolumeFraction';'gasDensity';'gasZ';'gasMolarMass';'gasEnthalpy';'gasEntropy';'gasWtFraction';'gasCp';'gasCv';'gasKappa';'gasViscosity';'gasThermalConductivity';'gasSoundSpeed';'gasJouleThomsonCoefficient'};
    Unitgas = {'-';'m3/mole';'-';'kg/m3';'-';'g/mol';'J/mol';'J/mol K';'-';'J/mol K';'J/mol K';'-';'Pa*s';'W/mK';'m/s';'K/Pa'};
    Variableoil = {'oilFractionc';'oilMolarVolume';'oilVolumeFraction';'oilDensity';'oilZ';'oilMolarMass';'oilEnthalpy';'oilEntropy';'oilWtFraction';'oilCp';'oilCv';'oilKappa';'oilViscosity';'oilThermalConductivity';'oilSoundSpeed';'oilJouleThomsonCoefficient'};
    Unitoil = {'-';'m3/mole';'-';'kg/m3';'-';'g/mol';'J/mol';'J/mol K';'-';'J/mol K';'J/mol K';'-';'Pa*s';'W/mK';'m/s';'K/Pa'};
    Variablewater = {'waterFractionc';'waterMolarVolume';'waterVolumeFraction';'waterDensity';'waterZ';'waterMolarMass';'waterEnthalpy';'waterEntropy';'waterWtFraction';'waterCp';'waterCv';'waterKappa';'waterViscosity';'waterThermalConductivity';'waterSoundSpeed';'waterJouleThomsonCoefficient'};
    Unitwater = {'-';'m3/mole';'-';'kg/m3';'-';'g/mol';'J/mol';'J/mol K';'-';'J/mol K';'J/mol K';'-';'Pa*s';'W/mK';'m/s';'K/Pa'};
    
    Variable = [VariableMix;Variablegas;Variableoil;Variablewater];
    Unit = [UnitMix;Unitgas;Unitoil;Unitwater];
    T = table(Variable,Unit,'RowNames',Variable);
    Output.outputUnits = T;

end

%% Extract mixture properties from PVTsim Flash output object

    Output.mixMolarMass(i) = FlashOutput.Mixture.Properties.MolecularWeight; %'g/mol'    
    Output.mixMassEnthalpy(i) = FlashOutput.Mixture.Properties.Enthalpy/Output.mixMolarMass(i);
    Output.mixMassEntropy(i) = FlashOutput.Mixture.Properties.Entropy/Output.mixMolarMass(i);
    %Output.numberOfPhases(i) = fluid.getNumberOfPhases();
    Output.mixMolarVolume(i) =  FlashOutput.Mixture.Properties.MolarVolume; %1.0 / fluid.getDensity("mol/m3"); 
    Output.mixDensity(i) = FlashOutput.Mixture.Properties.Density; %fluid.getDensity("kg/m3");
    Output.mixZfactor(i) = FlashOutput.Mixture.Properties.ZFactor;%fluid.getZ();
    Output.mixEnthalpy(i) = FlashOutput.Mixture.Properties.Enthalpy;%fluid.getEnthalpy("Jmol");
    Output.mixEntropy(i) = FlashOutput.Mixture.Properties.Entropy;%fluid.getEntropy("JmolK");
    Output.mixCp(i) = FlashOutput.Mixture.Properties.Cp;%fluid.getCp("J/molK");
    Output.mixCv(i) = FlashOutput.Mixture.Properties.Cv;%fluid.getCv("J/molK");
    Output.mixCp_Cv(i) = FlashOutput.Mixture.Properties.Kappa;%fluid.getKappa();
    Output.mixViscosity(i) = FlashOutput.Mixture.Properties.Viscosity;%fluid.getViscosity("kg/msec");
    Output.mixThermalConductivity(i) = FlashOutput.Mixture.Properties.ThermalConductivity;%fluid.getConductivity("W/mK");
        
%% Extract gas oil and water phase properties from PVTsim Flash output object

    for phase = 0:1:(FlashOutput.NumberofPhases-1)
        if strcmp(FlashOutput.Phase(phase).properties.phasetype,'PvtsFlashPhaseTypeEnum_PvtsFlashPhaseTypeGas')
        
        PhaseNumberGas = phase; %Used for interphase properties  
        hasPhaseTypeGas = true;
            

        Output.gasMoleFraction(i) = FlashOutput.Phase(phase).Properties.MolePercent/100; % fluid.getMoleFraction(phaseNumber) * 100; 
        Output.gasMolarVolume(i) = FlashOutput.Phase(phase).Properties.MolarVolume;% 1.0 / fluid.getPhase(phaseNumber).getDensity("mol/m3");
        Output.gasVolumeFraction(i) = FlashOutput.Phase(phase).Properties.VolumePercent/100; %fluid.getCorrectedVolumeFraction(phaseNumber) * 100; 
        Output.gasDensity(i) = FlashOutput.Phase(phase).Properties.Density; %kg/m3 %fluid.getPhase(phaseNumber).getDensity("kg/m3");
        Output.gasZ(i) = FlashOutput.Phase(phase).Properties.ZFactor; %fluid.getPhase(phaseNumber).getZ();
        Output.gasMolarMass(i)= FlashOutput.Phase(phase).Properties.MolecularWeight; %fluid.getPhase(phaseNumber).getMolarMass() * 1000;
        Output.gasEnthalpy(i) = FlashOutput.Phase(phase).Properties.Enthalpy; %J/mol %fluid.getPhase(phaseNumber).getEnthalpy("J/mol");
        Output.gasMassEnthalpy(i) =  Output.gasEnthalpy(i)/Output.gasMolarMass(i); %kJ/kg
        Output.gasEntropy(i) = FlashOutput.Phase(phase).Properties.Entropy;% J/mol K
        Output.gasWtFraction(i) = FlashOutput.Phase(phase).Properties.WeightPercent/100; % fluid.getWtFraction(phaseNumber) * 100;
        Output.gasCp(i) = FlashOutput.Phase(phase).Properties.Cp; %'J/mol K'
        Output.gasCv(i) = FlashOutput.Phase(phase).Properties.Cv; %'J/mol K'
        Output.gasCp_Cv(i) = FlashOutput.Phase(phase).Properties.Kappa; % fluid.getPhase(phaseNumber).getGamma();
        Output.gasViscosity(i) = FlashOutput.Phase(phase).Properties.Viscosity; %Pa*s  %fluid.getPhase(phaseNumber).getViscosity("kg/msec");
        Output.gasThermalConductivity(i) = FlashOutput.Phase(phase).Properties.ThermalConductivity; %W/mK % fluid.getPhase(phaseNumber).getConductivity("W/mK");
        Output.gasSoundSpeed(i) = FlashOutput.Phase(phase).Properties.VelocityOfSound; %m/s%fluid.getPhase(phaseNumber).getSoundSpeed();
        Output.gasJouleThomsonCoefficient(i) = FlashOutput.Phase(phase).Properties.JTCoefficient; %K/Pa ;% fluid.getPhase(phaseNumber).getJouleThomsonCoefficient() / 1e5;
        Output.gasKappa(i) = Output.gasSoundSpeed(i)^2*Output.gasMolarMass(i) /(Output.gasZ(i)*R*FlashOutput.Phase(phase).Properties.Temperature*1000);   
     end
    
    if strcmp(FlashOutput.Phase(phase).properties.phasetype,'PvtsFlashPhaseTypeEnum_PvtsFlashPhaseTypeLiquidHC')
     
          PhaseNumberOil = phase; %%Used for interphase properties 
          hasPhaseTypeOil = true;
          

        Output.oilMoleFraction(i) = FlashOutput.Phase(phase).Properties.MolePercent/100; % fluid.getMoleFraction(phaseNumber) * 100; 
        Output.oilMolarVolume(i) = FlashOutput.Phase(phase).Properties.MolarVolume;% 1.0 / fluid.getPhase(phaseNumber).getDensity("mol/m3");
        Output.oilVolumeFraction(i) = FlashOutput.Phase(phase).Properties.VolumePercent/100; %fluid.getCorrectedVolumeFraction(phaseNumber) * 100; 
        Output.oilDensity(i) = FlashOutput.Phase(phase).Properties.Density; %kg/m3 %fluid.getPhase(phaseNumber).getDensity("kg/m3");
        Output.oilZ(i) = FlashOutput.Phase(phase).Properties.ZFactor; %fluid.getPhase(phaseNumber).getZ();
        Output.oilMolarMass(i)= FlashOutput.Phase(phase).Properties.MolecularWeight; %fluid.getPhase(phaseNumber).getMolarMass() * 1000;
        Output.oilEnthalpy(i) = FlashOutput.Phase(phase).Properties.Enthalpy; %J/mol %fluid.getPhase(phaseNumber).getEnthalpy("J/mol");
        Output.oilMassEnthalpy(i) =  Output.oilEnthalpy(i)/Output.oilMolarMass(i); %kJ/kg
        Output.oilEntropy(i) = FlashOutput.Phase(phase).Properties.Entropy;% J/mol K
        Output.oilWtFraction(i) = FlashOutput.Phase(phase).Properties.WeightPercent/100; % fluid.getWtFraction(phaseNumber) * 100;
        Output.oilCp(i) = FlashOutput.Phase(phase).Properties.Cp; %'J/mol K'
        Output.oilCv(i) = FlashOutput.Phase(phase).Properties.Cv; %'J/mol K'
        Output.oilCp_Cv(i) = FlashOutput.Phase(phase).Properties.Kappa; % fluid.getPhase(phaseNumber).getGamma();
        Output.oilViscosity(i) = FlashOutput.Phase(phase).Properties.Viscosity; %Pa*s  %fluid.getPhase(phaseNumber).getViscosity("kg/msec");
        Output.oilThermalConductivity(i) = FlashOutput.Phase(phase).Properties.ThermalConductivity; %W/mK % fluid.getPhase(phaseNumber).getConductivity("W/mK");
        Output.oilSoundSpeed(i) = FlashOutput.Phase(phase).Properties.VelocityOfSound; %m/s%fluid.getPhase(phaseNumber).getSoundSpeed();
        Output.oilJouleThomsonCoefficient(i) = FlashOutput.Phase(phase).Properties.JTCoefficient; %K/Pa ;% fluid.getPhase(phaseNumber).getJouleThomsonCoefficient() / 1e5;
        

     end
    
     if strcmp(FlashOutput.Phase(phase).properties.phasetype,'PvtsFlashPhaseTypeEnum_PvtsFlashPhaseTypeAqueous')
        
        PhaseNumberAqueous = phase; %Used for interphase properties 
        hasPhaseTypeAqueous = true;
        
        Output.waterMoleFraction(i) = FlashOutput.Phase(phase).Properties.MolePercent/100; % fluid.getMoleFraction(phaseNumber) * 100; 
        Output.waterMolarVolume(i) = FlashOutput.Phase(phase).Properties.MolarVolume;% 1.0 / fluid.getPhase(phaseNumber).getDensity("mol/m3");
        Output.waterVolumeFraction(i) = FlashOutput.Phase(phase).Properties.VolumePercent/100; %fluid.getCorrectedVolumeFraction(phaseNumber) * 100; 
        Output.waterDensity(i) = FlashOutput.Phase(phase).Properties.Density; %kg/m3 %fluid.getPhase(phaseNumber).getDensity("kg/m3");
        Output.waterZ(i) = FlashOutput.Phase(phase).Properties.ZFactor; %fluid.getPhase(phaseNumber).getZ();
        Output.waterMolarMass(i)= FlashOutput.Phase(phase).Properties.MolecularWeight; %fluid.getPhase(phaseNumber).getMolarMass() * 1000;
        Output.waterEnthalpy(i) = FlashOutput.Phase(phase).Properties.Enthalpy; %J/mol %fluid.getPhase(phaseNumber).getEnthalpy("J/mol");
        Output.waterMassEnthalpy(i) =  Output.waterEnthalpy(i)/Output.waterMolarMass(i); %kJ/kg
        Output.waterEntropy(i) = FlashOutput.Phase(phase).Properties.Entropy;% J/mol K
        Output.waterWtFraction(i) = FlashOutput.Phase(phase).Properties.WeightPercent/100; % fluid.getWtFraction(phaseNumber) * 100;
        Output.waterCp(i) = FlashOutput.Phase(phase).Properties.Cp; %'J/mol K'
        Output.waterCv(i) = FlashOutput.Phase(phase).Properties.Cv; %'J/mol K'
        Output.waterCp_Cv(i) = FlashOutput.Phase(phase).Properties.Kappa; % fluid.getPhase(phaseNumber).getGamma();
        Output.waterViscosity(i) = FlashOutput.Phase(phase).Properties.Viscosity; %Pa*s  %fluid.getPhase(phaseNumber).getViscosity("kg/msec");
        Output.waterThermalConductivity(i) = FlashOutput.Phase(phase).Properties.ThermalConductivity; %W/mK % fluid.getPhase(phaseNumber).getConductivity("W/mK");
        Output.waterSoundSpeed(i) = FlashOutput.Phase(phase).Properties.VelocityOfSound; %m/s%fluid.getPhase(phaseNumber).getSoundSpeed();
        Output.waterJouleThomsonCoefficient(i) = FlashOutput.Phase(phase).Properties.JTCoefficient; %K/Pa ;% fluid.getPhase(phaseNumber).getJouleThomsonCoefficient() / 1e5;
    end
    end 
   
    if ~hasPhaseTypeGas
        Output.gasMoleFraction(i) = 0;
        Output.gasMolarVolume(i) = nan;
        Output.gasVolumeFraction(i) = 0;
        Output.gasDensity(i) = nan;
        Output.gasZ(i) = nan;
        Output.gasMolarMass(i)= nan;
        Output.gasEnthalpy(i) = nan;
        Output.gasMassEnthalpy(i) =  nan;
        Output.gasEntropy(i) = nan;
        Output.gasWtFraction(i) = 0;
        Output.gasCp(i) = nan;
        Output.gasCv(i) = nan;
        Output.gasCp_Cv(i) = nan;
        Output.gasViscosity(i) = nan;
        Output.gasThermalConductivity(i) = nan;
        Output.gasSoundSpeed(i) = nan;
        Output.gasJouleThomsonCoefficient(i) = nan;
        Output.gasKappa(i) = nan;

    end
    
    if ~hasPhaseTypeOil  
        Output.oilMoleFraction(i) = 0;
        Output.oilMolarVolume(i) = nan;
        Output.oilVolumeFraction(i) = 0;
        Output.oilDensity(i) = nan;
        Output.oilZ(i) = nan;
        Output.oilMolarMass(i)= nan;
        Output.oilEnthalpy(i) = nan;
        Output.oilMassEnthalpy(i) =  nan;
        Output.oilEntropy(i) = nan;
        Output.oilWtFraction(i) = 0;
        Output.oilCp(i) = nan;
        Output.oilCv(i) = nan;
        Output.oilCp_Cv(i) = nan;
        Output.oilViscosity(i) = nan;
        Output.oilThermalConductivity(i) = nan;
        Output.oilSoundSpeed(i) = nan;
        Output.oilJouleThomsonCoefficient(i) = nan;
        
    end
    
    if ~hasPhaseTypeAqueous
        
        Output.waterMoleFraction(i) = 0;
        Output.waterMolarVolume(i) = nan;
        Output.waterVolumeFraction(i) = 0;
        Output.waterDensity(i) = nan;
        Output.waterZ(i) = nan;
        Output.waterMolarMass(i)= nan;
        Output.waterEnthalpy(i) = nan;
        Output.waterMassEnthalpy(i) = nan;
        Output.waterEntropy(i) = nan;
        Output.waterWtFraction(i) = 0;
        Output.waterCp(i) = nan;
        Output.waterCv(i) = nan;
        Output.waterCp_Cv(i) = nan;
        Output.waterViscosity(i) = nan;
        Output.waterThermalConductivity(i) = nan;
        Output.waterSoundSpeed(i) = nan;
        Output.waterJouleThomsonCoefficient(i) = nan;
        
    end
        
    if hasPhaseTypeGas && hasPhaseTypeOil
        Output.gasoilSurfaceTension(i)  = FlashOutput.CrossPhase.SurfaceTension.PhasesSurfaceTension(PhaseNumberGas,PhaseNumberOil); %fluid.getInterfacialTension('gas','oil');
    else
        Output.gasoilSurfaceTension(i) = nan;
    end
            
    if hasPhaseTypeGas && hasPhaseTypeAqueous
        Output.gaswaterSurfaceTension(i)= FlashOutput.CrossPhase.SurfaceTension.PhasesSurfaceTension(PhaseNumberGas,PhaseNumberAqueous);  %fluid.getInterfacialTension('gas','aqueous');
    else
        Output.gaswaterSurfaceTension(i) = nan;
    end
        
    if hasPhaseTypeOil && hasPhaseTypeAqueous
        Output.oilwaterSurfaceTension(i) = FlashOutput.CrossPhase.SurfaceTension.PhasesSurfaceTension(PhaseNumberOil,PhaseNumberAqueous); %fluid.getInterfacialTension('oil','aqueous'); 
    else
        Output.oilwaterSurfaceTension(i) = nan;
    end

end