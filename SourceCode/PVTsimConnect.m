

%% connect to PVTsim_OpenStructure
Flash = actxserver('PVTsim_OpenStructure.Flash');
FlashGERG = actxserver('PVTsim_OpenStructure.Flash');
DataBase = actxserver('PVTsim_OpenStructure.FluidDepoWaxDataLayer');
Fluid = actxserver('PVTsim_OpenStructure.Fluid');
FluidGERG = actxserver('PVTsim_OpenStructure.Fluid');
FlashInput = actxserver('PVTsim_OpenStructure.FlashInput');
FlashOutput = actxserver('PVTsim_OpenStructure.FlashOutput');
FlashInputGERG = actxserver('PVTsim_OpenStructure.FlashInput');
FlashOutputGERG = actxserver('PVTsim_OpenStructure.FlashOutput');

% *** Database ***
% Connect to Database
%% Extract database fluids

    if isequal(Config.EOS,'SRK_HV')
        DataBase.Connect('PVTsimDatabaseSRK_HV.nfdb');
        % % Get Fluid
        DataBase.GetFluid1(Fluid,6); %mix oil, water and rich gas
        DataBase.GetFluid1(FluidGERG,5);%GERG fluid
    elseif isequal(Config.EOS,'SRK_CPA')
        DataBase.Connect('PVTsimDatabaseSRK_CPA.nfdb');
        % % Get Fluid
        DataBase.GetFluid1(Fluid,6); %mix oil, water and rich gas
        DataBase.GetFluid1(FluidGERG,5);%GERG fluid
    elseif isequal(Config.EOS,'PR_HV')
        DataBase.Connect('PVTsimDatabasePR_HV.nfdb');
        % % Get Fluid
        DataBase.GetFluid1(Fluid,6); %mix oil, water and rich gas
        DataBase.GetFluid1(FluidGERG,5);%GERG fluid
    elseif isequal(Config.EOS,'PR_CPA')
        DataBase.Connect('PVTsimDatabasePR_CPA.nfdb');
        % % Get Fluid
        DataBase.GetFluid1(Fluid,6); %mix oil, water and rich gas
        DataBase.GetFluid1(FluidGERG,5);%GERG fluid
    elseif isequal(Config.EOS,'SRK_Classic')
        DataBase.Connect('PVTsimDatabaseSRK_Classic.nfdb');
        % % Get Fluid
        DataBase.GetFluid1(Fluid,6); %mix oil, water and rich gas
        DataBase.GetFluid1(FluidGERG,5);%GERG fluid
    elseif isequal(Config.EOS,'PR_Classic')
        DataBase.Connect('PVTsimDatabasePR_Classic.nfdb');
        % % Get Fluid
        DataBase.GetFluid1(Fluid,6); %mix oil, water and rich gas
        DataBase.GetFluid1(FluidGERG,5);%GERG fluid
    end

    
DataBase.Disconnect; % Done with the database
Flash.Initialize; %Initialize Flash (this must be called once for the flash). %this line is time consuming
FlashGERG.Initialize;
Flash.SetReferenceFluid(Fluid); % Set reference fluid. %this line is time consuming
FlashGERG.SetReferenceFluid(FluidGERG);% Set reference fluid. %this line is time consuming

Ncomp = Fluid.Composition.NComponents; %Find the number of components in Fluid
NcompGERG = FluidGERG.Composition.NComponents;%Find the number of components in Fluid

Flash.ViscosityModel = 'PvtsFlashViscModelEnum_PvtsFlashViscModelCSP'; % Set Viscosity model to CSP
FlashGERG.ViscosityModel = 'PvtsFlashViscModelEnum_PvtsFlashViscModelCSP'; % Set Viscosity model to CSP

FlashInput.FlashMode = 'PvtsFlashModeEnum_PvtsFlashModePTAqueous'; % Set Flash type to aqueous PT-Flash
FlashInputGERG.FlashMode = 'PvtsFlashModeEnum_PvtsFlashModePTNonAqueous'; % Set Flash type to non-aqueous PT-Flash

FlashInput.CalculateDerivatives = 'PvtsFlashCalcDerivativesEnum_PvtsFlashCalcDerivativesNothing'; % Set no calculation of derivatives
FlashInputGERG.CalculateDerivatives = 'PvtsFlashCalcDerivativesEnum_PvtsFlashCalcDerivativesNothing'; % Set no calculation of derivatives

Fluid = Fluid.Clone;
Flash.SetReferenceFluid(Fluid); %this line is time consuming