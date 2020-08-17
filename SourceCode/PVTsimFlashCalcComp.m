function [Output] =  PVTsimFlashCalcComp(Output,FlashOutput,Ncomp,i)
%%Extract mixture and phase compositions

%% Extract mixture compositions
    Comp = zeros(1,Ncomp);
        for k=0:1:Ncomp-1
            Comp(k+1) =  FlashOutput.Mixture.ComponentMolePercent(k);    
        end
       
        Output.mixComposition{:,i} = (Comp/100)';

%% loop through phases and extract gas, oil and water phase compositions        
 hasPhaseTypeGas = false; 
 hasPhaseTypeOil = false;
 hasPhaseTypeAqueous = false;

    for phase = 0:1:(FlashOutput.NumberofPhases-1)
        if strcmp(FlashOutput.Phase(phase).properties.phasetype,'PvtsFlashPhaseTypeEnum_PvtsFlashPhaseTypeGas')
            hasPhaseTypeGas = true;
            Comp = zeros(1,Ncomp);
                for k=0:1:Ncomp-1
                    Comp(k+1) =  FlashOutput.Phase(phase).ComponentMolePercent(k);  
                end
            Output.gasComposition{:,i} = (Comp/100)';
 
        end
    
        if strcmp(FlashOutput.Phase(phase).properties.phasetype,'PvtsFlashPhaseTypeEnum_PvtsFlashPhaseTypeLiquidHC')
            hasPhaseTypeOil = true;  
            Comp = zeros(1,Ncomp);
                for k=0:1:Ncomp-1
                    Comp(k+1) =  FlashOutput.Phase(phase).ComponentMolePercent(k);  
                end
            Output.oilComposition{:,i} = (Comp/100)';
        end
    
        if strcmp(FlashOutput.Phase(phase).properties.phasetype,'PvtsFlashPhaseTypeEnum_PvtsFlashPhaseTypeAqueous')
            hasPhaseTypeAqueous = true;
            Comp = zeros(1,Ncomp);
                for k=0:1:Ncomp-1
                    Comp(k+1) =  FlashOutput.Phase(phase).ComponentMolePercent(k);  
                end
            Output.waterComposition{:,i} = (Comp/100)';
        end 
    end
    
    if ~hasPhaseTypeGas
         Output.gasComposition{:,i} = zeros(Ncomp,1);
    end
    
    if ~hasPhaseTypeOil
        Output.oilComposition{:,i} = zeros(Ncomp,1);
    end
    
    if ~hasPhaseTypeAqueous
        Output.waterComposition{:,i} = zeros(Ncomp,1);
    end
        
end
    




