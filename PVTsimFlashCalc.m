function [FlashOutput] =  PVTsimFlashCalc(Comp,P,T,Flash,FlashInput,FlashOutput)
    
    Ncomp = Flash.FluidInfo.NumberofComponents;
    %convert pressure and temperature from barg and C to Pascal and Kelvin
    P_Pa = (P + 1.01325)*1E5;
    T_K = T+273.15;
    %flash_input.Specification1 = 101325  # Pressure in Pascal
    FlashInput.Specification1 = P_Pa;
    %flash_input.Specification2 = 288.15  # Temperature in Kelvin
    FlashInput.Specification2 = T_K;
    %update flash object with fluid composition
    for i=0:Ncomp-1
    %setter komposition
        Flash.FluidInfo.set('MixtureMolePercent', i, Comp(i+1));
    end
    Flash.Run(FlashInput, FlashOutput);



