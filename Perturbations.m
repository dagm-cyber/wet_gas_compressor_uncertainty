function [IM] = Perturbations(n, IA_num,Config)
 
f = waitbar(0,'Perturbations progress...');

%% Perturbe Process parameters (all than overwrite compositional data later in this function)

waitbar(0,f,'Perturbations Process parameters progress...')
mu = IA_num(:,2)';
sigma = IA_num(:,4)';
pertRetry = 0;
IM = zeros(length(IA_num(:,1)),n);
for i = 1:n
    for j = 1:length(IA_num(:,1))
    IM(j,i) = normrnd(mu(j),sigma(j));%InputMatrix
        while IM(j,i)<0 %recalculate random numbers if any negative numbers arise.
            IM(j,i) = normrnd(mu(j),sigma(j));
            pertRetry = pertRetry + 1; %count number of pertubation retries
        end
    end
end

%% Perturbe GAS
waitbar(0.33,f,'Perturbations Gas composition progress...')

GasRange = Config.GasCompRange;
muRG = IA_num(GasRange,2);
muRG = muRG./100;

%According to ISO 6974-3
sigmaRG = exp(-5.57015709 + 0.71898517 * log(muRG')); %equation is for composition in fraction and is equivalent to S_R = exp(-4.28 + 0.715 * log(muRG_Perc)) for composition in percent
[~, C1Idx] = max(muRG);
sigmaRG(C1Idx) = 0; 

PercUncertaintyRG = 100*sigmaRG./muRG';
  

% create nonnormalized gas composition distribution matix
UnNorm = zeros(length(muRG),n);
for i = 1:n
    UnNorm(:,i) = mvnrnd(muRG,sigmaRG.^2);
    while any(UnNorm(:,i) < 0)
        UnNorm(:,i) = mvnrnd(muRG,sigmaRG.^2);%recalculate random numbers if any negative numbers arise.
    end
end

%generate normalized gas composition matrix
Norm = zeros(length(muRG),n);
for i=1:n
    Norm(:,i) = UnNorm(:,i)./sum(UnNorm(:,i));
end
IM_RG = 100*Norm;

%% Perturbe DEAD OIL
waitbar(0.67,f,'Perturbations Oil composition progress...')
OilRange = Config.OilCompRange;
muDeadOil = IA_num(OilRange,2);
muDeadOil = muDeadOil./100;

PercUncertaintyDeadOil = [1 0 1 1 1 1 1 1 1 1 1 2 2 5];
sigmaDeadOil = (PercUncertaintyDeadOil.*muDeadOil')./100;

% create nonnormalized oil composition distribution matix
UnNorm = zeros(length(muDeadOil),n); 

for i = 1:n
    UnNorm(:,i) = mvnrnd(muDeadOil,sigmaDeadOil.^2);
    while any(UnNorm(:,i) < 0)
        UnNorm(:,i) = mvnrnd(muDeadOil,sigmaDeadOil.^2);%recalculate random numbers if any negative numbers arise.
    end
end

%generate normalized oil composition matrix
Norm = zeros(length(muDeadOil),n);
if any(UnNorm(:)) %No normalization if no oil composition i present
    for i=1:n
        Norm(:,i) = UnNorm(:,i)./sum(UnNorm(:,i));
    end
end

%% complete oil
Norm(isnan(Norm)) = 0;%remove nan
IM_DeadOil = 100*Norm;

%% Concatinate results to InputMatrix
IM(GasRange,:) = IM_RG;
IM(OilRange,:) = IM_DeadOil;

%% Put expectation value into IM if stated in Config
    if Config.UseMeanInputAsFirstIteration 
        IM = [mu' IM];
    end
    
if Config.UseConstLoadComposition %simulate 0 uncertainty in composition
    IM(Config.GasCompRange,:) = repmat(IM(Config.GasCompRange,1),[1,length(IM(1,:))]);
    IM(Config.OilCompRange,:) = repmat(IM(Config.OilCompRange,1),[1,length(IM(1,:))]);
end

    PerturbationsPerOK_value = pertRetry/n; %#ok<NASGU>
    waitbar(1,f,'Perturbations complete...')
    close(f)
