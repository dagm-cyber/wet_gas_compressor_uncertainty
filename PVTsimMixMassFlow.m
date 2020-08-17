function newComp = PVTsimMixMassFlow(comps,massflows,Ncomp,Fluid)

%Allocate memory
Mw = zeros(1,Ncomp);

%Collect Mw(i) from Fluid
for i = 1:Ncomp
    Mw(i) = Fluid.Composition.ComponentParams(i).Mw;
end

%Calculate new mixed composition based on massflows
newComp = zeros(1,Ncomp);
for i=1:length(comps(1,:))
    curcomp = comps(:,i); %Collect i-th composition vector
    curcomp = MolToWt(curcomp,Mw); % convert to massPercent
        for j=1:Ncomp
            newComp(j) = newComp(j)+curcomp(j)*massflows(i); %Calculate wtP of component j based on massflow(i) of compnent(j)      
        end       
end
newComp = norm(newComp);
newComp = WtToMol(newComp,Mw);
newComp = 100*norm(newComp);
end

%% Convert from weigth to Mole fration    
    function [x,xP] = WtToMol(Wt,Mw)
    %Give same result both for Wt given in % and fraction
    
    if iscolumn(Wt)
        Wt=Wt';
    end
    if iscolumn(Mw)
        Mw=Mw';
    end
    
        WtDivMw = Wt./Mw;
        x = WtDivMw./sum(WtDivMw);
        xP = 100*x;
    end
    
%% Convert from mole to weigth fration
    function [Wt,WtP] = MolToWt(x,Mw)
    %Give same result both for x given in % and fraction
    
    if iscolumn(x)
        x=x';
    end
    if iscolumn(Mw)
        Mw=Mw';
    end
      
        XiMwi = x.*Mw;
        Wt = XiMwi./sum(XiMwi);
        WtP = 100.*Wt;
        
    end
 
 %% Normalization function
    function X = norm(X)
    X=X./sum(X);
    end


   