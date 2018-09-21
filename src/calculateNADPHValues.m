function [outflux, eqnsOut] = calculateATPValues(model, flux, direction)

if strcmp(direction, 'in')
    [outflux, eqnsOut] = ATPinflux(model, flux);
else
    [outflux, eqnsOut] = ATPoutflux(model, flux);
end

end

function [outflux, eqnsOut] = ATPoutflux(model, flux)
%       newAlbumin (protein synthesis)
%       HumanGrowth (biomass)
%       human_ATPMaintainance (maintainance)
%       Other
atpData = totalNADPH(model, flux)

eqnsOut = {'Biomass', 'Protein synthesis', 'Maintainance', 'other'};

GAM = sumRxns(model,flux, {'HumanGrowth'}, [150]);
prot = sumRxns(model,flux, {'newAlbumin'}, [4]);
maint = sumRxns(model,flux, {'human_ATPMaintainance'}, [1]);

other = atpData - GAM - prot - maint;

if other<0
    other = 0;
end

outflux = [ GAM;
            prot
            maint
            other
           ];
end

function [outflux, eqnsOut] = ATPinflux(model, flux)
%       HMR_3957+HMR_4456 glucose
%       HMR_5297  Total
%       HMR_3899 glutamine
%       Other
atpData = totalNADPH(model, flux);

eqnsOut = {'Glycolysis', 'Pyruvate->TCA', 'Glutamine->TCA', 'Other'};

glycolysis = sumRxns(model,flux, {'HMR_4358'; 'HMR_4368'; 'HMR_4379'; 'HMR_4394'}, [1;1;-1;-1]);
TCA1 = sumRxns(model,flux, {'HMR_3957'}, [1]);
TCA2 = sumRxns(model,flux, {'HMR_5297'; 'HMR_3957'}, [1;-1]);

OXPHOS = sumRxns(model,flux, {'HMR_6916'; 'HMR_4152'}, [1;1]);

TCA1 = TCA1/(TCA1 + TCA2);
TCA2 = 1-TCA1;

TCA1 = TCA1 * OXPHOS;
TCA2 = TCA2 * OXPHOS;

if glycolysis<0
   glycolysis = 0; 
end

other = atpData - glycolysis - TCA1 - TCA2;

if other<0
    other = 0;
end


outflux = [glycolysis;
            TCA1
            TCA2
            other
           ];
end

function value = sumRxns(model, flux, rxns, stochiometry)
value = 0;
for i = 1:length(rxns)
    rxn = findIndex(model.rxns, rxns{i});
    value = value + stochiometry(i) * flux(rxn);     
end

end

function NADPHmin = totalNADPH(model, flux)
    metNames=modifyMetNames(model);
    cytNADPH = findIndex(metNames, 'NADPH[c]');
    stochiometry = model.S(cytNADPH,:);
    NADPHFlux = flux' .* stochiometry;
    NADPHmin = sum(NADPHFlux(NADPHFlux>0));
end