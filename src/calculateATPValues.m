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
atpData = totalATP(model, flux);
glycolysisCost = sumRxns(model,flux, {'HMR_4394'; 'HMR_4379'}, [1;1]);
atpData = atpData - glycolysisCost;


eqnsOut = {'maintainance', 'protein synthesis', 'other sinks', 'GAM', 'secreted proteins'};

GAC = sumRxns(model,flux, {'human_GrowthMaintainance'}, [1]);
aminoacids = sumRxns(model,flux, {'human_proteinPool'}, [4.3]);
GAM = GAC - aminoacids;
prot = sumRxns(model,flux, {'newAlbumin'}, [4.3]);
maint = sumRxns(model,flux, {'human_ATPMaintainance'}, [1]);

other = atpData - GAC - prot - maint;

if other<0
    other = 0;
end

outflux = [ maint;
            aminoacids;
            other
            GAM;
            prot
           ];
end

function [outflux, eqnsOut] = ATPinflux(model, flux)
%       HMR_3957+HMR_4456 glucose
%       HMR_5297  Total
%       HMR_3899 glutamine
%       Other
atpData = totalATP(model, flux);
glycolysisCost = sumRxns(model,flux, {'HMR_4394'; 'HMR_4379'}, [1;1]);
atpData = atpData - glycolysisCost;

%eqnsOut = {'oxidation (PYR)', 'oxidation (GLU)', 'glycolysis', 'other'};
%glycolysis = sumRxns(model,flux, {'HMR_4358'; 'HMR_4368'; 'HMR_4379'; 'HMR_4394'}, [1;1;-1;-1]);
%TCAPyruvate = sumRxns(model,flux, {'HMR_3957', 'HMR_3958'}, [1 1])
%%TCAglutamate = sumRxns(model,flux, {'HMR_3802', 'HMR_3804', 'HMR_4109', 'HMR_3827', 'HMR_3807', 'HMR_4115', 'HMR_4852'}, [-1 -1 -1 -1 -1 1 1]);
%TCAglutamate = sumRxns(model,flux, {'HMR_5297', 'HMR_3957', 'HMR_3958'}, [1 -1 -1]);

%if TCAglutamate<0
%    TCAglutamate = 0; 
% end
% OXPHOS = sumRxns(model,flux, {'HMR_6328'}, [1]);
% 
% TCA1ratio = TCAPyruvate/(TCAPyruvate + TCAglutamate);
% TCA2ratio = 1-TCA1ratio;
% 
% TCA1 = TCA1ratio * OXPHOS;
% TCA2 = TCA2ratio * OXPHOS;

glycolysis = sumRxns(model,flux, {'HMR_4358'; 'HMR_4368'; 'HMR_4379'; 'HMR_4394'}, [1;1;-1;-1]);
OXPHOS = sumRxns(model,flux, {'HMR_6916'}, [1]);
TCA = sumRxns(model,flux, {'HMR_4147', 'HMR_4152'}, [-1 1]);
MITATP = sumRxns(model,flux, {'HMR_6328'}, [1 1]);

OXPHOSratio = OXPHOS/(OXPHOS + TCA);
TCAratio = 1-OXPHOSratio;

OXPHOS = OXPHOSratio * MITATP;
TCA = TCAratio * MITATP;

eqnsOut = {'OXPHOS', 'TCA', 'glycolysis', 'other'};


if glycolysis<0
   glycolysis = 0; 
end

%other = atpData - glycolysis - OXPHOS;

% if other<0
%     other = 0;
% end

outflux = [      
            OXPHOS
            TCA
            glycolysis
            %other
           ];

% outflux = [      
%             TCA1
%             TCA2
%             glycolysis
%             other
%            ];
end

function value = sumRxns(model, flux, rxns, stochiometry)
value = 0;
for i = 1:length(rxns)
    rxn = findIndex(model.rxns, rxns{i});
    contribution = stochiometry(i) * flux(rxn);
    value = value + contribution;
end

end

function atpmin = totalATP(model, flux)
    metNames=modifyMetNames(model);
    cytATP = findIndex(metNames, 'ATP[c]');
    stochiometry = full(model.S(cytATP,:));
    %printATPBalance(model, stochiometry, flux); 
    atpFlux = flux' .* stochiometry;
    atpFlux(atpFlux~=0);
    atpmin = sum(atpFlux(atpFlux>0));
end


function printATPBalance(model, stochiometry, flux)
    includedRxns = flux' .*stochiometry~=0;
    atpFlux = flux(includedRxns)'.*stochiometry(includedRxns);
    rxns = model.rxns(includedRxns);
    eqns = constructEquations(model, includedRxns);
    [atpFlux, indx] = sort(atpFlux);
    eqns = eqns(indx);
    rxns = rxns(indx);

    %print fluxes:
    printFlux = true;
    if printFlux
        for i = 1:length(eqns)
            fprintf('%s\t%s\t%2.5f\n', rxns{i}, eqns{i}, atpFlux(i));
        end
    end
end

