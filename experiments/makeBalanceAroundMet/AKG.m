clf

investigateMetabolite = 'acetyl-CoA';

[eqnsIn, eqnsOut, influx, outflux] = extractReactions(smallModel, smallSolution, investigateMetabolite);

% eqnsIn{1} = 'Glutamate';
% eqnsIn{2} = 'AKG';
% influx(2) = sum(influx(2:6))  + sum(influx(8:9));
% influx(3) = [];
% eqnsIn(3) = [];
% eqnsIn{3} = 'TCA';
% eqnsIn{4} = ''; %Nuclotide phospho transfer
% influx(4) = sum(influx(4:end));
% influx(5:end) = [];
% eqnsIn(5:end) = [];
% eqnsIn
% 
% eqnsOut{1} = 'Biomass';
% eqnsOut{2} = 'Maintenance';
% eqnsOut(3:4) = [];
% outflux(3:4) = [];
% eqnsOut{3} = 'Other';
% outflux(3) = sum(outflux(3:end));
% outflux(4:end) = [];
% eqnsOut(4:end) = [];
% 
% eqnsOut


plotBalance(eqnsIn, eqnsOut, influx, outflux, investigateMetabolite)