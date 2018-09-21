clf

investigateMetabolite = 'glutamate';

[inRxns, outRxns, influx, outflux] = extractReactions(smallModel, smallSolution, investigateMetabolite);
eqnsIn = constructEquations(model, inRxns);
eqnsOut = constructEquations(model, outRxns);

eqnsIn
eqnsOut
eqnsIn{1} = 'glutamine';
eqnsIn{2} = 'BCAA';
influx(2) = sum(influx(2:4));
eqnsIn(3:4) = [];
eqnsIn{3} = 'cysteine';
eqnsIn{4} = 'asparagine';
eqnsIn{5} = 'Nucleotides';
% influx(5) = sum(influx(5;));
% eqnsIn(5) = [];

% influx(3) = [];
% eqnsIn(3) = [];
% eqnsIn{3} = 'TCA';
% eqnsIn{4} = ''; %Nuclotide phospho transfer
% influx(4) = sum(influx(4:end));
% influx(5:end) = [];
% eqnsIn(5:end) = [];
% eqnsIn
% 
eqnsOut{1} = 'AKG';
eqnsOut{2} = 'production';
eqnsOut{3} = 'aspartate';
eqnsOut{4} = 'proline';
eqnsOut{5} = 'serine';
eqnsOut{6} = 'biomass';
% eqnsOut(3:4) = [];
% outflux(3:4) = [];
% eqnsOut{3} = 'Other';
% outflux(3) = sum(outflux(3:end));
% outflux(4:end) = [];
% eqnsOut(4:end) = [];
% 
% eqnsOut


plotBalance(eqnsIn, eqnsOut, influx, outflux, investigateMetabolite)