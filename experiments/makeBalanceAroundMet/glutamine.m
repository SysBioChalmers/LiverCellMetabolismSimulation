clf

investigateMetabolite = 'glutamine';

[eqnsIn, eqnsOut, influx, outflux] = extractReactions(smallModel, smallSolution, investigateMetabolite);

eqnsOut{1} = 'glutamate';
eqnsOut{2} = 'biomass';
eqnsOut{3} = 'asparagine'; %transamination
eqnsOut{4} = 'nucleotides'; %synthesis
outflux(4) = sum(outflux(4:end));
outflux(5:end) = [];
eqnsOut(5:end) = [];

eqnsOut
eqnsIn{1} = 'uptake';

plotBalance(eqnsIn, eqnsOut, influx, outflux, investigateMetabolite)