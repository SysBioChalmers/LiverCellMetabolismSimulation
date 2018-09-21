

investigateMetabolite = 'glutamine';

[eqnsIn, eqnsOut, influx, outflux] = extractReactions(smallModel, smallSolution, investigateMetabolite);

eqnsOut{7} = 'Glutamate';

plotBalance(eqnsIn, eqnsOut, influx, outflux)