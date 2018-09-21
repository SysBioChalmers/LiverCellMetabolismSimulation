clf
addpath('../')

model.subSystems{8128} = 'Maintenance';
model.subSystems{8131} = 'Growth';

investigateMetabolite = 'glutamate';

[inRxns, outRxns, influx, outflux] = extractReactions(model, solution.x, investigateMetabolite);
[eqnsIn, eqnsOut, influx, outflux] = mergeBySubsystem(model, inRxns, outRxns, influx, outflux);

sum(outflux)-sum(influx)

plotBalance(eqnsIn, eqnsOut, influx, outflux, investigateMetabolite)




