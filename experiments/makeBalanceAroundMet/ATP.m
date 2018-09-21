clf

tmpModel = smallModel;

investigateMetabolite = 'ATP';

[eqnsIn, eqnsOut, influx, outflux] = extractReactions(tmpModel, smallSolution, investigateMetabolite);
eqnsIn = constructEquations(smallModel, eqnsIn);
eqnsOut = constructEquations(smallModel, eqnsOut);


eqnsIn{1} = 'ATP synthase';
eqnsIn{2} = 'Glycolysis';
influx(2) = influx(2) + influx(3) - outflux(3) - outflux(4);
influx(3) = [];
eqnsIn(3) = [];
eqnsIn{3} = 'TCA';
eqnsIn{4} = ''; %Nuclotide phospho transfer
influx(4) = sum(influx(4:end));
influx(5:end) = [];
eqnsIn(5:end) = [];
eqnsIn

eqnsOut{1} = 'Biomass';
eqnsOut{2} = 'Maintenance';
eqnsOut(3:4) = [];
outflux(3:4) = [];
eqnsOut{3} = 'Other';
outflux(3) = sum(outflux(3:end));
outflux(4:end) = [];
eqnsOut(4:end) = [];

eqnsOut

sum(outflux)-sum(influx)

plotBalance(eqnsIn, eqnsOut, influx, outflux, investigateMetabolite)


