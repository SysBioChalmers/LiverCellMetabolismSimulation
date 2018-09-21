clc
load('workspace.mat') %load model and FBA solution.


%printExchangeFluxes(model, solution.x);
size(model.S)
[smallModel, smallSolution] = simplifyFluxResult(model, solution.x, true);
size(smallModel.S)

%Can be repeated for some additional compression...
[smallModel, smallSolution] = simplifyFluxResult(smallModel, smallSolution, false);
size(smallModel.S)
[smallModel, smallSolution] = simplifyFluxResult(smallModel, smallSolution, false);
size(smallModel.S)
[smallModel, smallSolution] = simplifyFluxResult(smallModel, smallSolution, false);
size(smallModel.S)

%printExchangeFluxes(smallModel, smallSolution);

plotHierarchicalModel(smallModel, smallSolution, {'alanine[c]', 'glycine[c]'});
xlim([1 15])

validateBalances = zeros(length(smallModel.mets),1);
for i = 1:length(smallModel.mets)
    a=smallSolution' .* smallModel.S(i,:);
    validateBalances(i) = sum(a);
end

sum(abs(validateBalances))