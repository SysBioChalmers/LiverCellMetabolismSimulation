%plotSmallModel(smallModel, metThreshold, rxnThreshold)
metThreshold = 7;

%[cMatrix, labels] = generateCmatrix(smallModel, metThreshold, rxnThreshold, false);
%[cMatrix, labels] = generateFilteredMatrix(smallModel, metThreshold, rxnThreshold);

model = smallModel;

for i = 1:length(model.rxns)
    subRxns = split(model.rxns{i},'+');
    model.rxns{i} = subRxns{1};
end

whiteList = {'homocysteine[c]'};
blackList = {};
biomass = {'HMR_6916', 'HMR_4396', 'human_DNAPool','human_RNAPool','PhosphatidylPool', 'lipidPool'};


[cMatrix, labels, rxnStart, exchangeMets, bioMets] = generateBiPartite(model, smallSolution, metThreshold, biomass, blackList, whiteList);

graph_to_dot('test.dot', cMatrix,  labels, rxnStart, exchangeMets, bioMets)


delete('test.dot.pdf')
commandStr = 'python displayNetwork/makePdf.py';
 [status, commandOut] = system(commandStr);
 commandOut
 if status==0
     fprintf('pass');
 else
     fprintf('fail');
 end

% test = biograph(cMatrix, labels);
% h = view(test);
% 
% for i = 1:length(h.Nodes)
%     set(h.Nodes(i),'FontSize',10)
%     set(h.Nodes(i),'Color',[0.2 0.7 1])    
%     set(h.Nodes(i),'LineColor',[0.5 0.5 0.5]) 
% end
% set(h, 'ShowArrows', 'off')
% set(h, 'Scale', 0.3)
% set(h, 'Scale', 0.5)


