function [mets, values] = loadFluxes(folder, celltype, condition)
    fileName = [folder '/' celltype '-' num2str(condition) 'mm-.txt'];
    aaData = importdata(fileName);
    mets = aaData.textdata;
    values = aaData.data;
end

