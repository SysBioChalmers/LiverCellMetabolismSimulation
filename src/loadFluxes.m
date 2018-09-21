function [mets, values] = loadFluxes(folder, inputFile)
    aaData = importdata([folder '/' inputFile]);
    mets = aaData.textdata;
    values = aaData.data;
end

