function saveFlux(condition, metlabels, predFlux, confFlux)
    metaboliteMap = IO('../data/metaboliteMap.tsv');
    metTranslation = containers.Map(metaboliteMap(:,1),metaboliteMap(:,2));

    filename = ['output/hepg2-' condition '.tsv'];

    fileID = fopen(filename,'w');
    fprintf(fileID,'Met\tFlux\tLb\tUb');

    for i = 1:length(metlabels)
        curMet = metlabels{i};
        curMet = metTranslation(curMet);
        fprintf(fileID,'\n');
        f = predFlux(i);
        fe1 = confFlux(i,1);
        fe2 = confFlux(i,2);
        fprintf(fileID,'%s\t%2.5f\t%2.5f\t%2.5f', curMet, f, fe1, fe2);
    end

    fclose(fileID);
end

