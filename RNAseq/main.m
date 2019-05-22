A = importdata('transcript_rna_celline.tsv/transcript_rna_celline.tsv');
categories = A.textdata(1,:);
geneNames = A.textdata(2:end,1);
transcriptNames = A.textdata(2:end,2);

match = ismember(categories, {'Hep G2.C13', 'Hep G2.C14'});
match(1:2) = [];
geneExpression = A.data(:,match);


fileID = fopen('seq.txt','w');
for i = 1:length(geneExpression)
    fprintf(fileID,'%s',geneNames{i});
    fprintf(fileID,'\t%s',transcriptNames{i});
    
    for j = 1:size(geneExpression,2)
        fprintf(fileID,'\t%f',geneExpression(i,j));
    end
    fprintf(fileID,'\r\n');
end

fclose(fileID);

%%
clear
A = importdata('transcript_rna_tissue.tsv/transcript_rna_tissue.tsv');
categories = A.textdata(1,:);
geneNames = A.textdata(2:end,1);
transcriptNames = A.textdata(2:end,2);

match = ismember(categories, {'liver.V108', 'liver.V110', 'liver.V111', 'liver.V348', 'liver.V349', 'liver.V350', 'liver.V351', 'liver.V358', 'liver.V362', 'liver.V363'});
match(1:2) = [];
geneExpression = A.data(:,match);

fileID = fopen('tisue.txt','w');
for i = 1:length(geneExpression)
    fprintf(fileID,'%s',geneNames{i});
    fprintf(fileID,'\t%s',transcriptNames{i});
    for j = 1:size(geneExpression,2)
        fprintf(fileID,'\t%f',geneExpression(i,j));
    end
    fprintf(fileID,'\r\n');
end

fclose(fileID);

%%
