function [data1, transcripts1] = extractGeneData(A, proteinCodingTranscripts, geneName)
    match1 = ismember(A.textdata(:,1), geneName);

    data1 = A.data(match1,:);

    transcripts1 = A.textdata(match1,2);

    transcriptMatch1 = ismember(transcripts1, proteinCodingTranscripts);

    data1 = data1(transcriptMatch1,:);
    
    transcripts1 = transcripts1(transcriptMatch1);  
end

