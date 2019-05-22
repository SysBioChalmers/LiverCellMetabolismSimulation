function [data1, data2, transcripts1, transcripts2] = getGeneData(A, B, proteinCodingTranscripts, geneName)
    [data1, transcripts1] = extractGeneData(A, proteinCodingTranscripts, geneName);
    [data2, transcripts2] = extractGeneData(B, proteinCodingTranscripts, geneName);   
end

