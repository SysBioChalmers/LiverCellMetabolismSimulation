function compareDistributions(proteomicsData)

meandist = zeros(length(proteomicsData),1);
stddist = zeros(length(proteomicsData),1);
nrofPoints = zeros(length(proteomicsData),1);

for i = 1:length(proteomicsData)
    subplot(3,2,i)
    hold all
    proteomicsFile = importdata(proteomicsData{i});
    protCount = proteomicsFile.data;
    histogram(log10(protCount));
    meandist(i) = mean(log10(protCount));
    stddist(i) = std(log10(protCount));
    nrofPoints(i)=length(protCount);
    histogram(randn(nrofPoints(i),1)*stddist(i) + meandist(i));
end
subplot(3,2,6)
    histogram(randn(mean(nrofPoints),1)*mean(stddist) + mean(meandist));
mean(meandist)
mean(stddist)
mean(nrofPoints)

end

