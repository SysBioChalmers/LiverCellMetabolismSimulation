addpath('src')
[num,txt,raw] = xlsread('input/rawdata20190415.xlsx');
raw(:,2) = []; %remove NaN column
raw(4,:) = []; %remove NaN row
raw(3,:) = []; %remove condition row

%Remove SSZ experiment
buis = [NaN cell2mat(raw(1,2:end))];
raw(:,buis>30) = [];

raw = translateMetaboliteNames(raw, 'input/nameConversion.tsv');

%Remove buis
raw(1,:) = [];

%Store data on format:
%AA sampleId time condition value
storeSampleData(raw, 'input/keyGrowth.tsv', 'output/growthData.tsv')
