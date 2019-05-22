addpath('src')
[num,txt,raw] = xlsread('input/rawdata20190415.xlsx');
raw(:,2) = []; %remove NaN column
raw(4,:) = []; %remove NaN row
raw(3,:) = []; %remove condition row

%Remove timeseries experiment
buis = [NaN cell2mat(raw(1,2:end))];
raw(:,buis<=30) = [];

%Uppdate sample ids (based on email from Jurgen 2019-04-16)
%The last 16 samples are 1-16, measured in this order, with 16 as the final
%sample I have checked this with Albert. You can see in the word file which
%sample is what. 13 and 16 are the t=0 samples as you can also see in the gln.
for i = 2:size(raw,2)
    raw{2,i} = [i-1];
end

%Translate amino acid names to english
raw = translateMetaboliteNames(raw, 'input/nameConversion.tsv');

%Remove buis
raw(1,:) = [];

%Store data on format:
%AA sampleId time condition value
storeSampleData(raw, 'input/keySSZ.tsv', 'output/SSZdata.tsv')

