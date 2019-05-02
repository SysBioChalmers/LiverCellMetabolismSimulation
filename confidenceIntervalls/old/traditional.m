clc
close all
addpath('../src')
addpath('src')
color = [93 155 211
         215 86 40
         238 178 32]/256;

conditions = {'glc'; 'glc-'};
growthData = importdata('growthData.txt');
data = growthData.data;
data(data(:,2)>0,2) = 1;

t = data(:,1);
x = log(data(:,3));
g = data(:,2);

dataTable = table(t, x, g);
fitresult1 = fitlm(t(g==0),x(g==0));
fitresult2 = fitlm(t(or(g==1,t==0)),x(or(g==1,t==0)));

inspectFits(t(g==0),x(g==0), fitresult1, color(1,:))
inspectFits(t(or(g==1,t==0)),x(or(g==1,t==0)), fitresult2, color(2,:))
figure()

t0 = [23.2500 30.2500 47.5000]';
mu = [fitresult1.Coefficients.Estimate(2) fitresult2.Coefficients.Estimate(2)];
muE = [fitresult1.Coefficients.SE(2) fitresult2.Coefficients.SE(2)];
X = exp([predict(fitresult1, t0), predict(fitresult2, t0)]);
DX = (X-X(1,:));
DXpMu = DX./mu;
DXpMu = [DXpMu DXpMu(:,2)];
plotFluxes(mu,muE,{'mu'},flip(conditions),color)

fitresult = fitlm(dataTable, 'x~t:(1+g)')
p = fitresult.Coefficients.pValue(3);
text(1.2, 0.038, sprintf('p=%2.3f', p))


%%
figure()
hold all
metData = importdata('metData.txt');
data = metData.data;
mets = metData.textdata(2:end,1);
volData = importdata('volume.txt');

metabolites = unique(mets);
wells{1} = [4 7 10 13 16];
wells{2} = 1+wells{1};
wells{3} = 1+wells{2};

linearModels = cell(length(metabolites),length(wells));

for i = 1:length(metabolites)
    subData = data(ismember(mets, metabolites{i}),:);
    S = getSestimate(subData,volData);
    subplot(2,2,i)
    for j = 1:length(wells)
        wellIndx = ismember(subData(:,2),wells{j});
        curX = getXestimate(subData(wellIndx,1),t0,DXpMu(:,j));
        curS = S(wellIndx);
        linearModels{i,j} = fitlm(curX,curS);
        inspectFits(curX, curS, linearModels{i,j}, color(j,:));
    end
end


figure();
hold all
conditions = {'0mM', '6 mM', '22mM'};
massPerCell = 1000 * 426.8 * 10^-12; %mg per cell
[flux, error] = extractAllFluxes(linearModels, massPerCell);


plotFluxes(flux,error,metabolites,conditions,color)

