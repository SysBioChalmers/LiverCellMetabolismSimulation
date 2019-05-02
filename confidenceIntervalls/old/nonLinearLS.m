clc
close all
addpath('src')
color = [93 155 211
         215 86 40
         238 178 32]/256;


volData = importdata('volume.txt');
volData(:,1) = volData(:,1)-volData(1,1);
volData(1,2) = volData(2,2);


%condition = '22mM';
%[dataX, dataY, metlabels] = makeDataStructure(condition);
condition = '22';
[dataX, dataY, metlabels] = makeDataStructureNew(condition);

tvals = unique(dataX(:,1))';
x0 = estimateInitialX(dataX, dataY, tvals);
[x,resnorm,residual,exitflag,output,lambda,jacobian] = lsqcurvefit(@(x,y) fitFunction(x,y,tvals,volData(:,2)),x0,dataX,dataY);
conf = nlparci(x,residual,'jacobian',jacobian);
x
figure()
hold all
predGrowth = x(1,2);
confGrowth = conf(size(x,1)+1,:);
bar(1,predGrowth)
errorbar(1, predGrowth, confGrowth(1)-predGrowth, confGrowth(2)-predGrowth,'k.');
ylim([0 0.04])

figure()
hold all
massPerCell = 426.8 * 10^-12; %g per cell
massPerMCell = massPerCell*10^6;
predFlux = x(2:end,2)/massPerMCell;
confFlux = conf((size(x,1)+2):end,:)/massPerMCell;
xLocation = 1:length(predFlux);
bar(xLocation,predFlux)
errorbar(xLocation, predFlux, confFlux(:,1)-predFlux, confFlux(:,2)-predFlux,'k.');
ylim([-200 200])
xticks(xLocation)
xticklabels(metlabels(2:end));
xtickangle(90)

figure()
[x,resnorm,residual,exitflag,output,lambda,jacobian] = lsqcurvefit(@(x,y) fitFunction(x,y,tvals,volData(:,2)),x0,dataX,dataY);

[x,R,J,CovB,MSE,ErrorModelInfo] = nlinfit(___) 
[x,residual,J,CovB,MSE] = nlinfit(x,y,modelfun,beta0,opts);
conf = nlparci(x,residual,'jacobian',jacobian);
[ypred,delta] = nlpredci(modelfun,xrange,x,residual,'Covar',CovB,'MSE',MSE,'SimOpt','on');
lower = ypred - delta;
upper = ypred + delta;
for i = 1:length(xLocation)
    subplot(6,4,i)
    
end