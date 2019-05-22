clc
close all
addpath('src')
celltype = 'hepg2';
condition = {'0', '6', '22'};

nrOfMets = 2;
k=0;

secColor  = [216 85 39]/256;

resultsIntersect = zeros(length(condition),2);
resultsMu = zeros(length(condition),2);
figure('DefaultAxesFontSize',12)
figure('DefaultAxesFontWeight','bold')

rawData = IO('additionalSamples/output/growthData.tsv');

showAdditionalData = false;

N = length(condition);
M = 1 + nrOfMets;

glutamineRate =  0.0023;

tvals = linspace(0, 100,1000);
glutamineDegradation = zeros(length(condition),length(tvals));

index = reshape(1:(N*M), N, M)';
for i = 1:length(condition)
    [expData, volumePoints, growthData] = loadExpdata('data', celltype, condition{i});
    growthData = meanAndError(growthData(:,1),growthData(:,2));
    
    growthData(end,:) = [];
    X = growthData';
    
    k = k+1;
    subplot(M,N,index(k))
    hold all
    
    
    Xvals = exp(interp1q(X(1,:)', log(X(2,:))',tvals'));
    
    area(tvals, Xvals, 'FaceColor', [0.9 0.9 0.9], 'EdgeColor', 'none');
    plot(tvals, Xvals, 'k', 'linewidth', 2)
    errorbar(X(1,:), X(2,:), X(3,:), 'kx')
    
%     km = polyfit(X(1,[1 2]), log(X(2,[1 2])),1);
%     mu = km(1);
%     text(10, max(Xvals)*0.6, sprintf('mu=%2.3f', mu));
%     plot(tvals, exp(km(2))*exp(mu*tvals), 'r--');
%     ylim([0 max(Xvals)*1.2])
    ylabel('cells')
    xlabel('time [h]')
    xlim([0 50])
    

    k = k+1;
    subplot(M,N,index(k))
    hold all
    Sdata = expData('pyruvate[s]');
    Sdata = meanAndError(Sdata(1,:)',Sdata(2,:)')';    
    
    Sdata(:,end) = [];
    Sdata(2:3,:) = Sdata(2:3,:)/1000; %to mM
    ymax= max(Sdata(2,:))*1.2;

    dataPoints = 3:5;
    T = Sdata(1,dataPoints)';
    X = interp1q(tvals', Xvals, T);
    S = Sdata(2,dataPoints)';
    svals1 = interpolateMetabolite(T, X, S, tvals);
    

    dataPoints = 1:2;
    T = Sdata(1,dataPoints)';
    X = interp1q(tvals', Xvals, T);
    S = Sdata(2,dataPoints)';
    svals2 = interpolateMetabolite(T, X, S, tvals);
    
    [x, y] = min(abs(svals2-svals1));
    %plot(tvals(y), svals1(y), 'ko');
    crossT = tvals(y);

    areaVals = min([svals1'; svals2']);
    area(tvals, areaVals, 'FaceColor', [0.9 0.9 0.9], 'EdgeColor', 'none');    
    
    plot(tvals, svals1, 'linewidth', 2)
    plot(tvals, svals2, 'linewidth', 2)
    errorbar(Sdata(1,:), Sdata(2,:), Sdata(3,:), 'kx')    
    plot(crossT * [1 1], [0 ymax], 'k--');
    text(crossT*1.1, 0.1*ymax, sprintf('t=%2.1f', crossT));
    ylabel('pyruvate [mM]')
    xlabel('time [h]')
    ylim([0 1.5])
    xlim([0 50])
    
    k = k+1;
    subplot(M,N,index(k))
    hold all
    Sdata = expData('glutamine[s]');
    Sdata(:,end) = [];
    Sdata(2:3,:) = Sdata(2:3,:)/1000; %to mM
    

    dataPoints = 2:3;
    T = Sdata(1,dataPoints)';
    X = interp1q(tvals', Xvals, T);
    S = Sdata(2,dataPoints)';
    svals1 = interpolateMetabolite(T, X, S, tvals);

    crossS = interp1q(tvals, svals1, crossT);
    T = [0; crossT];
    X = interp1q(tvals', Xvals, T);
    S = [Sdata(2,1); crossS];
    svals2 = interpolateMetabolite(T, X, S, tvals);
    
    
    areaVals = max([svals1'; svals2']);
    area(tvals, areaVals, 'FaceColor', [0.9 0.9 0.9], 'EdgeColor', 'none');
    
    plot(tvals, svals1, 'linewidth', 2)
    plot(tvals, svals2, 'linewidth', 2)
    errorbar(Sdata(1,:), Sdata(2,:), Sdata(3,:), 'kx')    
    %plot(crossT, crossS, 'ko', 'linewidth', 2)
    plot([0 max(tvals)], crossS * [1 1], 'k--');
    text(crossT+10, crossS*1.2, sprintf('S=%2.1f mM', crossS));
    ylabel('glutamine [mM]')
    xlabel('time [h]')
    ylim([0 inf])
    xlim([0 50])
    
    glutamineDegradation(i,:) = 1000*estimateGlutamineDegradation(areaVals, glutamineRate)./Xvals;
    
    if showAdditionalData
        curRaw = rawData;
        curRaw(not(ismember(curRaw(:,1),'Glutamine')),:) = [];
        curRaw = cell2nummat(curRaw(:,[3 4 6]));
        curRaw = curRaw(curRaw(:,2) == str2num(condition{i}),:);
        translucentScatter(curRaw(:,1), curRaw(:,3)/1000, secColor, 0.3, 0.8)
    end
    
    resultsIntersect(i,:) = [crossT crossS]; 
    
%      %serine
%      k = k+1;
%      subplot(M,N,index(k))
%      hold all
%      Sdata = expData('serine[s]');
%      Sdata(:,end) = [];
%      Sdata(2:3,:) = Sdata(2:3,:)/1000; %to mM
%     
%      dataPoints = 2:3;
%      T = Sdata(1,dataPoints)';
%      X = interp1q(tvals', Xvals, T);
%      S = Sdata(2,dataPoints)';
%      svals1 = interpolateMetabolite(T, X, S, tvals);
% 
%      crossS = interp1q(tvals, svals1, crossT);
%      T = [0; crossT];
%      X = interp1q(tvals', Xvals, T);
%      S = [Sdata(2,1); crossS];
%      svals2 = interpolateMetabolite(T, X, S, tvals);
% 
%      areaVals = max([svals1'; svals2']);
%      area(tvals, areaVals, 'FaceColor', [0.9 0.9 0.9], 'EdgeColor', 'none');
%      
%      plot(tvals, svals1, 'linewidth', 2)
%      plot(tvals, svals2, 'linewidth', 2)
%      errorbar(Sdata(1,:), Sdata(2,:), Sdata(3,:), 'kx')    
% 
%      ylabel('serine [mM]')
%      xlabel('time [h]')
%      ylim([0 inf])
%      xlim([0 50])
%      
%      curRaw = rawData;
%      curRaw(not(ismember(curRaw(:,1),'Serine')),:) = [];
%      curRaw = cell2nummat(curRaw(:,[3 4 6]));
%      curRaw = curRaw(curRaw(:,2) == str2num(condition{i}),:);
%      translucentScatter(curRaw(:,1), curRaw(:,3)/1000, secColor, 0.3, 0.8)
% 
%      %leucine
%      k = k+1;
%      subplot(M,N,index(k))
%      hold all
%      Sdata = expData('leucine[s]');
%      Sdata(:,end) = [];
%      Sdata(2:3,:) = Sdata(2:3,:)/1000; %to mM
%     
%      dataPoints = 2:3;
%      T = Sdata(1,dataPoints)';
%      X = interp1q(tvals', Xvals, T);
%      S = Sdata(2,dataPoints)';
%      svals1 = interpolateMetabolite(T, X, S, tvals);
% 
%      crossS = interp1q(tvals, svals1, crossT);
%      T = [0; crossT];
%      X = interp1q(tvals', Xvals, T);
%      S = [Sdata(2,1); crossS];
%      svals2 = interpolateMetabolite(T, X, S, tvals);
% 
%      areaVals = max([svals1'; svals2']);
%      area(tvals, areaVals, 'FaceColor', [0.9 0.9 0.9], 'EdgeColor', 'none');
%      
%      plot(tvals, svals1, 'linewidth', 2)
%      plot(tvals, svals2, 'linewidth', 2)
%      errorbar(Sdata(1,:), Sdata(2,:), Sdata(3,:), 'kx')    
% 
%      ylabel('leucine [mM]')
%      xlabel('time [h]')
%      ylim([0 inf])
%      xlim([0 50])
%      
%      curRaw = rawData;
%      curRaw(not(ismember(curRaw(:,1),'Leucine')),:) = [];
%      curRaw = cell2nummat(curRaw(:,[3 4 6]));
%      curRaw = curRaw(curRaw(:,2) == str2num(condition{i}),:);
%      translucentScatter(curRaw(:,1), curRaw(:,3)/1000, secColor, 0.3, 0.8)     
%      
%      %isoleucine
%      k = k+1;
%      subplot(M,N,index(k))
%      hold all
%      Sdata = expData('isoleucine[s]');
%      Sdata(:,end) = [];
%      Sdata(2:3,:) = Sdata(2:3,:)/1000; %to mM
%     
%      dataPoints = 2:3;
%      T = Sdata(1,dataPoints)';
%      X = interp1q(tvals', Xvals, T);
%      S = Sdata(2,dataPoints)';
%      svals1 = interpolateMetabolite(T, X, S, tvals);
% 
%      crossS = interp1q(tvals, svals1, crossT);
%      T = [0; crossT];
%      X = interp1q(tvals', Xvals, T);
%      S = [Sdata(2,1); crossS];
%      svals2 = interpolateMetabolite(T, X, S, tvals);
% 
%      areaVals = max([svals1'; svals2']);
%      area(tvals, areaVals, 'FaceColor', [0.9 0.9 0.9], 'EdgeColor', 'none');
%      
%      plot(tvals, svals1, 'linewidth', 2)
%      plot(tvals, svals2, 'linewidth', 2)
%      errorbar(Sdata(1,:), Sdata(2,:), Sdata(3,:), 'kx')    
% 
%      ylabel('isoleucine [mM]')
%      xlabel('time [h]')
%      ylim([0 inf])
%      xlim([0 50])
%      
%      curRaw = rawData;
%      curRaw(not(ismember(curRaw(:,1),'Iso-leucine')),:) = [];
%      curRaw = cell2nummat(curRaw(:,[3 4 6]));
%      curRaw = curRaw(curRaw(:,2) == str2num(condition{i}),:);
%      translucentScatter(curRaw(:,1), curRaw(:,3)/1000, secColor, 0.3, 0.8)     
%           
     
end


%%
figure()
hold on
massPerCell = 1000 * 426.8 * 10^-12; %mg per cell
glnflux = glutamineDegradation./massPerCell;
fill([23.2500   47.5000 47.5000 23.2500], [-20 -20 0   0], [0.8 0.8 0.8], 'edgecolor', 'none')
plot(tvals, glnflux)
filter = and(tvals>23.2500, tvals<47.5000);
average = mean(glnflux(:,filter),2);
for i = 1:length(average)
    plot([0 60], average(i) * [1 1], 'k--')
end
xlim([0 60])
ylim([-20 0])
xlabel('h')
ylabel('glutamine degradation')