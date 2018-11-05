function plotResults(t, yAbs, yconc, breakPoints, mu, expData, growthDat, plotMets)
global massPerCell
breakPoints = [0;  breakPoints];
breakPoints(end) = [];
primColor = [0    0.4470    0.7410];
secColor = [0.8500    0.3250    0.0980];
    
fillColors = [199 227 187
              232 208 190
              190 209 223
              233 233 233
            ]/256;

    plotAbsolute = false;
    tMax = max(t);
    

    
    %remove datapoints outside of time range
    growthDat(growthDat(:,1)>tMax,:) = [];
    
    
    %growth
    hold all

    modelY = interp1q(t,yAbs(:,1),breakPoints);
    
    maxVal = max((growthDat(:,2) +  growthDat(:,3)) * massPerCell);
    maxVal = max(maxVal, max(yAbs(:,1)));
    
    plotBreakpoints(t, yAbs(:,1), breakPoints, fillColors)

    plot(t,yAbs(:,1), 'linewidth', 2, 'color', primColor, 'linewidth', 3);
    a1 = gca;
    x = growthDat(:,1);
    y = growthDat(:,2)*massPerCell;
    e = growthDat(:,3)*massPerCell;
    errorbar(x,y,e, 'o', 'color', secColor, 'linewidth', 1.3)
    yMax = round(1.1*(max(y)+max(e)));
    
    if plotAbsolute
        ydata = y(:,2:end);
    else
        ydata = yconc(:,2:end);
    end
    
    
    ylabel('cell count')
    
    set(gca,'FontSize', 15, 'FontWeight', 'bold');
    xlim([0 tMax])
    ylim([0 yMax])
        

    
%     for i = 1:length(breakPoints)
%         txtLabel = sprintf('µ=%2.3f', mu(i));
%         text(breakPoints(i), 0.03 * maxVal, txtLabel, 'fontsize', 15, 'FontWeight', 'bold', 'Rotation',35)
%     end    
    

    xlabel('h')
    ylabel('mg cdw')  
    
    
    box off
    % Create second Y axes on the right.
    a2 = axes('YAxisLocation', 'Right');
    % Hide second plot.
    set(a2, 'color', 'none')
    set(a2, 'XTick', [])
    % Set scala for second Y.
    set(a2, 'YLim', [0 yMax/massPerCell])
    xlabel(' ')
    set(gca,'FontSize', 15, 'FontWeight', 'bold');
    
    addlistener(a2,'MarkedClean',@(varargin)set(a1,'Position',get(a2,'Position')));
    
    
    figure()
    %plot AA
    dict = expData.keys;
    subplotDimentionsX = ceil(sqrt(length(plotMets)));

     subplotDimentionsX = 4;
     subplotDimentionsY = ceil(length(plotMets)/subplotDimentionsX);
%     subplotDimentionsX = 1;
%     subplotDimentionsY = 4;    
 


    lastRow = (subplotDimentionsY-1) * (subplotDimentionsX)+1;
    

    for k = 1:length(plotMets)
        i = find(ismember(dict, plotMets{k}));
        subplot(subplotDimentionsY,subplotDimentionsX,k);        
        hold all        
        data = expData(dict{i});
        modelY = interp1q(t,ydata(:,i),breakPoints);

        maxVal = max(data(2,data(1,:)<tMax)+data(3,data(1,:)<tMax));
        maxVal = max(maxVal, max(ydata(:,i)));
        maxVal = 2^ceil(log2(1.1*maxVal));

        plotBreakpoints(t, ydata(:,i), breakPoints, fillColors)

        plot(t,ydata(:,i), 'linewidth', 3, 'color', primColor);
        errorbar(data(1,:), data(2,:), data(3,:), 'o', 'color', secColor, 'linewidth', 1.3);
        titleStr = strrep(dict{i}, '[s]', '');
        text(tMax*0.1,maxVal*0.85,titleStr, 'FontSize', 15, 'FontWeight', 'bold', 'color', 0.4*[1 1 1])

        if not(lastRow<=k)
            set(gca,'XTick',[]);
        %else
            %xlabel('h')
        end
%         if mod(k,subplotDimentionsX) == 1
%             ylabel('µmol')
%         end
        ylim([0 maxVal])
        xlim([0 tMax])
        set(gca,'FontSize', 15, 'FontWeight', 'bold');
    end
    

%     figure()
%     hold all
%     plot(t, sum(yconc,2), '-')
%     ylim([0 inf])
%     ylabel('osmolarity [µm]')
%     xlabel('time [h]')
    
end


function plotBreakpoints(xdata, ydata, breakPoints, fillColors)
breakPoints = [breakPoints; max(xdata)];

    for i = 2:length(breakPoints)
        %plot([breakPoints(i) breakPoints(i)], [0 modelY(i)], '-','Color', breakColor);
        curData = and((xdata>=breakPoints(i-1)), xdata<=breakPoints(i));
        if sum(curData)>0
            fillColors(i-1,:);
            area(xdata(curData), ydata(curData), 'FaceColor', fillColors(i-1,:), 'EdgeColor', 'none');
        end
    end     
end
