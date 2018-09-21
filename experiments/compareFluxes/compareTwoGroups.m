function [ output_args ] = compareTwoGroups(dataA, dataB, labels, showOutliers)
serialized = [dataA(:); dataB(:)];
lb = min(serialized);
ub = max(serialized);
lb = lb - 0.05 * abs(lb);
ub = ub + 0.05 * abs(ub);



    %Significant difference
    similarity = zeros(length(dataB),1);
   for i = 1:length(dataB)
       [pdca,gn,gl] = fitdist(dataA(:,i),'Kernel','Kernel','epanechnikov')
       similarity(i) = cdf(pdca,dataB(i),mu,sigma);
   end

    [crap, idx] = sort(median(dataA), 'descend');
    dataA = dataA(:,idx);
    dataB = dataB(idx);
    labels = labels(idx);
    similarity = similarity(idx);
    


   orientation = 'horizontal';    
   
   hold all
   plot([0 0], [0 (length(labels)+1)], 'k-')
   boxplot(dataA, 'Orientation',  orientation,  'Widths', 0.5, 'OutlierSize', 4, 'Symbol', 'o', 'Colors', [17 115 187]/256);


   
    set(findobj(gcf,'LineStyle','--'),'LineStyle','-')

    h = findobj(gca,'Tag','Box');
    for j=1:length(h)
        a=patch(get(h(j),'XData'),get(h(j),'YData'), [17 115 187]/256);
        uistack(a,'bottom');
    end    
    uistack(a,'bottom');
%     c = get(gca, 'Children');
%     hleg1 = legend(c(1:2), groups);
    
    
    if not(showOutliers)
        h=findobj(gca,'tag','Outliers');
        delete(h)
    end
    
    if strcmp(orientation, 'horizontal')
        set(gca,'ytick', 1:length(labels))
        set(gca,'yticklabel', labels)      
    else
        set(gca,'xtick', centeredPosition)
        set(gca,'xticklabel', labels)
        set(gca, 'XTickLabelRotation', 45);
    end    

    plot(dataB, 1:length(dataB), 'x', 'linewidth', 2, 'color', [216 85 39]/256);
       
    lines = findobj(gcf, 'type', 'line', 'Tag', 'Median');
    set(lines, 'Color', 'w');

    xlim([lb ub])
    
    similarity = abs(similarity);
   for i = 1:length(dataB)
       if similarity(i)>0.5
           similarity(i) = 1 - similarity(i);
       end
       %convert to paired test
       similarity(i) = similarity(i)*2;
       
       if similarity(i)>0.05
           text(ub*0.85, i, sprintf('%1.2f',similarity(i)), 'color', 0.6*[1 1 1]);
       elseif similarity(i)>0.01
           text(ub*0.85, i, sprintf('%1.2f',similarity(i)), 'color', [216 85 39]/256);
       else
           text(ub*0.85, i, sprintf('%1.2e',similarity(i)), 'color', [216 85 39]/256);
       end
   end

   xlabel('Normalized Flux')
    
end

