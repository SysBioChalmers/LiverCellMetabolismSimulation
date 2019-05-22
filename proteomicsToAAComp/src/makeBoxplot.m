function makeBoxplot(samples, labels, showOutliers)

    boxplot(samples, 'Labels', labels, 'Whisker', 100, 'Widths', 0.8, 'Orientation', 'horizontal', 'Colors', [17 115 187]/256, 'Symbol', 'o');
    set(findobj(gcf,'LineStyle','--'),'LineStyle','-')
    h = findobj(gca,'Tag','Box');
    for j=1:length(h)
        a=patch(get(h(j),'XData'),get(h(j),'YData'), [17 115 187]/256);
        uistack(a,'bottom');
    end       
    uistack(a,'bottom');
    
    if not(showOutliers)
        h=findobj(gca,'tag','Outliers');
        delete(h)
    end    
    
    lines = findobj(gcf, 'type', 'line', 'Tag', 'Median');
    set(lines, 'Color', 'w');
    
    xlabel('mol ratio')
end

