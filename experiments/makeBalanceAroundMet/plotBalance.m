function plotBalance(eqnsIn, eqnsOut, influx, outflux, metabolite)

    exMap = [67 116 160
             80 137 188
             91 155 213
             151 185 224
             190 209 234]/255;
         
    subplot('Position',[0.1 0.5 0.8 0.4])
    hold on    
    matlabWorkaround = [outflux'; zeros(1,length(outflux))];
    bar(matlabWorkaround,'stack', 'LineStyle', 'none');
    l=legend(eqnsOut,'location', 'ne');
    l.PlotChildren = l.PlotChildren(length(l.PlotChildren):-1:1);
    legend boxoff
    xlim([0 2.5])
    set(gca,'xtick',[])
    
    subplot('Position',[0.1 0.1 0.8 0.4])
    hold on
    colormap(exMap) 
    matlabWorkaround = [-influx'; zeros(1,length(influx))];
    bar(matlabWorkaround,'stack', 'LineStyle', 'none');
    legend(eqnsIn,'location', 'se');
    legend boxoff
    xlim([0 2.5])

   

    plot([0 2],[0,0], 'k')
end