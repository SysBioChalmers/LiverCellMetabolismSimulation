function inspectFits(Xdata, Ydata, a, color)
hold all
    xvals = linspace(min(Xdata), max(Xdata));
    [yvals, aCi] = predict(a, xvals');
    scatter(Xdata, Ydata, 'filled', 'markerfacecolor', color, 'markeredgecolor', 'none')
    fillArea(xvals, aCi, color)
    plot(xvals, yvals, 'color', color, 'linewidth', 2)
end

