function p = plot_with_shade(x, Y, prc, alphaa, color)
    p = plot(x, mean(Y, 2), 'linewidth', 3, 'color', color);
    x2 = [x(:)', flip(x(:)')];
    fill(x2, [prctile(Y, prc, 2)', flip(prctile(Y, 100-prc, 2))'], color, 'LineStyle', 'none');
    alpha(alphaa)
end