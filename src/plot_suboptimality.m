%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Semi-Disrete Optimal Transport: Hardness, Regularization and Numerical
% Solution
% Authors: Bahar Taskesen, Soroosh Shafieezadeh-Abadeh, Daniel Kuhn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the suboptimality gap of the Algorithm~1 for the original, entopic 
% regularized and chi2 divergence regularized optimal transport problems.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('non-reg_reg_sub_opt.mat')
non_reg_sub_opt = sub_opt;
load('entropic_reg_sub_opt.mat')
entropic_reg_sub_opt = sub_opt;
load('entropic2_reg_sub_opt.mat')
entropic2_reg_sub_opt = sub_opt;
load('chi2_reg_sub_opt.mat')
tikhonov_reg_sub_opt = sub_opt;
load('hyperbolic_reg_sub_opt.mat')
hyperbolic_reg_sub_opt = sub_opt;
run params.m
% ind = (1:9)' * logspace(0, log10(param.K)-1, log10(param.K));
% ind = [ind(:)', param.K];
ind = (1:9)' * logspace(0, log10(K)-1, log10(K));
ind = [ind(:)', K];
colors = [0, 0.45, 0.75; 0.85, 0.325, 0.01; 0.925, 0.70, 0.125; 0.50, 0.20, 0.55];
fig = figure;
% set(fig, 'Units', 'normalized', 'Position', [0.35, 0.25, 0.4, 0.55])
hold on

% tiledlayout(2,2, 'Padding', 'none', 'TileSpacing', 'compact'); 


% p1 = plot_with_shade(ind, sub_opt_main, 5, 0.1, colors(1,:));
p_non = plot_with_shade(ind, abs(non_reg_sub_opt), 5, 0.1, colors(1,:));
p_tik = plot_with_shade(ind, abs(tikhonov_reg_sub_opt), 5, 0.1, colors(2,:));
p_ent = plot_with_shade(ind, abs(entropic_reg_sub_opt), 5, 0.1, colors(3,:));
% p_ent2 = plot_with_shade(ind, abs(entropic2_reg_sub_opt), 5, 0.1, colors(2,:));

p_hyp = plot_with_shade(ind, abs(hyperbolic_reg_sub_opt), 5, 0.1, colors(4,:));


% ylim([10^-7, 10^0])
% xlim([min(ind), max(ind)])
set(gca, 'XScale', 'log', 'YScale', 'log');
set(gca, 'FontSize', 16);
xlabel('$\#$ of iterations ($T$)', 'FontSize', 20, 'interpreter','latex');
ylabel('Suboptimality','FontSize', 20, 'interpreter','latex')
grid on
% lgd = legend([p_ent, p_ent2], 'Entropic opt with agd',...
%     'Entropic opt with longer sgd',...
%     'Location', 'southwest', 'interpreter','latex');
lgd = legend([p_non, p_ent, p_tik, p_hyp], 'No regularization', 'Entropic regularization',...
    '$\chi^2$-divergence regularization', 'Hyperbolic divergence regularization',...
    'Location', 'southwest', 'interpreter','latex');
lgd.FontSize = 16;
ylim([1e-7, 1e-1])
xticks(logspace(0,5, 6))
yticks(logspace(-7, -1, 7))
remove_border()


saveas(gcf, 'convergence', 'svg')
% saveas(gcf, 'convergence', 'png')

