%% fusome plots

cell_fus_vol_fracs_mean = [0.214 0.176 0.088 0.091 0.059 0.058 0.051 0.054 0.027 0.027 0.022 0.028 0.024 0.027 0.024 0.030];

cell_fus_vol_fracs_stdev = [0.0376 .0388 0.0282 0.0311 0.0169 0.0213 0.0174 0.0161 0.0200 0.0112 0.0116 0.0179 0.0198 0.0253 0.0142 0.0193];

figure; hold on; box on;

errorbar(1:16, cell_fus_vol_fracs_mean, cell_fus_vol_fracs_stdev, 'ko', 'MarkerSize',10, 'LineWidth',1.5)
axis([0.5 16.5 0 0.3])
h = gca;
h.FontSize = 20;
xlabel('Cell $\#$','interpreter','latex')
ylabel('Fusome Vol. Frac.','interpreter','latex')

%% fusome volume - orb correlation plots

fus_volfrac_12 = [0.600 0.516 0.448 0.548 0.571 0.623 0.456 0.430 0.498 0.508];

orb_frac_12 = [1 1 0.989 0.972 1 0.704 0.903 0.890 1 1];

figure; box on;
plot(orb_frac_12, fus_volfrac_12./(1-fus_volfrac_12),'ko','LineWidth',2,'MarkerSize',10); hold on;
plot(linspace(0,1.5,100), ones(100,1),'k--','LineWidth',1.5)
axis([0.5 1.05 0.7  1.7])
ylabel('``Winner" Fusome Volume Ratio','interpreter','latex')
xlabel('``Winner" \textit{orb} Fraction','interpreter','latex')
axis square;
h = gca;
h.FontSize = 20;