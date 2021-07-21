clear
load('preload_variation_01.mat');
load('preload_variation_008.mat');
load('preload_variation_004.mat');

load('preload_variation_05.mat');
load('preload_variation_075.mat');

clf
plot(Vve_01(30000:end), Pve_01(30000:end), 'Linewidth', 2.5)
hold on
plot(Vve_075(30000:end), Pve_075(30000:end), 'Linewidth', 2.5)
hold on
plot(Vve_05(30000:end), Pve_05(30000:end), 'Linewidth', 2.5)
hold on
Vve_points = [Vve_01(index_01) Vve_075(index_075) Vve_05(index_05)];
Pve_points = [Pve_01(index_01) Pve_075(index_075) Pve_05(index_05)];
plot(Vve_points,Pve_points, 'o', 'Color', '#777777', 'Linewidth', 2.5)
% hold on
% plot(Vve_01(index_01), Pve_01(index_01), '.w')
% hold on
% plot(Vve_008(index_008), Pve_008(index_008), '.w')
% hold on
% plot(Vve_004(index_004), Pve_004(index_004), '.w')
xlim([0 130])
ylim([0 120])
grid on
set(gca,'fontname','Times New Roman')
xlabel('$V_{ve}$ (ml)','interpreter','latex')
ylabel('$P_{ve}$ (mmHg)','interpreter','latex')
legend('$R_s = 1.0$', '$R_s = 0.75$',   '$R_s = 0.5$', 'Aortic Valve Close', 'ESPVR', 'Location', 'Northwest', 'interpreter','latex')
set(gca,'FontSize', 21)
%%
Vve_points = [Vve_01(index_01) Vve_008(index_008) Vve_004(index_004)];
Pve_points = [Pve_01(index_01) Pve_008(index_008) Pve_004(index_004)];
preload_figura(Vve_01(30000:end), Pve_01(30000:end), Vve_008(30000:end), Pve_008(30000:end), Vve_004(30000:end), Pve_004(30000:end), Vve_points, Pve_points)