clear all 
clc 
run('~/aspire/initpath.m')
addpath(genpath('../src/'))
addpath '../EMDB_Data/'

load('result_N_1e4.mat')
load('result_N_1e5.mat')
load('result_N_1e6.mat')



res2(1,1) = res2_1e4(1);
res2(1,2) = res2_1e4(2);
res2(1,3) = res2_1e4(3);
res2(2,1) = res2_1e5(1);
res2(2,2) = res2_1e5(2);
res2(2,3) = res2_1e5(3);
res2(3,1) = res2_1e6(1);
res2(3,2) = res2_1e6(2);
res2(3,3) = res2_1e6(3);


res3(1,1) = res3_1e4(1);
res3(1,2) = res3_1e4(2);
res3(1,3) = res3_1e4(3);
res3(2,1) = res3_1e5(1);
res3(2,2) = res3_1e5(2);
res3(2,3) = res3_1e5(3);
res3(3,1) = res3_1e6(1);
res3(3,2) = res3_1e6(2);
res3(3,3) = res3_1e6(3);

N = [1e4,1e5,1e6];


% figure 
% semilogx(N, res2(1:3,1), '-*', 'LineWidth', 2)
% hold on 
% semilogx(N, res3(1:3,1), '-o', 'LineWidth', 2)
% grid on 
% set(gca,'fontsize',20)
% xlabel('$N$','Interpreter','latex', 'FontSize', 30)
% ylabel('Resolutions($\textup{\AA}$)','Interpreter','latex', 'FontSize', 30)
% saveas(gcf,'SNR=1.pdf')


% figure 
% semilogx(N, res2(1:3,2), '-*', 'LineWidth', 2)
% hold on 
% semilogx(N, res3(1:3,2), '-o', 'LineWidth', 2)
% grid on 
% set(gca,'fontsize',20)
% xlabel('$N$','Interpreter','latex', 'FontSize', 30)
% ylabel('Resolutions($\textup{\AA}$)','Interpreter','latex', 'FontSize', 30)
% saveas(gcf,'SNR=0.1.pdf')



% figure 
% semilogx(N, res2(1:3,3), '-*', 'LineWidth', 2)
% hold on 
% semilogx(N, res3(1:3,3), '-o', 'LineWidth', 2)
% grid on 
% set(gca,'fontsize',20)
% xlabel('$N$','Interpreter','latex', 'FontSize', 30)
% ylabel('Resolutions($\textup{\AA}$)','Interpreter','latex', 'FontSize', 30)
% saveas(gcf,'SNR=0.01.pdf')


% figure 
% semilogx(N, res2(1:3,1), '-*', 'LineWidth', 2)
% hold on 
% semilogx(N, res2(1:3,2), '-o', 'LineWidth', 2)
% grid on 
% semilogx(N, res2(1:3,3), '-+', 'LineWidth', 2)
% grid on
% set(gca,'fontsize',20)
% xlabel('$N$','Interpreter','latex', 'FontSize', 30)
% ylabel('Resolutions($\textup{\AA}$)','Interpreter','latex', 'FontSize', 30)
% legend('SNR=1','SNR=0.1','SNR=0.01')
% saveas(gcf,'m2_res.pdf')



figure 
semilogx(N, res3(1:3,1), '-*', 'LineWidth', 2)
hold on 
semilogx(N, res3(1:3,2), '-o', 'LineWidth', 2)
hold on 
semilogx(N, res3(1:3,3), '-+', 'LineWidth', 2)
grid on
set(gca,'fontsize',20)
xlabel('$N$','Interpreter','latex', 'FontSize', 30)
ylabel('Resolutions($\textup{\AA}$)','Interpreter','latex', 'FontSize', 30)
legend({'$\mathrm{SNR}=1$','$\mathrm{SNR}=0.1$','$\mathrm{SNR}=0.01$'},'Interpreter','latex', 'FontSize', 30)
set(gcf, 'PaperPositionMode', 'auto');
fig = gcf;
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
% print(gcf, 'm3_res.pdf', '-dpdf', '-r0');
print('m3_res', '-depsc');



V = double(ReadMRC('V.mrc'));
V_alms = double(ReadMRC('V_alms.mrc'));
pixel_size = 1.04*196/64;
res_best = get_fsc(V,V_alms,pixel_size);

pixelsize=  pixel_size;



fsc1 = FSCorr(V,V2_1e6{3});
fsc1(1)=1;
fsc2 = FSCorr(V,V3_1e6{3});
fsc2(1)=1;
fsc3 = FSCorr(V,V_alms);
fsc3(1)=1;
n = numel(fsc1);

plot(1:n,fsc1,'r-*','LineWidth',2); % Plot FSC
hold on;
plot(1:n,fsc2,'b-*','LineWidth',2); % Plot FSC
hold on;
plot(1:n,fsc3,'k-*','LineWidth',2); % Plot FSC
hold off;

xlim([1 n]);
ylim([-0.1 1.05]);
grid on

y=ones(1,n)*.5;
hold on; 
plot(1:n,y,'k--','Linewidth',1.5);
hold off;

j1 = fscres(fsc1,.5);
j2 = fscres(fsc2,.5);
j3 = fscres(fsc3,.5);

res1=2*pixelsize*n/j1;
res2=2*pixelsize*n/j2;
res3=2*pixelsize*n/j3;

xticks=get(gca,'XTick');
df=1/(2*pixelsize*n);
xticks=xticks*df;
set(gca,'XTickLabel',sprintf('%7.3f\n',xticks));
fontsize(gcf,16,"points");
xlabel('$1/\textup{\AA}$', 'Interpreter', 'latex', 'FontSize', 30)
ylabel('FSC', 'Interpreter', 'latex', 'FontSize', 30)
set(gca,'FontSize',20)
legend({strcat([num2str(res1,'%.2f'),'$\textup{\AA}$']), strcat([num2str(res2,'%.2f'), '$\textup{\AA}$']), strcat([num2str(res3,'%.2f'), '$\textup{\AA}$'])},'Interpreter','latex','FontSize',30)


set(gcf, 'PaperPositionMode', 'auto');
fig = gcf;
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
% print(gcf, 'FSC.pdf', '-dpdf', '-r0');
print('FSC', '-depsc');