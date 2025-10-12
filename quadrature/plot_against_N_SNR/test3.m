clear 
clc 

addpath(genpath('../src/'))
addpath '../EMDB_Data/'

run('~/aspire/initpath.m')
V = double(ReadMRC('emd_34948.map'));
V = V/norm(V(:));
n = 64;
V = vol_downsample(V,n);
L = 10;
c = 0.5;
[L,S] = get_truncate_limit(L, n, c, 1);
N_basis = get_num_basis(L,S);
vol_coef = vol_t_vol_coeffs(V, L, S, c, 1);
V_alms = vol_coeffs_t_vol(vol_coef, n, c, L, S, 1);
pixelsize = 1.04*196/n;


% WriteMRC(V, 1, 'V.mrc')


load('result_N_1e6.mat')
V1 = V3_1e6{1}; 
WriteMRC(V1, 1, 'V3_1e6_1.mrc')
V2 = V3_1e6{2};
WriteMRC(V2, 1, 'V3_1e6_0.1.mrc')
V3 = V3_1e6{3};
WriteMRC(V3, 1, 'V3_1e6_0.01.mrc')

fsc1 = FSCorr(V_alms,V1);
fsc1(1)=1;
fsc2 = FSCorr(V_alms,V2);
fsc2(1)=1;
fsc3 = FSCorr(V_alms,V3);
fsc3(1)=1;
n = numel(fsc1);

plot(1:n,fsc1,'r-*','LineWidth',2); % Plot FSC

hold on;
plot(1:n,fsc2,'b-*','LineWidth',2); % Plot FSC

hold on;
plot(1:n,fsc3,'g-*','LineWidth',2); % Plot FSC
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
xlabel('$1/\textup{\AA}$', 'Interpreter', 'latex', 'FontSize', 20)
set(gca,'FontSize',20)
% legend({strcat([num2str(res1,'%.2f'),'$\textup{\AA}$']), strcat([num2str(res2,'%.2f'), '$\textup{\AA}$'])},'Interpreter','latex','FontSize',30)
legend({'SNR=1','SNR=0.1','SNR=0.01'},'Interpreter','latex','FontSize',30)
title('$N=1\times 10^6$','Interpreter', 'latex', 'FontSize', 30)

saveas(gcf, 'fsc_1e6.pdf')