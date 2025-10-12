clear 
clc 


load('test4.mat','res2_all', 'res3_all')
L = [4,6,8,10,12];

figure 
plot(L,res2_all,'-*','LineWidth',2)
hold on 
plot(L,res3_all,'-*','LineWidth',2)
grid on
set(gca,'fontsize',20)
xlabel('$L$','Interpreter','latex', 'FontSize', 30)
ylabel('Resolutions($\textup{\AA}$)','Interpreter','latex', 'FontSize', 30)
set(gcf, 'PaperPositionMode', 'auto');
fig = gcf;
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
% print(gcf, 'resolutions_L.pdf', '-dpdf', '-r0');
% saveas(gcf, 'resolutions_L.eps');
print('resolutions_L', '-depsc');

% set(gcf, 'PaperSize', [6.25 7.5]);
% set(gcf, 'PaperPositionMode', 'manual');
% set(gcf, 'PaperPosition', [0 0 6.25 7.5]);

% set(gcf, 'PaperUnits', 'inches');
% set(gcf, 'PaperSize', [6.25 7.5]);
% set(gcf, 'PaperPositionMode', 'manual');
% set(gcf, 'PaperPosition', [0 0 6.25 7.5]);

% set(gcf, 'renderer', 'painters');
% print(gcf, '-dpdf', 'resolutions_L.pdf');


load('test3_L=4.mat')
x2_time(1) = sum(out.x2_time);
x3_time(1) = out.x3_time;


load('test3_L=6.mat')
x2_time(2) = sum(out.x2_time);
x3_time(2) = out.x3_time;


load('test3_L=8.mat')
x2_time(3) = sum(out.x2_time);
x3_time(3) = out.x3_time;


load('test3_L=10.mat')
x2_time(4) = sum(out.x2_time);
x3_time(4) = out.x3_time;


load('test3_L=12.mat')
x2_time(5) = sum(out.x2_time);
x3_time(5) = out.x3_time;

figure 
plot(L,x2_time/3600,'-*','LineWidth',2)
hold on 
plot(L,x3_time/3600,'-*','LineWidth',2)
grid on
set(gca,'fontsize',20)
xlabel('$L$','Interpreter','latex', 'FontSize', 30)
ylabel('hours','Interpreter','latex', 'FontSize', 30)
set(gcf, 'PaperPositionMode', 'auto');
fig = gcf;
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
% print(gcf, 'runningtime_L.pdf', '-dpdf', '-r0');
print('runningtime_L', '-depsc');




