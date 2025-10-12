function res = get_fsc(V1,V2,pixelsize)

fsc = FSCorr(V1,V2);
fsc(1)=1;
n = numel(fsc);

plot(1:n,fsc,'r-*','LineWidth',2); % Plot FSC

xlim([1 n]);
ylim([-0.1 1.05]);
grid on

y=ones(1,n)*.5;
hold on; 
plot(1:n,y,'k--','Linewidth',1.5);
hold off;

j = fscres(fsc,.5);



res=2*pixelsize*n/j;


xticks=get(gca,'XTick');
df=1/(2*pixelsize*n);
xticks=xticks*df;
set(gca,'XTickLabel',sprintf('%7.3f\n',xticks));
fontsize(gcf,16,"points");
xlabel('$1/\textup{\AA}$', 'Interpreter', 'latex', 'FontSize', 20)
set(gca,'FontSize',20)
legend({strcat([num2str(res,'%.2f'), '$\textup{\AA}$'])},'Interpreter','latex','FontSize',30)



end