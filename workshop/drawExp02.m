% collect and plot for hw2
close all;
FILE = {'./exp02result'};
fg = 1;

figure(fg); fg = fg + 1;
dat = dlmread(FILE{1});
for j = 2:size(dat, 2)
    ln(j-1) = line(dat(:,1), dat(:,j),'Marker','o','LineWidth',2,'MarkerFaceColor','none','MarkerSize',8,'Color',[j/size(dat, 2) 0 0]);
end
set(gca, 'yScale', 'log', 'yMinorTick','on');
L = legend('Far end user', 'Near end user');
set(L,'FontSize',14);
xlabel('Distance (m)','FontSize',14,'Color','k');
ylabel('Psr (%)','FontSize',14,'Color','k');
set(gca,'fontsize',14);
grid on;
box on;
saveas(gcf,'fig1','epsc');
% 
% figure(fg); fg = fg + 1;
% dat = dlmread(FILE{2});
% for j = 2:size(dat, 2)
%     ln(j-1) = line(dat(:,1), dat(:,j)*100,'Marker','o','LineWidth',2,'MarkerFaceColor','none','MarkerSize',8,'Color',[j/size(dat, 2) 0 0]);
% end
% set(gca, 'XLim', [20 200]);
% L = legend('6Mbps','9Mbps','12Mbps','18Mbps','24Mbps','36Mbps','48Mbps','54Mbps');
% set(L,'FontSize',14);
% xlabel('Distance (m)','FontSize',14,'Color','k');
% ylabel('Psr (%)','FontSize',14,'Color','k');
% set(gca,'fontsize',14);
% grid on;
% box on;
% saveas(gcf,'fig2','epsc');
% 
% figure(fg); fg = fg + 1;
% dat = dlmread(FILE{3});
% for j = 2:size(dat, 2)
%     ln(j-1) = line(dat(:,1), dat(:,j)*100,'Marker','o','LineWidth',2,'MarkerFaceColor','none','MarkerSize',8,'Color',[j/size(dat, 2) 0 0]);
% end
% set(gca, 'XLim', [20 120]);
% L = legend('6Mbps','9Mbps','12Mbps','18Mbps','24Mbps','36Mbps','48Mbps','54Mbps');
% set(L,'FontSize',14);
% xlabel('Distance (m)','FontSize',14,'Color','k');
% ylabel('Psr (%)','FontSize',14,'Color','k');
% set(gca,'fontsize',14);
% grid on;
% box on;
% saveas(gcf,'fig3','epsc');
% 
% figure(fg); fg = fg + 1;
% dat = dlmread(FILE{4});
% for j = 2:size(dat, 2)
%     ln(j-1) = line(dat(:,1), dat(:,j)*100,'Marker','o','LineWidth',2,'MarkerFaceColor','none','MarkerSize',8,'Color',[j/size(dat, 2) 0 0]);
% end
% set(gca, 'XLim', [80 350]);
% L = legend('6Mbps','9Mbps','12Mbps','18Mbps','24Mbps','36Mbps','48Mbps','54Mbps');
% set(L,'FontSize',14);
% xlabel('Distance (m)','FontSize',14,'Color','k');
% ylabel('Psr (%)','FontSize',14,'Color','k');
% set(gca,'fontsize',14);
% grid on;
% box on;
% saveas(gcf,'fig4','epsc');
