%% ideal model plots
close all;
clear all;

Ps = 4;
N0 = 4e-15;
h1_list = [4.24591e-13];
h2_list = [1.86826e-14];
real_700to_1900 = [ 4.5 4.0 3.0 3.0 3.0 2.6 2.0 2.0 2.0 1.5 1.5 1.5 1.3];
real_rate = [1.213675];
color = hsv(length(h1_list));

fg= 1;
figure(fg);
for idx=1:length(h1_list)
    h1 = h1_list(idx);
    h2 = h2_list(idx);
    GRatio(idx) = h1/h2;
    fprintf('g1/g2 = %5.2f \n', h1/h2);
    alpha = 0.00:0.01:1.00;
    ratio = 4; %2
    
    % plot
    R1_max = log2(1+Ps*h1/ N0);
    R2_max = log2(1+Ps*h2/ N0);
    for i = 1:length(alpha)
        P1 = Ps * alpha(i); % user 1 near
        P2 = Ps * (1-alpha(i)); % user 2 far
        if (P2>=P1)
            R1(i) = log2(1+P1*h1 / N0);
            R2(i) = log2(1+P2*h2 / (P1*h2 + N0));
        end
    end
    line(0.00:0.01:0.50 ,R1/R1_max+R2/R2_max, 'Color',color(idx,:),'LineWidth',2,'Marker','x');
    hold on;
    plot([0.09 0.1 0.11], [1.213675 1.213675 1.213675],'Color','b','LineWidth',2,'Marker','x');
    hold on;
end
set(gca, 'yScale', 'linear',...
    'yMinorTick','on',...
    'yLim', [1 1.6]);
L = legend('Theoretical', 'Simulation');
set(L,'FontSize',14);
xlabel('Power allocation factor (\alpha)','FontSize',14,'Color','k');
ylabel('Gain of Weighted Sum','FontSize',14,'Color','k');
set(gca,'fontsize',14);
box on;
grid on;
saveas(gcf,'alphaVSgain','epsc');