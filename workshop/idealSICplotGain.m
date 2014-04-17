%% ideal model plots
close all;
clear all;

Ps = 5e-3;
N0 = 4e-15;
dist_1 = 150;%150;
dist_2 = [250.2 417.5 696.5 1161.5];%220;

fg= 1;
figure(fg);
for distIndex=1:length(dist_2)
    h1 = (dist_1^-4.5);
    h2 = (dist_2(distIndex)^-4.5);
    GRatio(distIndex) = h1/h2;
    fprintf('g1/g2 = %5.2f \n', h1/h2);
    alpha = 0.00:0.01:1.00;
    ratio = 4; %2
    R1_12far = [];
    % R2_12far = [];
    % R1_12near = [];
    % R2_12near = [];
    % R1_21far = [];
    R2_21far = [];
    % R1_21near = [];
    % R2_21near = [];

    % plot
    for i = 1:length(alpha)
        P1 = Ps * alpha(i);
        P2 = Ps * (1-alpha(i));
        if (P1>P2)
            snr1 = P1*h1/(P2*h1+N0);
            snr2 = P2*h2/N0;
            if (P1>ratio*P2)
                R1_12far(end+1) = log2(1+snr1);
    %             R2_12far(end+1) = log2(1+snr2);
            else
    %             R1_12near(end+1) = log2(1+snr1);
    %             R2_12near(end+1) = log2(1+snr2);
            end
        else
            R1(i) = log2(1 + P1*h1 / N0);
            R2(i) = log2(1 + P2*h2 / (P1*h2 + N0));
            if (P2>ratio*P1)
    %             R1_21far(end+1) = log2(1 + P1*h1 / N0);
                R2_21far(end+1) = log2(1 + P2*h2 / (P1*h2 + N0));
            else
    %             R1_21near(end+1) = log2(1 + P1*h1 / N0);
    %             R2_21near(end+1) = log2(1 + P2*h2 / (P1*h2 + N0));
            end
        end
    end
    line(0.00:0.01:0.50 ,R1/R1_12far(end)+R2/R2_21far(1), 'Color',[1/distIndex^0.5 0 0],'LineWidth',2,'Marker','x');
end
set(gca, 'yScale', 'linear',...
    'yMinorTick','on',...
    'yLim', [1 1.6]);
L = legend('PL Ratio = 10', '100', '1000', '10000');
set(L,'FontSize',14);
xlabel('\alpha','FontSize',14,'Color','k');
ylabel('Gain of Weighted Sum','FontSize',14,'Color','k');
set(gca,'fontsize',14);
box on;
grid on;
saveas(gcf,'alphaVSgain','epsc');