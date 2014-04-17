%% ideal model plots
close all;
clear all;

Ps = 5e-3;
N0 = 4e-15;
dist_1 = 150;%150;
dist_2 = 220;%220;
h1 = (dist_1^-4.5);
h2 = (dist_2^-4.5);
alpha = 0.0:0.01:1.0;

% plot
for i = 1:length(alpha)
    P1 = Ps * alpha(i);
    P2 = Ps * (1-alpha(i));
    if (P1>P2)
        snr1 = P1*h1/(P2*h1+N0);
        snr2 = P2*h2/N0;
        R1(i) = log2(1+snr1);
        R2(i) = log2(1+snr2);
    else
        R1(i) = log2(1 + P1*h1 / N0);
        R2(i) = log2(1 + P2*h2 / (P1*h2 + N0));
    end
end

scatter(R1(:),R2(:))
