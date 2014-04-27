%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All rights reserved by Krishna Sankar, http://www.dsplog.com
% The file may not be re-distributed without explicit authorization
% from Krishna Sankar.
% Checked for proper operation with Octave Version 3.0.0
% Author        : Krishna Sankar M 
% Email         : krishna@dsplog.com
% Version       : 1.0
% Date          : 24 January 2010
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Script for computing the BER for BPSK modulation in 3 tap ISI 
% channel. Minimum Mean Square Error (MMSE) equalization with 7 tap 
% and the BER computed (and is compared with Zero Forcing equalization)

clear
N  = 10^6; % number of bits or symbols
Eb_N0_dB = [0:15]; % multiple Eb/N0 values
K = 3;

for ii = 1:length(Eb_N0_dB)

   % Transmitter
   ip = rand(1,N)>0.5; % generating 0,1 with equal probability
   s = 2*ip-1; % BPSK modulation 0 -> -1; 1 -> 0 

   % Channel model, multipath channel
   nTap = 3;
   ht = [1 0.8 0.64]; 
   L  = length(ht);

   chanOut = conv(s,ht);  
   n = 1/sqrt(2)*[randn(1,N+length(ht)-1) + j*randn(1,N+length(ht)-1)]; % white gaussian noise, 0dB variance 
   
   % Noise addition
   y = chanOut + 10^(-Eb_N0_dB(ii)/20)*n; % additive white gaussian noise

   
   % zero forcing equalization
   hM = toeplitz([ht([2:end]) zeros(1,2*K+1-L+1)], [ ht([2:-1:1]) zeros(1,2*K+1-L+1) ]);
   d  = zeros(1,2*K+1);
   d(K+1) = 1;
   c_zf  = [inv(hM)*d.'].';
   yFilt_zf = conv(y,c_zf);
   yFilt_zf = yFilt_zf(K+2:end); 
   yFilt_zf = conv(yFilt_zf,ones(1,1)); % convolution
   ySamp_zf = yFilt_zf(1:1:N);  % sampling at time T
 
   % mmse equalization
   hAutoCorr = conv(ht,fliplr(ht));
   hM = toeplitz([hAutoCorr([3:end]) zeros(1,2*K+1-L)], [ hAutoCorr([3:end]) zeros(1,2*K+1-L) ]);
   hM
   hM = hM + 1/2*10^(-Eb_N0_dB(ii)/10)*eye(2*K+1);
   d  = zeros(1,2*K+1);
   d([-1:1]+K+1) = fliplr(ht);
   c_mmse  = [inv(hM)*d.'].';
   yFilt_mmse = conv(y,c_mmse);
   yFilt_mmse = yFilt_mmse(K+2:end); 
   yFilt_mmse = conv(yFilt_mmse,ones(1,1)); % convolution
   ySamp_mmse = yFilt_mmse(1:1:N);  % sampling at time T
  
   % receiver - hard decision decoding
   ipHat_zf = real(ySamp_zf)>0;
   ipHat_mmse = real(ySamp_mmse)>0;

   % counting the errors
   nErr_zf(1,ii) = size(find([ip- ipHat_zf]),2);
   nErr_mmse(1,ii) = size(find([ip- ipHat_mmse]),2);

   

end

simBer_zf = nErr_zf/N; % simulated ber
simBer_mmse = nErr_mmse/N; % simulated ber
theoryBer = 0.5*erfc(sqrt(10.^(Eb_N0_dB/10))); % theoretical ber

% plot
close all
figure
semilogy(Eb_N0_dB,simBer_zf(1,:),'bs-','Linewidth',2);
hold on
semilogy(Eb_N0_dB,simBer_mmse(1,:),'gd-','Linewidth',2);
axis([0 14 10^-5 0.5])
grid on
legend('sim-zf', 'sim-mmse');
xlabel('Eb/No, dB');
ylabel('Bit Error Rate');
title('Bit error probability curve for BPSK in ISI with MMSE equalizer');







