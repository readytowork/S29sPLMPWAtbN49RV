%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All rights reserved by Krishna Pillai, http://www.dsplog.com 
% The file may not be re-distributed without explicit authorization
% from Krishna Pillai.
% Checked for proper operation with Octave Version 3.0.0
% Author	: Krishna Pillai
% Email		: krishna@dsplog.com
% Version	: 1.0
% Date		: 01 January 2012
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Computing the Symbol Error probability for a general M-QAM 
% using OFDM modulation

function compute_symbol_error_rate_qam_ofdm_awgn()

close all; figure
EsN0dB 	= [0:33]; % symbol to noise ratio
M 	= [16 64 256]; % 16QAM/64QAM and 256 QAM  
color_vec1 = {'b-','m-','g-'};    
color_vec2 = {'ks-','rx-','cd-'};    
for (jj= 1:length(M))
	k      = sqrt(1/((2/3)*(M(jj)-1))); 
	simSer(jj,:) = compute_symbol_error_rate(EsN0dB, M(jj));
	theorySer(jj,:) = 2*(1-1/sqrt(M(jj)))*erfc(k*sqrt((10.^(EsN0dB/10)))) ...
	              - (1-2/sqrt(M(jj)) + 1/M(jj))*(erfc(k*sqrt((10.^(EsN0dB/10))))).^2;
	semilogy(EsN0dB,theorySer(jj,:),color_vec1(jj));
	hold on
	semilogy(EsN0dB,simSer(jj,:),   color_vec2(jj));
end
axis([0 33 10^-5 1])
grid on
legend('theory-16QAM', 'sim-16QAM', 'theory-64QAM', 'sim-64QAM', 'theory-256QAM', 'sim-256QAM');
xlabel('Es/No, dB')
ylabel('Symbol Error Rate')
title('Symbol error probability curve for 16QAM/64QAM/256QAM using OFDM');
return ;

function [simSer] = compute_symbol_error_rate(EsN0dB, M);

% ofdm specifications
nFFT = 64; % fft size
nDSC = 52; % number of data subcarriers
nConstperOFDMsym = 52; % number of bits per OFDM symbol (same as the number of subcarriers for BPSK)
nOFDMsym = 10^4; % number of ofdm symbols

% modulation 
k = sqrt(1/((2/3)*(M-1))); % normalizing factor
m = [1:sqrt(M)/2]; % alphabets
alphaMqam = [-(2*m-1) 2*m-1]; 

EsN0dB_eff = EsN0dB  + 10*log10(nDSC/nFFT) + 10*log10(64/80); % accounting for the used subcarriers and cyclic prefix

for ii = 1:length(EsN0dB)

   % Transmitter
   ipMod = randsrc(1,nConstperOFDMsym*nOFDMsym,alphaMqam) + j*randsrc(1,nConstperOFDMsym*nOFDMsym,alphaMqam);
   ipMod_norm = k*reshape(ipMod,nConstperOFDMsym,nOFDMsym).'; % grouping into multiple symbolsa

   % Assigning modulated symbols to subcarriers from [-26 to -1, +1 to +26]
   xF = [zeros(nOFDMsym,6) ipMod_norm(:,[1:nConstperOFDMsym/2]) zeros(nOFDMsym,1) ipMod_norm(:,[nConstperOFDMsym/2+1:nConstperOFDMsym]) zeros(nOFDMsym,5)] ;
   
   % Taking FFT, the term (nFFT/sqrt(nDSC)) is for normalizing the power of transmit symbol to 1 
   xt = (nFFT/sqrt(nDSC))*ifft(fftshift(xF.')).';

   % Appending cylic prefix
   xt = [xt(:,[49:64]) xt];

   % Concatenating multiple symbols to form a long vector
   xt = reshape(xt.',1,nOFDMsym*80);

   % Gaussian noise of unit variance, 0 mean
   nt = 1/sqrt(2)*[randn(1,nOFDMsym*80) + j*randn(1,nOFDMsym*80)];

   % Adding noise, the term sqrt(80/64) is to account for the wasted energy due to cyclic prefix
   yt = sqrt(80/64)*xt + 10^(-EsN0dB_eff(ii)/20)*nt;

   % Receiver
   yt = reshape(yt.',80,nOFDMsym).'; % formatting the received vector into symbols
   yt = yt(:,[17:80]); % removing cyclic prefix

   % converting to frequency domain
   yF = (sqrt(nDSC)/nFFT)*fftshift(fft(yt.')).'; 
   yMod = sqrt(64/80)*yF(:,[6+[1:nConstperOFDMsym/2] 7+[nConstperOFDMsym/2+1:nConstperOFDMsym] ]); 

   % demodulation
   y_re = real(yMod)/k;
   y_im = imag(yMod)/k;
   % rounding to the nearest alphabet
   % 0 to 2 --> 1
   % 2 to 4 --> 3
   % 4 to 6 --> 5 etc
   ipHat_re = 2*floor(y_re/2)+1;
   ipHat_re(find(ipHat_re>max(alphaMqam))) = max(alphaMqam);
   ipHat_re(find(ipHat_re<min(alphaMqam))) = min(alphaMqam);
            
   % rounding to the nearest alphabet
   % 0 to 2 --> 1
   % 2 to 4 --> 3
   % 4 to 6 --> 5 etc
   ipHat_im = 2*floor(y_im/2)+1;
   ipHat_im(find(ipHat_im>max(alphaMqam))) = max(alphaMqam);
   ipHat_im(find(ipHat_im<min(alphaMqam))) = min(alphaMqam);
    
   ipHat = ipHat_re + j*ipHat_im; 

   % converting to vector 
   ipHat_v = reshape(ipHat.',nConstperOFDMsym*nOFDMsym,1).';

   % counting the errors
   nErr(ii) = size(find(ipMod - ipHat_v ),2);

end
simSer = nErr/(nOFDMsym*nConstperOFDMsym);

return;







