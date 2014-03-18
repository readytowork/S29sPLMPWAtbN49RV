%% ------------------- RECEIVER ------------------------
function [ RCVD_DATA ] = sicReceiver( SIGNAL_IN, CONFIG_IN, SIC_LEVEL)
%% Removing Cyclic Extension
for i=1:64
    rxed_sig(i)=SIGNAL_IN(i+16);
end

%% FFT
ff_sig=fft(rxed_sig,64);

%% Pilot Sync
for i=1:52
    synched_sig1(i)=ff_sig(i+6);
end

k=1;
for i=(1:13:CONFIG_IN.nPilot*13)
    for j=(i+1:i+12);
        synched_sig(k)=synched_sig1(j);
        k=k+1;
    end
end
% scatterplot(synched_sig)

%% Demodulation
divided_sig = 0;
% accRcvdPower = abs(synched_sig).^2;
% avgRcvdPower_ref = sum(accRcvdPower) / length(accRcvdPower);
% SIC_LEVEL
% avgRcvdPower_ref
% synched_sig = synched_sig / (avgRcvdPower^0.5) * (20^0.5); % Normalized, 10 is the average power for default 16QAM
% accRcvdPower = abs(synched_sig).^2;
% avgRcvdPower = sum(accRcvdPower) / length(accRcvdPower)
% return
for i = 1:SIC_LEVEL
    synched_sig = synched_sig - divided_sig;
    
% accRcvdPower = abs(synched_sig).^2;
% avgRcvdPower = sum(accRcvdPower) / length(accRcvdPower);
% synched_sig = synched_sig / (avgRcvdPower^0.5) * (avgRcvdPower_ref^0.5); % Normalized, 10 is the average power for default 16QAM

    dem_data = qamdemod(synched_sig, 16, 0, 'gray');
    divided_sig = qammod(dem_data, 16, 0, 'gray');
end

%% Decimal to binary conversion
bin = de2bi(dem_data','left-msb');
bin = bin';

%% De-Interleaving
if CONFIG_IN.isInterleaved
    deintlvddata = matdeintrlv(bin,2,2); % De-Interleave
    deintlvddata=deintlvddata';
    deintlvddata=deintlvddata(:)';
else
    deintlvddata=bin';
    deintlvddata=deintlvddata(:)';
end

%% Decoding data
if CONFIG_IN.isCoded
    n=6;
    k=3;
    decodedata = vitdec(deintlvddata,trellis,5,'trunc','hard');  % decoding datausing veterbi decoder
    rxed_data=decodedata;
else
    rxed_data=deintlvddata;
end

RCVD_DATA=rxed_data;
end

