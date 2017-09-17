clc
clear all
close all

Fs = 32000;          % sample frequency of simulation (Hz)
dataRate = 1600;     % data rate in bps
symbols = 2;        % symbols per 2-bit (1-bit, per symbol)(symbol periods)

A = 10;             %  scale factor
f = Fs/4;           % carrier frequency
cosine=[1 0 -1 0];  % cos(2*pi*f*n)
%cosine = cos(2*pi*f.*[1 2 3 4]);
counter = 1;        % used to get a new data bit

% calculated terms
time = 1;
numberOfSamples = Fs*time;
samplesPerSymbol = Fs/dataRate;

% create rised-cosine filter
beta = 0.25;
 B = (rcosfir(beta, symbols/2, samplesPerSymbol, 1/Fs));

% Generate the BPSK transmitter's signal
[BPSKsignal, dataArray] = impModBPSK(time);


% PLL Loop Filter settings
alphaPLL = 0.010;   % PLL loop filter "alpha"
betaPLL = 0.002;    % PLL loop filter "beta"

N = samplesPerSymbol;
% All down this, add realism to simulation ensuring initial errors

% PLL inits
phaseAccumPLL = rand(1);    % init phase accumulator's value
VCOphaseError = 2*pi*rand(1); % phase error random init
VCOrestFrequencyError = rand(1);    % freq error init
Fcarrier = f;
phi = VCOphaseError;    % init the VCO's phase
analyticSignal = hilbert(BPSKsignal);
%vcoOutput = [exp(-j*VCOphaseError) zeros(1,N)];
vcoOutput = exp(-j*VCOphaseError);
loopFilterOutput = zeros(1, N+1);
Ts = 1/Fs;
B_PLL = [(alphaPLL + betaPLL) -alphaPLL];
A_PLL = [1 -1];
Zi_pll = 0;
loopFilterOutputSummary = [];
                
% Maximum Likelyhood criteria for Timing Recovery
alphaML = 0.0050; % loop filter "alpha"
betaML = 0.0002;  % loop filter "beta"
alpha = 0.25;   % RRC roll-off factor
phaseIncML = pi/10;
B_ML = [(alphaML + betaML) -alphaML];
A_ML = [1 -1];

% init the ML NCO
phaseAccumML = 2*pi*rand(1);
symbolsPerSecond = dataRate + 1 ;   % symbol rate w/ offset from 1600


delay = samplesPerSymbol;
decisionSummary = []; % zeros(1,length(BPSKsignal)+delay);
errorSummary = [];
rcvData = [];
decision = [];

% Real Time Simulation
for i = 1:length(BPSKsignal)
    % ********** PLL *******************
    phaseDetectorOutput = analyticSignal(i)*vcoOutput(i);
    m = 7*real(phaseDetectorOutput);
    q = real(phaseDetectorOutput) * imag(phaseDetectorOutput);
    [loopFilterOutputPLL, Zf_pll] = filter(B_PLL, A_PLL, q);
    loopFilterOutputSummary = [loopFilterOutputSummary loopFilterOutputPLL];
    phi = mod(phi + loopFilterOutputPLL + 2*pi*Fcarrier*Ts, 2*pi);
    vcoOutput(i+1) = exp(-j*phi);
    % **********************************
    
    % Maximum Likelyhood for Time recovery **********
    [MFoutput, Zf_MF] = filter(B, 1, m);
    [diffMFoutput, Zf_diff] = filter([1 0 -1], 1, MFoutput);
    phaseAccumML = phaseAccumML + phaseIncML;
    
    delayedMFoutput = MFoutput;

    % Time recovery loop
    if phaseAccumML >= 2*pi
       phaseAccumML = phaseAccumML - 2*pi;
       decision = sign(delayedMFoutput);
       [error, Zi_ML_loop] = filter(B_ML, A_ML, decision*diffMFoutput);
       phaseAccumML = phaseAccumML - error;
       errorSummary = [errorSummary error]; % plot
       decisionSummary = [decisionSummary decision]; % plotc
    else
        errorSummary = [errorSummary 0];    % plot
    end
    % ***********************************************
    %delayedMFoutput = MFoutput;
end
figure
plot(real(vcoOutput))
Xr = abs(fft(vcoOutput));
figure, plot(fftshift(Xr))
title('Spectrum')
figure, plot(decisionSummary);
title('Decision Summary');
figure, plot(errorSummary);
title('errorSummary');


[decisionSummary(1:10)', downsample(dataArray(1:200),samplesPerSymbol)']
%dataArray(downsample(1:50,samplesPerSymbol))
%decisionSummary(downsample(1 + delay : 50 +delay,samplesPerSymbol))








