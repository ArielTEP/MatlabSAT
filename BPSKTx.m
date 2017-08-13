clc
clear all
close all
% this is sth new hello MB
% Simulation of an impulse modulated raised-cosine BPSK signal generator
% input terms
Fs = 32000;          % sample frequency of simulation (Hz)
dataRate = 16000;     % data rate in bps
beta = 0.25;        % raised-cosine rolloff factor
symbols = 2;        % symbols per 2-bit (1-bit, per symbol)(symbol periods)
time = 3;       % length of signal in seconds
A = 1000;             %  scale factor
f = Fs/4;           % carrier frequency
cosine=[1 0 -1 0];  % cos(2*pi*f*n)
counter = 1;        % used to get a new data bit

% calculated terms
numberOfSamples = Fs*time;
samplesPerSymbol = Fs/dataRate;

% create rised-cosine filter
B = (rcosfir(beta, symbols/2, samplesPerSymbol, 1/Fs));
%B = rcosdesign(beta, symbols, samplesPerSymbol, 'sqrt');
Zi = zeros(1, (length(B) - 1));
%plot(B)
xlabel('Samples')

inputSamples = [1,1,1,1, 0,1,1,0,0,1,1,0,1 ,0,0,0, 1,1,1,1, 0,1,1,0,0,1,1,0,1 ,0,0,0, 1,1,1,1, 0,1,1,0,0,1,1,0,1 ,0,0,0, 1,1,1,1, 0,1,1,0,0,1,1,0,1 ,0,0,0, 1,1,1,1, 0,1,1,0,0,1,1,0,1 ,0,0,0, 1,1,1,1, 0,1,1,0,0,1,1,0,1 ,0,0,0, 1,1,1,1, 0,1,1,0,0,1,1,0,1 ,0,0,0, 1,1,1,1, 0,1,1,0,0,1,1,0,1 ,0,0,0, 1,1,1,1, 0,1,1,0,0,1,1,0,1 ,0,0,0, 1,1,1,1, 0,1,1,0,0,1,1,0,1 ,0,0,0];

% Real-Time (In serie simulation)
k=1;
j=1;
dataArray = inputSamples;
for i = 1:numberOfSamples
    % get the new data bit at the beginning of a symbol period
    if(counter == 1)
        data = A*(2*(rand > 0.5)-1);
        %data = A*(2*(inputSamples(j))-1);
        j=j+1;
        %dataArray = [dataArray data/A];
    else
        data = 0;
    end
    %sound(data,Fs)
    % pulse modulated signal
    [impulseModulatedData, Zf] = filter(B, 1, data, Zi);
    Zi = Zf;
    BPSKsignal(k) = impulseModulatedData*cosine(mod(i ,4) + 1) ;
    
    % reset at the end of a symbol period
    if(counter == samplesPerSymbol)
        counter = 0;
    end
    %sound(output,Fs)
    % increment counter
    counter = counter + 1;
    k=k+1;
end
plot((0:1/Fs:time-1/Fs),BPSKsignal)
X = abs(fft(BPSKsignal, 1024));
figure, plot(fftshift(X))
sound(BPSKsignal,Fs)

% --- Receiver

% simulation inputs PLL
alphaPLL = 0.010;   % PLL's loop filter parameter 'alpha'
betaPLL = 0.002;
N = samplesPerSymbol;
Fs = 32000;         % simulation sample frequency
phaseAccumPLL = randn(1);    % current phase accumulator's value
VCOphaseError = 2*pi*rand(1); % random phase error
VCOrestFrequencyError = randn(1); % error in the VCO's rest freq.
Fcarrier = 8000;
phi = VCOphaseError;    % initializing VCO's phase

% simulation inputs ... Timing Recovery
alphaML = 0.0040;   %ML loop filter parameter 'alpha'
betaML = 0.0002;    %ML loop filter parameter 'beta'
alpha = 0.25;       % root raised-cosine rolloff factor
symbols = 2;

phaseAccumML = 2*pi*rand(1);    % initializing the ML NCO
symbolsPerSecond = 16001;

% real time loop
for i = 1:length(BPSKsignal)
    % processing the data by the PLL
%     phaseDetectorOutput = analyticSignal(i)*vcoOutput;
%     m = 6 * real( phaseDetectorOutput ) ; % scale for a max value
%     q = real(phaseDetectorOutput)*imag(phaseDetectorOutput);
%     [loopFilterOutputPLL,Zi_pll]=filter(B_PLL, A_PLL, q, Zi_pll); 
%     loopFilterOutputPLLSummary = ...
%     [loopFilterOutputPLLSummary loopFilterOutputPLL] ; % plot 
%     phi = mod(phi + loopFilterOutputPLL + 2*pi*Fcarrier*T, 2*pi); 
%     vcoOutput = exp(-j*phi);
%     % processing the data by the ML-based receiver
%     [MFoutput, Zi_MF] = filter(B_MF, 1, m, Zi_MF); 
%     [diffMFoutput,Zi_diff]=filter([1 0 -1], 1, MFoutput,Zi_diff);
%    
%     phaseAccumML = phaseAccumML + phaseIncML; 
%     
%     if phaseAccumML >= 2*pi
%         phaseAccumML = phaseAccumML - 2*pi; 
%         decision = sign ( delayedMFoutput ) ;
%         [ error , Zi_ML_loop ] = filter (B_ML , A_ML , ...
%             decision*diffMFoutput , Zi_ML_loop ) ; 
%         phaseAccumML = phaseAccumML - error ;
%         errorSummary = [errorSummary error] ; %plot
%         decisionSummary = [ decisionSummary decision ] ;
%     else
%         errorSummary = [errorSummary 0];
%     end
%     
%     delayedMFoutput = MFoutput; % accounts for group delay

end

% % state storage for plotting ... not part of the ISR 
% delayedMFoutputSummary = . . .
% errorSummary = decisionSummary
% % plot % plot
%    errorSummary = [ errorSummary 0 ] ;
% % plot
%        [ delayedMFoutputSummary delayedMFoutput ] ; decisionSummaryHoldOn = [ decisionSummaryHoldOn
% decision ] ;
%   
