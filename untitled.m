%% --- Receiver

% input terms
Fs = 32000;								% sample frequency of simulation (Hz)
dataRate = 16000;     					% data rate in bps
beta = 0.25;        					% raised-cosine rolloff factor
symbols = 2;        					% symbols per 2-bit (1-bit, per symbol)(symbol periods)
time = 3;       						% length of signal in seconds
A = 10;             					% scale factor
f = Fs/4;           					% carrier frequency
cosine=[1 0 -1 0];  					% cos(2*pi*f*n)
%cosine = cos(2*pi*f.*[1 2 3 4]);
counter = 1;        					% used to get a new data bit

% calculated terms
numberOfSamples = Fs*time;
samplesPerSymbol = Fs/dataRate;                        
B = (rcosfir(beta, symbols/2, samplesPerSymbol, 1/Fs));

filename = 'prueba.wav'; 						% Nombre archivo
recorder = audiorecorder(Fs,16,1); 				% Creacion del objeto de grabacion
%msgbox('Empezando Grabacion',' Grabadora '); 	% Mensaje de informacion
recordblocking(recorder,4); 					% Grabacion del sonido
%msgbox('Terminando Grabacion',' Grabadora '); 	% Mensaje de informacion
audio = getaudiodata(recorder, 'single'); 		% Paso los valores del objeto a una se?al
audiowrite(filename,audio,Fs); 					% Grabamos y guardamos la se?al
[x, Fs] = audioread(filename); 					% Leemos el sonido

X = fft(x); 							% Transformada rapida de fourier de x
X = X(1:end/2);
f = linspace(0,Fs/2,length(X)); 		% Frecuencia
t = 0:1/Fs:(length(x)-1)/Fs; 			% Vector Tiempo

%%

Fs = 32000;         % sample frequency of simulation (Hz)
dataRate = 16000;   % data rate in bps
beta = 0.25;        % raised-cosine rolloff factor
symbols = 2;        % symbols per 2-bit (1-bit, per symbol)(symbol periods)

A = 10;             %  scale factor
f = Fs/4;           % carrier frequency
cosine=[1 0 -1 0];  % cos(2*pi*f*n)
%cosine = cos(2*pi*f.*[1 2 3 4]);
counter = 1;        % used to get a new data bit

% calculated terms
numberOfSamples = Fs*time;
samplesPerSymbol = Fs/dataRate;

% create rised-cosine filter
B = (rcosfir(beta, symbols/2, samplesPerSymbol, 1/Fs));

% generate the BPSK transmitter's signal
[BPSK, dataArray] = impModBPSK(0.1);

alphaPLL = 0.010;       % PLL's loop filter parameter "alpha"
betaPLL = 0.002;        % PLL's loop filter parameter "beta"
N = samplesPerSymbol;

% PLL inits
phaseAccumPLL = rand(1);
VCOphaseError = 2*pi*rand(1);       % Phase error random init
VCOrestFrecuencyError = rand(1);    % Frequency error init
Fcarrier = 8000;                    % Carrier freq of the transmitter
phi = VCOphaseError;                % initializing the VCO's phase

% simulation inputs ... ML timing recovery
alphaML = 0.0040;       % ML loop filter parameter "alpha"
betaML = 0.0002;        % ML loop filter parameter "beta"
alpha = 0.25;           % root raised-cosine rolloff factor

phaseAccumML = 2*pi*rand(1);        % initializing the ML NCO
symbolsPerSecond = dataRate + 1;    % symbol rate w/ offset from 16000

