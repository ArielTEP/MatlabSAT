function [BPSKsignal, dataArray] = impModBPSK(time)
    % this is sth new hello MB
    % Simulation of an impulse modulated raised-cosine BPSK signal generator
    % input terms
    Fs = 32000;          % sample frequency of simulation (Hz)
    dataRate = 16000;     % data rate in bps
    beta = 0.25;        % raised-cosine rolloff factor
    symbols = 2;        % symbols per 2-bit (1-bit, per symbol)(symbol periods)
    %time = 3;       % length of signal in seconds
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
    %B = rcosdesign(beta, symbols, samplesPerSymbol, 'sqrt');
    Zi = zeros(1, (length(B) - 1));
    %plot(B)
    xlabel('Samples')

    inputSamples = [1,1,1,1, 0,1,1,0,0,1,1,0,1 ,0,0,0, 1,1,1,1, 0,1,1,0,0,1,1,0,1 ,0,0,0, 1,1,1,1, 0,1,1,0,0,1,1,0,1 ,0,0,0, 1,1,1,1, 0,1,1,0,0,1,1,0,1 ,0,0,0, 1,1,1,1, 0,1,1,0,0,1,1,0,1 ,0,0,0, 1,1,1,1, 0,1,1,0,0,1,1,0,1 ,0,0,0, 1,1,1,1, 0,1,1,0,0,1,1,0,1 ,0,0,0, 1,1,1,1, 0,1,1,0,0,1,1,0,1 ,0,0,0, 1,1,1,1, 0,1,1,0,0,1,1,0,1 ,0,0,0, 1,1,1,1, 0,1,1,0,0,1,1,0,1 ,0,0,0];
    isLength = length(inputSamples);
    addSamples = numberOfSamples - isLength;
    inputSamples = [inputSamples zeros(1,addSamples)];
    % Real-Time (In serie simulation)
    k=1;
    j=1;
    dataArray = inputSamples;
    for i = 1:numberOfSamples
        % get the new data bit at the beginning of a symbol period
        if(counter == 1)
            %data = A*(2*(rand > 0.5)-1);
            data = A*(2*(inputSamples(j))-1);
            j=j+1;
            dataArray = [dataArray data/A];
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
    xlim([0 0.001])
    X = abs(fft(BPSKsignal, 1024));
    figure, plot(fftshift(X))
    sound(BPSKsignal,Fs)

end