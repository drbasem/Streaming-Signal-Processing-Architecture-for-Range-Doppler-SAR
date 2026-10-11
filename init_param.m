function [rawData,H_ref,H_az1,H_az2] = init_param
%% ================= Radar & Platform Parameters =================
fc = 10e9;                  % Carrier frequency (Hz)
c  = 3e8;
lambda = c / fc;
fs  = 100e6;                % ADC sampling rate
prf = 1000;                 % Pulse Repetition Frequency
Tp  = 1e-6;                 % Chirp duration
B   = 50e6;                 % Chirp bandwidth
K   = B / Tp;               % Chirp rate
v   = 100;                  % Platform velocity (m/s)
numPulses  = 512;           % Azimuth samples
numSamples = 1024;          % Range samples
targets = [1000 0;1200 50]; % [range (m), azimuth (m)]
%% ================= Transmitted Chirp =================
t = (0:1/fs:Tp-1/fs);
txChirp = exp(1j*pi*K*t.^2);
%% ================= Raw ADC Data Simulation =================
rawData = zeros(numSamples, numPulses);
for p = 1:numPulses
    azPos = (p - numPulses/2) * v / prf;
    for k = 1:size(targets,1)
        R   = sqrt(targets(k,1)^2 + (azPos - targets(k,2))^2);
        tau = 2 * R / c;
        idx = round(tau * fs) + 1;
        if idx + length(txChirp) - 1 <= numSamples
            phase = exp(-1j * 4 * pi * fc * R / c);
            rawData(idx:idx+length(txChirp)-1, p) = ...
                rawData(idx:idx+length(txChirp)-1, p) + phase * txChirp.';
        end
    end
end
rawData = rawData + 0.01 * (randn(size(rawData)) + 1j * randn(size(rawData)));
%%
H_ref = fft(conj(fliplr(txChirp)), numSamples);
fd = (0:numPulses-1) * (prf / numPulses);
fd(fd > prf/2) = fd(fd > prf/2) - prf;
alphaR = c / (2 * fs);
R0 = (1:numSamples) * alphaR;
Ka = -2 * v^2 ./ (lambda * R0);
H_az1 = zeros(numSamples,numPulses);
H_az2 = zeros(numSamples,numPulses);
for r = 1:numSamples
    h = exp(-1j * pi * (fd.^2) / Ka(r));
    H_az1(r,:) = h;
    H_az2(r,:) = h(bitrevorder(1:512));
end


