%% Beginner's Tutorial: Streaming Range-Doppler SAR (MATLAB Path)

% Run this live script from the repository folder (it needs init_param.m, stream_fft.m, stream_fft_dit.m, fft_dif.m, fft_dit.m on the path, and the Signal Processing Toolbox for bitrevorder). Run it section by section with Ctrl+Enter. The companion document TUTORIAL_MATLAB (Word/PDF/Markdown) explains the theory in more detail.

%% 1. The big picture

% A SAR records one column of echoes per radar pulse. Rows are range (fast time), columns are azimuth (slow time). Focusing is two compressions: range compression (FFT, multiply by a matched filter, IFFT down each column) and azimuth compression (the same along each row). The classic path uses MATLAB's fft/ifft. The streaming path rebuilds the azimuth FFTs from 9 hand-written butterfly stages, the way an FPGA would.

%% 2. Build the simulated radar data

% init_param simulates two point targets and returns the raw echoes plus three filters: H_ref (range), H_az1 (azimuth, natural order) and H_az2 (azimuth, bit-reversed order for the streaming path).

clear; clc;
[rawData,H_ref,H_az1,H_az2] = init_param;
size(rawData)   % 1024 range samples x 512 pulses

% Raw data looks like noise: the echoes of the two targets are long chirps smeared over many range samples.

imagesc(abs(rawData)); xlabel('Pulse (azimuth)'); ylabel('Range sample'); title('Raw data magnitude'); colorbar;

%% 3. Classic path

% Step 1-3: range compression, column by column. H_ref is a row vector, so .' turns it into a column before the element-wise multiply.

x1 = zeros(size(rawData));
for i = 1:size(rawData,2)
    x1(:,i) = fft(rawData(:,i));      % range FFT
end
x2 = zeros(size(x1));
for i = 1:size(x1,2)
    x2(:,i) = x1(:,i) .* H_ref.';     % matched filter
end
x3 = zeros(size(x2));
for i = 1:size(x2,2)
    x3(:,i) = ifft(x2(:,i));          % range IFFT
end

% After range compression each target is a sharp line in range (expect rows near 766 and 900: the echo start plus about one chirp length of filter lag).

imagesc(abs(x3)); xlabel('Pulse (azimuth)'); ylabel('Range sample'); title('After range compression'); colorbar;

% Step 4-6: azimuth compression, row by row, each row with its own filter.

x4 = zeros(size(x3));
for i = 1:size(x3,1)
    x4(i,:) = fft(x3(i,:));           % azimuth FFT
end
x5 = zeros(size(x4));
for i = 1:size(x3,1)
    x5(i,:) = x4(i,:) .* H_az1(i,:);  % azimuth matched filter
end
x6 = zeros(size(x5));
for i = 1:size(x3,1)
    x6(i,:) = ifft(x5(i,:));          % azimuth IFFT
end
x7 = abs(x6);

%% 4. Streaming path

% The range steps are identical, so they are repeated here. The azimuth FFT is built from 9 DIF stages (natural order in, bit-reversed order out), so the filter must be bit-reversed (H_az2). Then 9 DIT stages (bit-reversed in, natural out) finish the job without any reordering memory in between.

xx3 = x3;   % range compression is the same in both paths
xx4 = zeros(size(xx3));
for i = 1:size(xx3,1)
    [xx4(i,:),w1,w2,w3,w4,w5,w6,w7,w8,w9] = stream_fft(xx3(i,:));
end
xx5 = xx4 .* H_az2;
xx6 = zeros(size(xx5));
for i = 1:size(xx5,1)
    [xx6(i,:),ww1,ww2,ww3,ww4,ww5,ww6,ww7,ww8,ww9] = stream_fft_dit(xx5(i,:));
end
xx7 = abs(xx6)/512;

% Sanity check of the DIF chain alone: its output must equal MATLAB's fft in bit-reversed order.

x = randn(1,512) + 1j*randn(1,512);
[y,~] = stream_fft(x);
ref = fft(x);
max(abs(y - ref(bitrevorder(1:512))))   % should be about 1e-12

%% 5. Compare the two paths

subplot(2,1,1); mesh(x7);        title('Classic');
subplot(2,1,2); mesh(xx7);       title('Streaming (/512)');

% The last streaming stage uses forward twiddles, so it is a forward FFT where the classic path has an inverse FFT. The result is scaled by 512 (removed above) and reversed along azimuth. A Python port of this pipeline gave about 0.59 relative error when compared directly, and about 2e-16 after flipping the azimuth axis. Check it in MATLAB:

errDirect  = norm(x7 - xx7,'fro') / norm(x7,'fro')
xxFlip     = xx7(:, [1 512:-1:2]);
errFlipped = norm(x7 - xxFlip,'fro') / norm(x7,'fro')

% If errFlipped is tiny and errDirect is not, the streaming output is a mirrored copy of the classic one. A fix is to use exp(+1j*2*pi/m) as the twiddle in fft_dit.m, or to conjugate the input and output of the last FFT.

%% 6. Exercises

% 1) Remove the second target row from targets in init_param.m and rerun. 2) Move the first target from 1000 m to 1100 m: the range peak should shift by about 67 rows (100 m / 1.5 m). 3) Time both paths with tic/toc and explain why the hand-written stages are slower in MATLAB. 4) Apply the twiddle fix and confirm errDirect becomes tiny.
