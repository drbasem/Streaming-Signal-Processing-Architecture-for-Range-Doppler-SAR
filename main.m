clc; clear;
[rawData,H_ref,H_az1,H_az2] = init_param;
%% classic
%% fft
x1 = zeros(size(rawData));
for i = 1:size(rawData,2)
    x1(:,i) = fft(rawData(:,i));
end
%% mul
x2 = zeros(size(x1));
for i = 1:size(x1,2)
    x2(:,i)  = x1(:,i) .* H_ref.';
end
%% ifft
x3 = zeros(size(x2));
for i = 1:size(x2,2)
    x3(:,i) = ifft(x2(:,i));
end
%% fft
x4 = zeros(size(x3));
for i = 1:size(x3,1)
    x4(i,:) = fft(x3(i,:));
end
%% mul
x5 = zeros(size(x4));
for i = 1:size(x3,1)
    x5(i,:) = x4(i,:) .* H_az1(i,:);
end
%% ifft
x6 = zeros(size(x5));
for i = 1:size(x3,1)
    x6(i,:) = ifft(x5(i,:));
end
%%
x7 = abs(x6);
%% porposed
%% fft
xx1 = zeros(size(rawData));
for i = 1:size(rawData,2)
    xx1(:,i) = fft(rawData(:,i));
end
%% mul
xx2 = zeros(size(xx1));
for i = 1:size(xx1,2)
    xx2(:,i) = xx1(:,i) .* H_ref.';
end
%% ifft
xx3 = zeros(size(xx2));
for i = 1:size(rawData,2)
    xx3(:,i) = ifft(xx2(:,i));
end
%% fft 
xx4 = zeros(size(xx3));
for i = 1:size(rawData,1)
    [xx4(i,:),w1,w2,w3,w4,w5,w6,w7,w8,w9] = stream_fft(xx3(i,:));
end
%% mul
xx5 = zeros(size(xx4));
for i = 1:size(xx4,1)
    xx5(i,:) = xx4(i,:).* H_az2(i,:);
end
%% fft
xx6 = zeros(size(xx5));
for i = 1:size(xx5,1)
    [xx6(i,:),ww1,ww2,ww3,ww4,ww5,ww6,ww7,ww8,ww9] = stream_fft_dit(xx5(i,:));
end
%%
xx7 = abs(xx6);
%%
subplot(2,1,1)
mesh(x7)
subplot(2,1,2)
mesh(xx7/512)
