function [azFFT,w1,w2,w3,w4,w5,w6,w7,w8,w9] = stream_fft(xx)

    [y9,w9] = fft_dif(xx,9);
    [y8,w8] = fft_dif(y9,8);   
    [y7,w7] = fft_dif(y8,7);
    [y6,w6] = fft_dif(y7,6);
    [y5,w5] = fft_dif(y6,5);
    [y4,w4] = fft_dif(y5,4);
    [y3,w3] = fft_dif(y4,3);
    [y2,w2] = fft_dif(y3,2);
    [y1,w1] = fft_dif(y2,1);
    azFFT = y1;