function [y9, w1, w2, w3, w4, w5, w6, w7, w8, w9] = stream_fft_dit(xx)
    [y1, w1] = fft_dit(xx, 1); 
    [y2, w2] = fft_dit(y1, 2);
    [y3, w3] = fft_dit(y2, 3);
    [y4, w4] = fft_dit(y3, 4);
    [y5, w5] = fft_dit(y4, 5);
    [y6, w6] = fft_dit(y5, 6);
    [y7, w7] = fft_dit(y6, 7);
    [y8, w8] = fft_dit(y7, 8);
    [y9, w9] = fft_dit(y8, 9);
end