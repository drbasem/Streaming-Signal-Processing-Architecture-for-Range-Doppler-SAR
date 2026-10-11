function [y, ww] = fft_dit(x, stage)
    m = 2^stage;
    y = zeros(1, m);
    ww = ones(1, m/2);
    wm = exp(-1j * 2 * pi / m);
    for k = 0:m:512-1  
        w = 1;
        for j = 0:(m/2 - 1)
            ww(j+1) = w;
            u = x(k + j + 1      );
            v = x(k + j + 1 + m/2) * w;  
            y(k + j + 1      ) = u + v;
            y(k + j + 1 + m/2) = u - v;
            w = w * wm;
        end
    end
end


