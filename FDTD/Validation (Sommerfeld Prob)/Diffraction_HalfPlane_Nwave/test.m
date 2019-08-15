function test()

N_duration = 0.5;
N_amplitude = 1;
N_slope = N_duration/2/N_amplitude;
tt = -1:0.1:2;

figure, plot(tt,Signature(tt));

function X = Heaviside(x)
    X = ones(size(x));
    I = find(x<0);
    X(I) = 0;
end

function res = Signature(S)
    res = Heaviside(N_duration-S).*(N_amplitude-S/N_slope);
end

end