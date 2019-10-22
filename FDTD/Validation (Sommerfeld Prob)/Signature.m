function res = Signature(S,fs,input_wave)
    res = zeros(size(S));
    I = find(S>0);
    res(I) = input_wave(ceil(S(I)*fs));
end
