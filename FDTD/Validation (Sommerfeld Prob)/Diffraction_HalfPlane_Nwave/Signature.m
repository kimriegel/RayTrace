function res = Signature(S,N_duration,N_amplitude,N_slope)
    res = Heaviside(S).*Heaviside(N_duration-S).*(N_amplitude-S/N_slope);
end