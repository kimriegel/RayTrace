function X = Heaviside(x)
    X = ones(size(x));
    I = find(x<0);
    X(I) = 0;
end
