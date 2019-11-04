function Theta = atan2_mod(X,Y)
    Theta = atan2(X,Y);
    I = find(Theta<0);
    Theta(I) = Theta(I)+2*pi;
end
