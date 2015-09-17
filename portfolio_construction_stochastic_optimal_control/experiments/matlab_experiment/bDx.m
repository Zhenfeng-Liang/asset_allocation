function y = bDx(x, sig)
    syms x;
    y = diff(b(x, sig), x);
end