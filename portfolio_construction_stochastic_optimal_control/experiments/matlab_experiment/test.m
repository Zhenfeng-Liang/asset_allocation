y = @(x) x^2-2;

error = 1;
tol = 0.0001;
yOld = -2;

while error>tol
     yNew = feval(y, yOld);
     error = abs(yNew - yOld);
     yOld = yNew;
end

yOld;

fh = @(y) y-y^2+2;

fplot(fh,[-5,5]);











