function stdv = EMstd(x)

y = x;
iter=5;
for i = 1:iter
    z = zscore(y);
    mask = abs(z) > 3;
    y(mask) = mean(y(~mask));
end

stdv = std(y);