% Sanity check
% Last updated: 03/21/2021

% Make a normal distribution
pd = makedist('Normal');

% Generate distribution values for specific percentiles
percentiles = 0.01:0.01:0.99;
p = zeros(size(percentiles));
for i = 1:length(percentiles)
    p(i) = icdf(pd, percentiles(i));
end

% "Convert" to log normal with mean/std dev of interest
mean = 10;
stddev = 2;
pValues = zeros(size(percentiles));
for i = 1:length(percentiles)
    pValues(i) = 10^(p(i)*log10(stddev)+log10(mean));
end
