% [out1 out2] = psvr2(target, candidate)
%
% inputs:
% three columns (depth, dec, inc)
% depth vector for target should be primary depth
% depth vector for candidate should be warped vector
%
% outputs:
% out1 = cacualted with sum of squares of x y and z components
% out2 = calculated with cosine distances
%
% requires code xyz, decinc, and angdiff
%
% Brendan Reilly, Oregon State, January 2019



function [out1, out2] = psvr2(target, candidate)

%Organize Data
cdep = candidate(:, 1);
cdec = candidate(:, 2);
cinc = candidate(:, 3);

tdep = nan(size(cdep));
tdec = nan(size(cdep));
tinc = nan(size(cdep));

%Pair target depths with candidate depths
for n = 1:length(cdep)
    i = find(target(:, 1) == cdep(n)); %%Why values that don't match??? Check with Cedric!!
    if isempty(i)
        tdep(n) = NaN;
        tdec(n) = NaN;
        tinc(n) = NaN;
    else
    tdep(n) = target(i, 1);
    tdec(n) = target(i, 2);
    tinc(n) = target(i, 3);
    end
end

% Catch NaNs...
i = ~isnan(tdep);
tdep = tdep(i);
tdec = tdec(i);
tinc = tinc(i);
cdep = cdep(i);
cdec = cdec(i);
cinc = cinc(i);

t = nan(length(tdep), 3);
c = nan(length(cdep), 3);

%Calculate vector components
[t(:,1), t(:, 2), t(:, 3)] = xyz(tdec, tinc, ones(size(tdec)));
[c(:, 1), c(:, 2), c(:, 3)] = xyz(cdec, cinc, ones(size(cdec)));

% Preallocate
tss = [nan nan nan];
rss = tss;
tbar = tss;
cbar = tss;
m = tss;
b = tss;
e = nan(length(cdec), 3);
mod = e;

clf

for n = 1:3
    
    %caculate means
    tbar(n) = mean(t(:, n));
    cbar(n) = mean(c(:, n));
    
    %calculate slope
    m(n) = sum((t(:, n)-tbar(n)).*(c(:, n)-cbar(n)))/sum((t(:, n)-tbar(n)).^2);
    
    %calculate intercept
    b(n) = cbar(n) - m(n)*tbar(n);
    
    %calcualte model
    mod(:, n) = m(n)*t(:, n) + b(n);
    
    %calculate residuals
    e(:, n) = c(:, n) - b(n) - (m(n)*t(:, n));
    
    %total sum of squares
    tss(n) = sum((c(:, n) - cbar(n)).^2);
    
    %residual sum of squares
    rss(n) = sum(e(:, n).^2);
end

%calucalte model inc and dec
[moddec, modinc] = decinc(mod(:, 1), mod(:, 2), mod(:, 3));

%calcualte mean candidate inc and dec
[mdec, minc] = decinc(cbar(1), cbar(2), cbar(3));

%calculate the residual cosine distance and total cosine distance
rcosd = sum(1-cosd(angdiff(moddec, modinc, cdec, cinc)));
tcosd = sum(1-cosd(angdiff(ones(size(cdec))*mdec, ones(size(cinc))*minc, cdec, cinc)));

% sum euclidean
out1 = 1 - (rss(1) + rss(2) + rss(3))/(tss(1) + tss(2) + tss(3));

%cosine distance
out2 = 1 - rcosd/tcosd;