function [x y z] = xyz(dec, inc, int);

x = int.*cosd(inc).*cosd(dec);
y = int.*cosd(inc).*sind(dec);
z = int.*sind(inc);

[l1 l2] = size(inc);
if l1 == 1
    x = x';
    y = y';
    z = z';
end
