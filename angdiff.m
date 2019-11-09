function angdiff = angdiff(dec1, inc1, dec2, inc2)

if length(inc1(:, 1)) == 1
    inc1 = inc1';
    dec1 = dec1';
    inc2 = inc2';
    dec1 = dec1';
end

x1 = cosd(inc1).*cosd(dec1);
y1 = cosd(inc1).*sind(dec1);
z1 = sind(inc1);

x2 = cosd(inc2).*cosd(dec2);
y2 = cosd(inc2).*sind(dec2);
z2 = sind(inc2);

t1 = cross([x1 y1 z1], [x2 y2 z2]);
t2 = (t1(:, 1).^2 + t1(:, 2).^2 + t1(:, 3).^2).^.5;
t3 = dot([x1 y1 z1]', [x2 y2 z2]');
angdiff = atan2d(t2, t3');