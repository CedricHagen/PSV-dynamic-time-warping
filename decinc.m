function [dec, inc, int] = decinc(x, y, z);

int = sqrt(x.^2+y.^2+z.^2);

dec = atan2(y, x) * 180/pi();

inc = asind(z./int);

