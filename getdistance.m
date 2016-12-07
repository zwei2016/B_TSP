function d = getdistance (lat1, lon1, lat2, lon2)
% calculat the distance in the earth between two cities 
% method inspired by good map

EARTH_RADIUS = 6378.137; %km
a = lat1*pi/180 - lat2*pi/180;
b = lon1*pi/180 - lon2*pi/180;
d = EARTH_RADIUS*2*asin(sqrt(sin(a/2)^2) + cos(lat1*pi/180)*cos(lat2*pi/180)*sin(b/2)^2);

end

