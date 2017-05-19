% Script for loading and resampling the Yosemite DEM

[Z,R] = usgsdem('mariposa-w',1);

lat0 = 37;
lon0 = -120;

dlat = 1/1200;
dlon = 1/1200;

lon = lon0:dlon:lon0+1;
lat = lat0:dlat:lat0+1;

dlat2 = 0.025;
dlon2 = 0.025;

lon2 = lon0+.025/2:dlon2:lon0+1-.025/2;
lat2 = lat0+.025/2:dlat2:lat0+1-.025/2;

[X,Y] = meshgrid(lon,lat);
[Xi,Yi] = meshgrid(lon2,lat2);

Z2 = interp2(X,Y,Z,Xi,Yi);
