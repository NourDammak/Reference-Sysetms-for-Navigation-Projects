
% Mapping of all points
 
gn=importdata('data_to_anomalies.txt');
phi=gn.data(:,1);
lambda=gn.data(:,2);
H=gn.data(:,3);
g=gn.data(:,4);

%map of all points
figure('Name','map')
geoscatter(phi, lambda, H, 'filled','MarkerFaceAlpha', 0.3);
c = colorbar;
c.Location = 'eastoutside'
c.Label.String = 'H [m]';
title('all points locations')
colormap(jet);
geobasemap topographic;
geolimits([48.5 55],[13 25])


%map of sample points (n°3,n°53,n°103,n°183 and n°233)

gn=importdata('test_pnts.txt');
phi=gn.data(:,1);
lambda=gn.data(:,2);
H=gn.data(:,3);
g=gn.data(:,4);

figure('Name','map')
geoscatter(phi, lambda, H, 'filled','MarkerFaceAlpha', 0.3);
c = colorbar;
c.Location = 'eastoutside'
c.Label.String = 'H [m]';
title('test points locations')
colormap(jet);
geobasemap topographic;
geolimits([48.5 55],[13 25])


%mapping of point n°103

phi=54.70619;
lambda=17.72920;
H=120.964;
g=981448.392;

figure('Name','map')
geoscatter(phi, lambda, H, 'filled','MarkerFaceAlpha', 0.3);
c = colorbar;
c.Location = 'eastoutside'
c.Label.String = 'H [m]';
title('point 103 location')
colormap(jet);
geobasemap topographic;
geolimits([48.5 55],[13 25])




 
 



