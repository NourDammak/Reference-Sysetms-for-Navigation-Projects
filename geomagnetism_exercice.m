%%
clear

%load the data from igrf13coeffs_data.txt file containing dipole model g1_0, g1_1 and h1_1 exctacting from igrf13coeffs.txt
igrf=importdata('igrf13coeffs_data.txt');
year=igrf(1,:);
g1_0=igrf(2,:);
g1_1=igrf(3,:);
h1_1=igrf(4,:);

%lets suppose λ=21 and ϕ=52 --> θ=90-52=38
lat=52
co_lat=90-52
lon=21

%compute magnetic intensity X, Y and Z for point ϕ and λ as a function of year
for i=1:length(g1_0)
X(i,1)=-g1_0(i)*sin(co_lat*pi/180)+(((g1_1(i)*cos(lon*pi/180))+(h1_1(i)*sin(lon*pi/180)))*cos(co_lat*pi/180))
Y(i,1)=(g1_1(i)*sin(lon*pi/180))-(h1_1(i)*cos(lon*pi/180))
Z(i,1)=-2*((g1_0(i)*cos(co_lat*pi/180))+(((g1_1(i)*cos(lon*pi/180))+(h1_1(i)*sin(lon*pi/180)))*sin(co_lat*pi/180)))
end

%compute the Declination
for i=1:length(year)
D(i,1)=atand(Y(i)/X(i))
end

%compute the Inclination
for i=1:length(year)
I(i,1)=atand(Z(i)/sqrt((X(i)^2)+(Y(i)^2)))
end


%visualisations of magnetic intensity (3 components) and declination, and inclination

figure('Name','Plots of magnetic intensity (X, Y and Z componement)')
 subplot(1,3,1)
 plot(year,X);
 title(' X = f(year)')
 xlabel(' year ')
 ylabel(' X [nanoTesla] ') 
 grid on
  subplot(1,3,2)
 plot(year,Y);
 title(' Y = f(year) ')
 xlabel(' year ')
 ylabel(' Y [nanoTesla] ') 
 grid on
  subplot(1,3,3)
 plot(year,Z); 
 title(' Z = f(year)')
 xlabel(' year ')
 ylabel('Z [nanoTesla]') 
 grid on

figure('Name','Plots of declination and inclination')
subplot(1,2,1)
plot(year,D)
title(' Declination = f(year) ')
xlabel(['year'])
ylabel('D [°]')
grid on

subplot(1,2,2)
plot(year,I)
title(' Inclination = f(year) ')
xlabel('year')
ylabel('I [°]')
grid on
  