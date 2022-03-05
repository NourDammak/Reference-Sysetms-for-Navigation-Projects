format long
clear

%upload time series from file with phi lam h in radians
%Convert to degrees


ts=load('phi_lam_h_Time_Series.txt');
time=ts(:,1);
phi=ts(:,2);
la=ts(:,3);
h=ts(:,4);


for i=1:length(time)
    la_deg(i) = rad2deg(la(i));
    phi_deg(i) = rad2deg(phi(i));
end








