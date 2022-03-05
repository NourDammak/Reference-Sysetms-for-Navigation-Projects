clear
%load the data from gravimetric measurements, we need phi, lambda,H and
%gravity, as the separately vectors

gn=importdata('data_to_anomalies.txt');
fi=gn.data(:,1);
la=gn.data(:,2);
H=gn.data(:,3);
g=gn.data(:,4);

% load anomalies from analysed models
m1_B=importdata('EIGEN-6C4_Bouguer.dat')
B_m1=m1_B.data(:,4);

m2_B=importdata('ITSG-Grace2018s_Bouguer.dat')
B_m2=m2_B.data(:,4);

m3_B=importdata('XGM2019e_2159_Bouguer.dat')
B_m3=m3_B.data(:,4);

m1_Fa=importdata('EIGEN-6C4_Free_Air.dat')
Fa_m1=m1_B.data(:,4);

m2_Fa=importdata('ITSG-Grace2018s_Free_Air.dat')
Fa_m2=m2_B.data(:,4);

m3_Fa=importdata('XGM2019e_2159_Free_Air.dat')
Fa_m3=m3_B.data(:,4);




% modeled density of the Earth's surface topographic mass
sigma=2.67

%normal gravity formula GRS80 depending phi
for i=1:length(fi)
gamma(i,1)=978032.66*(1+0.0053024*(sin(fi(i)*pi/180)^2-0.00000585*(sin(2*fi(i)*pi/180)^2)));
end

%compute gravimetric reductions Rfa, RBg and anomalies
for i=1:length(H)
Rfa(i,1)=0.30855*H(i);
RBg(i,1)=-0.04192*sigma*H(i);
Afa(i,1)=g(i)+Rfa(i)-gamma(i);
ABg(i,1)=g(i)+Rfa(i)+RBg(i)-gamma(i);
end

%compute of residuals of models
res_B1(:,1)=ABg(:,1)-B_m1(:,1);
res_B2(:,1)=ABg(:,1)-B_m2(:,1);
res_B3(:,1)=ABg(:,1)-B_m3(:,1);
res_Fa1(:,1)=Afa(:,1)-Fa_m1(:,1);
res_Fa2(:,1)=Afa(:,1)-Fa_m2(:,1);
res_Fa3(:,1)=Afa(:,1)-Fa_m3(:,1);


%maps of residuals - countur plot and 3D plot
x = min(la(:,1)):0.1:max(la(:,1)); 
y = min(fi(:,1)):0.1:max(fi(:,1)); 

[XI,YI] = meshgrid(x,y);

RB1 =  griddata(la(:,1),fi(:,1),res_B1(:,1),XI,YI,'cubic');
RB2 =  griddata(la(:,1),fi(:,1),res_B2(:,1),XI,YI,'cubic');
RB3 =  griddata(la(:,1),fi(:,1),res_B3(:,1),XI,YI,'cubic');

RFa1 =  griddata(la(:,1),fi(:,1),res_Fa1(:,1),XI,YI,'cubic');
RFa2 =  griddata(la(:,1),fi(:,1),res_Fa2(:,1),XI,YI,'cubic');
RFa3 =  griddata(la(:,1),fi(:,1),res_Fa3(:,1),XI,YI,'cubic');



stat_B_m1 = datastats(res_B1(:,1))
stat_B_m2 = datastats(res_B2(:,1))
stat_B_m3 = datastats(res_B3(:,1))

stat_Fa_m1 = datastats(res_Fa1(:,1))
stat_Fa_m2 = datastats(res_Fa2(:,1))
stat_Fa_m3 = datastats(res_Fa3(:,1))

iso1=ceil(stat_B_m1.range);
iso2=ceil(stat_B_m2.range);
iso3=ceil(stat_B_m3.range);

iso4=ceil(stat_Fa_m1.range);
iso5=ceil(stat_Fa_m2.range);
iso6=ceil(stat_Fa_m3.range);


figure('Name','Bouguer anomalies residuals')
 subplot(1,3,1)
 contourf(XI,YI,RB1,iso1), colorbar('southoutside'), hold on
 title('model EIGEN-6C4 residuals')
 xlabel(['\lambda', '[deg]'])
 ylabel(['\phi', '[deg]'])

  subplot(1,3,2)
 contourf(XI,YI,RB2,iso2), colorbar('southoutside'), hold on
 title('model ITSG-Grace2018s residuals')
 xlabel(['\lambda', '[deg]'])
 ylabel(['\phi', '[deg]'])

   subplot(1,3,3)
 contourf(XI,YI,RB3,iso3), colorbar('southoutside'), hold on
 title('model XGM2019e_2159 residuals')
 xlabel(['\lambda', '[deg]'])
 ylabel(['\phi', '[deg]'])

 
 
 figure('Name','Bouguer anomalies residuals')
 subplot(1,3,1)
 surfc(XI,YI,RB1)
 title('model EIGEN-6C4 residuals')
 xlabel(['\lambda', '[deg]'])
 ylabel(['\phi', '[deg]'])
 
 subplot(1,3,2)
 surfc(XI,YI,RB2)
 title('model ITSG-Grace2018s residuals')
 xlabel(['\lambda', '[deg]'])
 ylabel(['\phi', '[deg]'])

  subplot(1,3,3)
 surfc(XI,YI,RB3)
 title('model XGM2019e_2159 residuals')
 xlabel(['\lambda', '[deg]'])
 ylabel(['\phi', '[deg]'])




  figure('Name','Free-air anomalies residuals')
 subplot(1,3,1)
 contourf(XI,YI,RFa1,iso4), colorbar('southoutside'), hold on
 title('model EIGEN-6C4 residuals')
 xlabel(['\lambda', '[deg]'])
 ylabel(['\phi', '[deg]'])

  subplot(1,3,2)
 contourf(XI,YI,RFa2,iso5), colorbar('southoutside'), hold on
 title('model ITSG-Grace2018s residuals')
 xlabel(['\lambda', '[deg]'])
 ylabel(['\phi', '[deg]'])

   subplot(1,3,3)
 contourf(XI,YI,RFa3,iso6), colorbar('southoutside'), hold on
 title('model XGM2019e_2159 residuals')
 xlabel(['\lambda', '[deg]'])
 ylabel(['\phi', '[deg]'])

 
 
   figure('Name','Free-air anomalies residuals')
 subplot(1,3,1)
 surfc(XI,YI,RFa1)
 title('model EIGEN-6C4 residuals')
 xlabel(['\lambda', '[deg]'])
 ylabel(['\phi', '[deg]'])
 
 subplot(1,3,2)
 surfc(XI,YI,RFa2)
 title('model ITSG-Grace2018s residuals')
 xlabel(['\lambda', '[deg]'])
 ylabel(['\phi', '[deg]'])

  subplot(1,3,3)
 surfc(XI,YI,RFa3)
 title('model XGM2019e_2159 residuals')
 xlabel(['\lambda', '[deg]'])
 ylabel(['\phi', '[deg]'])

