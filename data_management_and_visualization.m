clear
%load the data from gravimetric measurements, we need phi, lambda,H and
%gravity, as the separately vectors

gn=importdata('data_to_anomalies.txt');
fi=gn.data(:,1);
la=gn.data(:,2);
H=gn.data(:,3);
g=gn.data(:,4);

% modeled density of the Earth's surface topographic mass
sigma=2.25


%normal gravity formula GRS80 depending phi
for i=1:length(fi)
gamma(i,1)=978032.66*(1+0.0053024*(sin(fi(i)*pi/180)^2-0.00000585*(sin(2*fi(i)*pi/180)^2)));
end


%compute gravimetric reductions Rfa, RBg, RPP

for i=1:length(H)
    Rfa(i,1)=0.30855*H(i)
    RBg(i,1)=-0.04192*sigma*H(i)
    RPP(i,1)=(0.30888-0.08384*sigma)*H(i)
    
end
    
%compute gravity anomalies Afa,ABg, APP

for i=1:length(H)
    Afa(i,1)=g(i)+Rfa(i)-gamma(i)
    ABg(i,1)=g(i)+RBg(i)+Rfa(i)-gamma(i)
    APP(i,1)=g(i)+RPP(i)-gamma(i)
    end

%visualisations of the anomalies as function of H and determine empirical
%function of dependency Ag=f(H) al linear interpolation model

figure('Name',' gravity anomalies vs Height')
subplot(2,2,1)
plot(fi,gamma,'x')
title(' gamma = f(fi) ')
xlabel([' \phi', '[deg]'])
ylabel('gravity [mGals]')

subplot(2,2,2)
p1 = polyfit(H,Afa,1);
y1 = polyval(p1,H);
plot(H,Afa,'o')
hold on
plot(H,y1)
hold off
title([' Afa = f(H)  ',num2str(p1(1,1),'%.4f'),' mGal/m'])
xlabel(' H ')
ylabel(' Afa ')                 

subplot(2,2,3)
p2 = polyfit(H,ABg,1);
y2 = polyval(p2,H);
plot(H,ABg,'o')
hold on
plot(H,y2)
hold off
title([' ABg = f(H)  ',num2str(p2(1,1),'%.4f'),' mGal/m'])
xlabel(' H ')
ylabel(' ABg ') 


subplot(2,2,4)
p3 = polyfit(H,APP,1);
y3 = polyval(p3,H);
plot(H,APP,'o')
hold on
plot(H,y3)
hold off
title([' APP = f(H) ',num2str(p3(1,1),'%.4f'),' mGal/m'])
xlabel(' H ')
ylabel(' APP ') 



%corelation coefficient between H and anomalies
cF=corrcoef(H,Afa)
cBg=corrcoef(H,ABg)
cPP=corrcoef(H,APP)


%maps - countur plot and 3D plot
x = min(la(:,1)):0.1:max(la(:,1)); 
y = min(fi(:,1)):0.1:max(fi(:,1)); 

[XI,YI] = meshgrid(x,y);

Zfa =  griddata(la(:,1),fi(:,1),Afa(:,1),XI,YI,'cubic');
ZBg =  griddata(la(:,1),fi(:,1),ABg(:,1),XI,YI,'cubic');
ZPp =  griddata(la(:,1),fi(:,1),APP(:,1),XI,YI,'cubic');

figure('Name','Free-air anomalies')
 subplot(1,2,1)
 contourf(XI,YI,Zfa,20), colorbar, hold on
 xlabel(['\lambda', '[deg]'])
 ylabel(['\phi', '[deg]'])

 subplot(1,2,2)
 surfc(XI,YI,Zfa)
 xlabel(['\lambda', '[deg]'])
 ylabel(['\phi', '[deg]'])

 figure('Name','Bouguer anomaly')
 subplot(1,2,1)
 contourf(XI,YI,ZBg,20), colorbar, hold on
 xlabel(['\lambda', '[deg]'])
 ylabel(['\phi', '[deg]'])

 subplot(1,2,2)
 surfc(XI,YI,ZBg)
 xlabel(['\lambda', '[deg]'])
 ylabel(['\phi', '[deg]'])

  figure('Name','Poincare-Prey anomaly')
 subplot(1,2,1)
 contourf(XI,YI,ZPp,20), colorbar, hold on
 xlabel(['\lambda', '[deg]'])
 ylabel(['\phi', '[deg]'])

 subplot(1,2,2)
 surfc(XI,YI,ZPp)
 xlabel(['\lambda', '[deg]'])
 ylabel(['\phi', '[deg]'])
 

 % interpolation methods
 % sample point number 3 in the test_points
 
 fi_3=54.28769
 la_3=19.08789
 H_3=1.216
 g_3=981437.974



 % 1. using of empirical H relationship
 Afa_3_i1=p1(1,1)*H_3+p1(1,2)
 ABg_3_i1=p2(1,1)*H_3+p2(1,2)
 APP_3_i1=p3(1,1)*H_3+p3(1,2)

 
 

 
 % 2 spatial interpolation from phi and lambda
 Afa_3_i2=griddata(la(:,1),fi(:,1),Afa(:,1),la_3,fi_3,'linear')
 ABg_3_i2=griddata(la(:,1),fi(:,1),ABg(:,1),la_3,fi_3,'linear')
 APP_3_i2=griddata(la(:,1),fi(:,1),APP(:,1),la_3,fi_3,'linear')

 
 % measured gravity anomalies at test point
 Afa_3=g_3+0.30855*H_3-(978032.66*(1+0.0053024*((sin(fi_3*pi/180))^2-0.00000585*((sin(2*fi_3*pi/180))^2))))
 ABg_3=g_3+(-0.04192*sigma*H_3)+(0.30855*H_3)-(978032.66*(1+0.0053024*((sin(fi_3*pi/180))^2-0.00000585*((sin(2*fi_3*pi/180))^2))))
 APP_3=g_3+((0.30888-0.08384*sigma)*H_3)-(978032.66*(1+0.0053024*((sin(fi_3*pi/180))^2-0.00000585*((sin(2*fi_3*pi/180))^2))))
 
 % residuals of two method interpolation
 resfa1_3=Afa_3-Afa_3_i1
 resfa2_3=Afa_3-Afa_3_i2

 resB1_3=ABg_3-ABg_3_i1
 resB2_3=ABg_3-ABg_3_i2

 resPP1_3 =APP_3-APP_3_i1
 resPP2_3 =APP_3-APP_3_i2

 
 % interpolation methods
 % sample point number 53 in the test_points
 
 fi_53=54.53464
 la_53=17.42375
 H_53=94.961
 g_53=981437.948



 % 1. using of empirical H relationship
 Afa_53_i1=p1(1,1)*H_53+p1(1,2)
 ABg_53_i1=p2(1,1)*H_53+p2(1,2)
 APP_53_i1=p3(1,1)*H_53+p3(1,2)

 
 

 
 % 2 spatial interpolation from phi and lambda
 Afa_53_i2=griddata(la(:,1),fi(:,1),Afa(:,1),la_53,fi_53,'linear')
 ABg_53_i2=griddata(la(:,1),fi(:,1),ABg(:,1),la_53,fi_53,'linear')
 APP_53_i2=griddata(la(:,1),fi(:,1),APP(:,1),la_53,fi_53,'linear')

 
 % measured gravity anomalies at test point
 Afa_53=g_53+0.30855*H_53-(978032.66*(1+0.0053024*((sin(fi_53*pi/180))^2-0.00000585*((sin(2*fi_53*pi/180))^2))))
 ABg_53=g_53+(-0.04192*sigma*H_53)+(0.30855*H_53)-(978032.66*(1+0.0053024*((sin(fi_53*pi/180))^2-0.00000585*((sin(2*fi_53*pi/180))^2))))
 APP_53=g_53+((0.30888-0.08384*sigma)*H_53)-(978032.66*(1+0.0053024*((sin(fi_53*pi/180))^2-0.00000585*((sin(2*fi_53*pi/180))^2))))
 
 % residuals of two method interpolation
 resfa1_53=Afa_53-Afa_53_i1
 resfa2_53=Afa_53-Afa_53_i2

 resB1_53=ABg_53-ABg_53_i1
 resB2_53=ABg_53-ABg_53_i2

 resPP1_53=APP_53-APP_53_i1
 resPP2_53=APP_53-APP_53_i2

 % interpolation methods
 % sample point number 103 in the test_points
 
 fi_103=54.70619
 la_103=17.72920
 H_103=120.964
 g_103=981448.392



 % 1. using of empirical H relationship
 Afa_103_i1=p1(1,1)*H_103+p1(1,2)
 ABg_103_i1=p2(1,1)*H_103+p2(1,2)
 APP_103_i1=p3(1,1)*H_103+p3(1,2)

 
 

 
 % 2 spatial interpolation from phi and lambda
 Afa_103_i2=griddata(la(:,1),fi(:,1),Afa(:,1),la_103,fi_103,'linear')
 ABg_103_i2=griddata(la(:,1),fi(:,1),ABg(:,1),la_103,fi_103,'linear')
 APP_103_i2=griddata(la(:,1),fi(:,1),APP(:,1),la_103,fi_103,'linear')

 
 % measured gravity anomalies at test point
 Afa_103=g_103+0.30855*H_103-(978032.66*(1+0.0053024*((sin(fi_103*pi/180))^2-0.00000585*((sin(2*fi_103*pi/180))^2))))
 ABg_103=g_103+(-0.04192*sigma*H_103)+(0.30855*H_103)-(978032.66*(1+0.0053024*((sin(fi_103*pi/180))^2-0.00000585*((sin(2*fi_103*pi/180))^2))))
 APP_103=g_103+((0.30888-0.08384*sigma)*H_103)-(978032.66*(1+0.0053024*((sin(fi_103*pi/180))^2-0.00000585*((sin(2*fi_103*pi/180))^2))))
 
 % residuals of two method interpolation
 resfa1_103=Afa_103-Afa_103_i1
 resfa2_103=Afa_103-Afa_103_i2

 resB1_103=ABg_103-ABg_103_i1
 resB2_103=ABg_103-ABg_103_i2

 resPP1_103 =APP_103-APP_103_i1
 resPP2_103 =APP_103-APP_103_i2

  % interpolation methods
 % sample point number 103 in the test_points
 
 fi_183=51.22722
 la_183=16.04245
 H_183=162.618
 g_183=981141.500



 % 1. using of empirical H relationship
 Afa_183_i1=p1(1,1)*H_183+p1(1,2)
 ABg_183_i1=p2(1,1)*H_183+p2(1,2)
 APP_183_i1=p3(1,1)*H_183+p3(1,2)

 
 

 
 % 2 spatial interpolation from phi and lambda
 Afa_183_i2=griddata(la(:,1),fi(:,1),Afa(:,1),la_183,fi_183,'linear')
 ABg_183_i2=griddata(la(:,1),fi(:,1),ABg(:,1),la_183,fi_183,'linear')
 APP_183_i2=griddata(la(:,1),fi(:,1),APP(:,1),la_183,fi_183,'linear')

 
 % measured gravity anomalies at test point
 Afa_183=g_183+0.30855*H_183-(978032.66*(1+0.0053024*((sin(fi_183*pi/180))^2-0.00000585*((sin(2*fi_183*pi/180))^2))))
 ABg_183=g_183+(-0.04192*sigma*H_183)+(0.30855*H_183)-(978032.66*(1+0.0053024*((sin(fi_183*pi/180))^2-0.00000585*((sin(2*fi_183*pi/180))^2))))
 APP_183=g_183+((0.30888-0.08384*sigma)*H_183)-(978032.66*(1+0.0053024*((sin(fi_183*pi/180))^2-0.00000585*((sin(2*fi_183*pi/180))^2))))
 
 % residuals of two method interpolation
 resfa1_183=Afa_183-Afa_183_i1
 resfa2_183=Afa_183-Afa_183_i2

 resB1_183=ABg_183-ABg_183_i1
 resB2_183=ABg_183-ABg_183_i2

 resPP1_183 =APP_183-APP_183_i1
 resPP2_183 =APP_183-APP_183_i2


 % interpolation methods
 % sample point number 3 in the test_points
 
 fi_233=51.49563
 la_233=18.52562
 H_233=190.169
 g_233=981154.398



 % 1. using of empirical H relationship
 Afa_233_i1=p1(1,1)*H_233+p1(1,2)
 ABg_233_i1=p2(1,1)*H_233+p2(1,2)
 APP_233_i1=p3(1,1)*H_233+p3(1,2)

 
 

 
 % 2 spatial interpolation from phi and lambda
 Afa_233_i2=griddata(la(:,1),fi(:,1),Afa(:,1),la_233,fi_233,'linear')
 ABg_233_i2=griddata(la(:,1),fi(:,1),ABg(:,1),la_233,fi_233,'linear')
 APP_233_i2=griddata(la(:,1),fi(:,1),APP(:,1),la_233,fi_233,'linear')

 
 % measured gravity anomalies at test point
 Afa_233=g_233+0.30855*H_233-(978032.66*(1+0.0053024*((sin(fi_233*pi/180))^2-0.00000585*((sin(2*fi_233*pi/180))^2))))
 ABg_233=g_233+(-0.04192*sigma*H_233)+(0.30855*H_233)-(978032.66*(1+0.0053024*((sin(fi_233*pi/180))^2-0.00000585*((sin(2*fi_233*pi/180))^2))))
 APP_233=g_233+((0.30888-0.08384*sigma)*H_233)-(978032.66*(1+0.0053024*((sin(fi_233*pi/180))^2-0.00000585*((sin(2*fi_233*pi/180))^2))))
 
 % residuals of two method interpolation
 resfa1_233=Afa_233-Afa_233_i1
 resfa2_233=Afa_233-Afa_233_i2

 resB1_233=ABg_233-ABg_233_i1
 resB2_233=ABg_233-ABg_233_i2

 resPP1_233=APP_233-APP_233_i1
 resPP2_233=APP_233-APP_233_i2

figure('Name','plots')
    x = [3 53 103 183 233];
    y = [resfa1_3 resfa2_3 resB1_3 resB2_3 resPP1_3 resPP2_3
        resfa1_53 resfa2_53 resB1_53 resB2_53 resPP1_53 resPP2_53
        resfa1_103 resfa2_103 resB1_103 resB2_103 resPP1_103 resPP2_103
        resfa1_183 resfa2_183 resB1_183 resB2_183 resPP1_183 resPP2_183
        resfa1_233 resfa2_233 resB1_233 resB2_233 resPP1_233 resPP2_233];
    b = bar(x,y);
    xtips2 = b(2).XEndPoints;
    ytips2 = b(2).YEndPoints;
    labels2 = string(b(2).YData);
    text(xtips2,ytips2,labels2,'HorizontalAlignment','center',...
        'VerticalAlignment','bottom')

%%
