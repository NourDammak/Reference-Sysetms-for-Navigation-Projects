format long
clear

%upload time series from file with Txyz

ts=load('Budapest_Time_Serie_Data.txt');
time=ts(:,1);
x=ts(:,2);
y=ts(:,3);
z=ts(:,4);

% visualisation of time series xyz separately,

subplot(3,3,1)
plot(time,x,'o')
title('change of X')
xlabel('Time [year]')
ylabel('x [m]')

subplot(3,3,2)
plot(time,y,'o')
title('change of Y')
xlabel('Time [year]')
ylabel('y [m]')

subplot(3,3,3)
plot(time,z,'o')
title('change of Z')
xlabel('Time [year]')
ylabel('z [m]')

% linear trend fit (regresion), velocity determination..

subplot(3,3,4)
ax=polyfit(time,x,1);
xreg=polyval(ax,time);
plot(time,x,'o',time,xreg,'-')
title([' X and vx = ',num2str(ax(1,1),'%.4f'),' m/year'])
xlabel('Time [year]')
ylabel('x [m]')

subplot(3,3,5)
ay=polyfit(time,y,1);
yreg=polyval(ay,time);
plot(time,y,'o',time,yreg,'-')
title([' y and vy = ',num2str(ay(1,1),'%.4f'),' m/year'])
xlabel('Time [year]')
ylabel('y [m]')


subplot(3,3,6)
az=polyfit(time,z,1);
zreg=polyval(az,time);
plot(time,z,'o',time,zreg,'-')
title([' z and vz = ',num2str(az(1,1),'%.4f'),' m/year'])
xlabel('Time [year]')
ylabel('z [m]')


%time for bar plot with residuals, differences between observed coordinate
%and linear model
subplot(3,3,7)
res_x=x-polyval(ax,time);
bar(time,res_x*1000)
title('residual in mm')
xlabel('Time [year]')
ylabel('x res [mm]')

subplot(3, 3, 8)
res_y=y-polyval(ay,time);
bar(time,res_y*1000)
title('residual in mm')
xlabel('Time [year]')
ylabel('y res [mm]')


subplot(3, 3, 9)
res_z=z-polyval(az,time);
bar(time,res_z*1000)
title('residual in mm')
xlabel('Time [year]')
ylabel('z res [mm]')

