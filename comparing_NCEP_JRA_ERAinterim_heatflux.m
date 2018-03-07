%comparing NCEP and JRA-55 heat fluxes

%December 2015
clear all
clc

load ../NCEP_data/NCEP_total_hflx.mat
%make a control volume 

ii = min(nearest(latt,-29));
ff = min(nearest(latt,-34));
iii = min(nearest(lont,173));
fff = min(nearest(lont,180));
tind = nearest(time_hflx,datenum(1980,01,01));
tfin = nearest(time_hflx,datenum(2015,01,01));
time2 = time_hflx(tind:tfin);
total_htflx_NZ = -total_hflx(iii:fff,ii:ff,tind:tfin);%negative sign is for the ocean
lon = lont(iii:fff);
lat = latt(ii:ff);

area_htflx = reshape(total_htflx_NZ,length(lon)*length(lat),length(time2),1);
mean_htflx = nan_mean(area_htflx,1);

load ../JRA-55/JRA-55_netHF.mat
%make a control volume 

ii = min(nearest(lattj,-29));
ff = min(nearest(lattj,-34));
iii = min(nearest(lontj,173));
fff = min(nearest(lontj,180));
tind = nearest(timej,datenum(1980,01,01));
tfin = nearest(timej,datenum(2015,01,01));
time2j = timej(tind:tfin);
total_htflx_NZ_JRA = netHF(iii:fff,ff:ii,tind:tfin);%negative sign is for the ocean
lon = lontj(iii:fff);
lat = lattj(ff:ii);

area_htflx_JRA = reshape(total_htflx_NZ_JRA,length(lon)*length(lat),length(time2j),1);
mean_htflx_JRA = nan_mean(area_htflx_JRA,1);


load ./ERA_interim_hflx.mat
%make a control volume 

ii = min(nearest(lat,-29));
ff = min(nearest(lat,-34));
iii = min(nearest(lon,173));
fff = min(nearest(lon,180));
tind = nearest(time_ERA,datenum(1980,01,01));
tfin = nearest(time_ERA,datenum(2015,01,01));
time2e = time_ERA(tind:tfin);
total_htflx_NZ_ERA = total_hflx_ERA(iii:fff,ii:ff,tind:tfin);%negative sign is for the ocean
lone = lon(iii:fff);
late = lat(ii:ff);

area_htflx_ERA = reshape(total_htflx_NZ_ERA,length(lone)*length(late),length(time2e),1);
mean_htflx_ERA = nan_mean(area_htflx_ERA,1);


figure(55)
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 20 8]);
%[h1] = plot(time,10^(-8)*cosfilt(ts_OHC,7,1),'color',[0.3 0.3 0.3],'LineWidth',2);hold on
[h1] = plot(time2,cosfilt(mean_htflx,7,1),'color',[1 0 0],'LineWidth',2);hold on
[h1b] = plot(time2j,cosfilt(mean_htflx_JRA,7,1),'color',[0 0 0],'LineWidth',1);hold on
[h1c] = plot(time2e,cosfilt(mean_htflx_ERA,7,1),'color',[0 0 1],'LineWidth',1);hold on

plot(time2,0,'color',[0 0 0],'LineWidth',2);hold on
%ylim([120 260]);
ax1 = gca;
set(ax1,'XColor','k','YColor','k','XTickLabel',1980:10:2015,...
           'XTick',[datenum([1980:10:2015],1,1)],'FontSize',12)
%xlim([731900 735600]);
ylabel('Surface Heat Flux [W m^{-2}]','FontSize',12); 
legend([h1 h1b h1c],'NCEP','JRA-55','ERA Interim');
%[rho,p] = corrcoef(cosfilt(mean_htflx,1,1),cosfilt(mean_htflx_JRA,1,1))

%deseason mean heat flux NCEP
[yy, mm, tmp, tmp, tmp, tmp]=datevec(time2);

for ii=1:12
  ind=find(mm==ii);
  htflx_ave(ii) = mean(mean_htflx(ind));
  htflx_std(ii) = std(mean_htflx(ind));
  htflx_ste(ii) = htflx_std(ii)./sqrt(length(ind));
  
end;
%deseason surface heat flux NCEP
for ii=1:length(time2)
  [yy,mm,tmp,tmp,tmp,tmp]=datevec(time2(ii));
  ts_htflx_ds(ii) = mean_htflx(ii) - htflx_ave(mm);
end;

%deseason mean heat flux JRA-55
[yy, mm, tmp, tmp, tmp, tmp]=datevec(time2j);

for ii=1:12
  ind=find(mm==ii);
  htflx_ave_JRA(ii) = mean(mean_htflx_JRA(ind));
  htflx_std_JRA(ii) = std(mean_htflx_JRA(ind));
  htflx_ste_JRA(ii) = htflx_std_JRA(ii)./sqrt(length(ind));
  
end;
%deseason surface heat flux JRA-55
for ii=1:length(time2j)
  [yy,mm,tmp,tmp,tmp,tmp]=datevec(time2j(ii));
  ts_htflx_JRA_ds(ii) = mean_htflx_JRA(ii) - htflx_ave_JRA(mm);
end;

%deseason mean heat flux and ERA interim
[yy, mm, tmp, tmp, tmp, tmp]=datevec(time2e);

for ii=1:12
  ind=find(mm==ii);
  htflx_ave_ERA(ii) = mean(mean_htflx_ERA(ind));
  htflx_std_ERA(ii) = std(mean_htflx_ERA(ind));
  htflx_ste_ERA(ii) = htflx_std_ERA(ii)./sqrt(length(ind));
  
end;
%deseason surface heat flux ERA
for ii=1:length(time2e)
  [yy,mm,tmp,tmp,tmp,tmp]=datevec(time2e(ii));
  ts_htflx_ERA_ds(ii) = mean_htflx_ERA(ii) - htflx_ave_ERA(mm);
end;

figure(56)
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 20 8]);

[h1] = plot(time2,cosfilt(ts_htflx_ds,7,1),'color',[1 0 0],'LineWidth',2);hold on
[h1b] = plot(time2j,cosfilt(ts_htflx_JRA_ds,7,1),'color',[0 0 0],'LineWidth',2);hold on
[h1c] = plot(time2e,cosfilt(ts_htflx_ERA_ds,7,1),'color',[0 1 0],'LineWidth',2);hold on

plot(time2,0,'color',[0 0 0],'LineWidth',2);hold on
%ylim([120 260]);
ax1 = gca;
set(ax1,'XColor','k','YColor','k','XTickLabel',1980:10:2015,...
           'XTick',[datenum([1980:10:2015],1,1)],'FontSize',12)
%xlim([731900 735600]);
ylabel('Surface Heat Flux [W m^{-2}]','FontSize',12); 
legend([h1 h1b h1c],'NCEP','JRA-55','ERA-Interim');
%[rho,p] = corrcoef(cosfilt(ts_htflx_ds,1,1),cosfilt(ts_htflx_JRA_ds,1,1))
% %% Differences
% 
% figure(57)
% set(gcf, 'PaperPositionMode', 'manual');
% set(gcf, 'PaperUnits', 'centimeters');
% set(gcf, 'PaperPosition', [0 0 20 8]);
% 
% [h1] = plot(time2,cosfilt(ts_htflx_ds - ts_htflx_JRA_ds,1,1),'color',[1 0 0],'LineWidth',2);hold on
% [h1b] = plot(time2j,cosfilt(ts_htflx_JRA_ds - ts_htflx_ERA_ds,1,1),'color',[0 0 0],'LineWidth',2);hold on
% [h1c] = plot(time2e,cosfilt(ts_htflx_ds - ts_htflx_ERA_ds,1,1),'color',[0 1 0],'LineWidth',2);hold on
% 
% plot(time2,0,'color',[0 0 0],'LineWidth',2);hold on
% %ylim([120 260]);
% ax1 = gca;
% set(ax1,'XColor','k','YColor','k','XTickLabel',1980:1:2015,...
%            'XTick',[datenum([1980:1:2015],1,1)],'FontSize',12)
% %xlim([731900 735600]);
% ylabel('Surface Heat Flux [W m^{-2}]','FontSize',12); 
% legend([h1 h1b h1c],'NCEP - JRA ','JRA-55 - ERA-Interim','NCEP - ERA-Interim');