%Script that calculates heat content (or heat storage)
%Fernandez August2015

clear all
clc
%This part calculates OHC using cp, rho and DTdz
load('ARGO_thickness_STMW_new.mat','DTdzrange2','lona','lata','time','pressure','mld');
load('monthly_ARGO_STMW_new.mat','T','S');

%make a control volume 
% ii = min(nearest(lata,-20));
% ff = min(nearest(lata,-40));
% iii = min(nearest(lona,150));
% fff = min(nearest(lona,180));

ii = min(nearest(lata,-28));
ff = min(nearest(lata,-34));
iii = min(nearest(lona,173));
fff = min(nearest(lona,180));

T = T(iii:fff,ff:ii,:,:);
S = S(iii:fff,ff:ii,:,:);
lon = lona(iii:fff);
lat = lata(ff:ii);

cp = nan*T;
for t = 1:length(time)
for j = 1:length(lon)
    for i = 1:length(lat)
cp(j,i,:,t) =  sw_cp(squeeze(S(j,i,:,t)),squeeze(T(j,i,:,t)),pressure);
    end 
end
end

% rho = 1025;
rho_z = nan*T;
for t = 1:length(time)
for j = 1:length(lon)
    for i = 1:length(lat)
        for p = 1:length(pressure)
%H(j,i,p,t) =  rho*squeeze(cp(j,i,p,t))*squeeze(T(j,i,p,t));%Units [Joules]
rho_z(j,i,p,t) = sw_pden(S(j,i,p,t),T(j,i,p,t),pressure(p),0);
        end  
    end 
end
end

%find mean depth of isotherm of 14C
ind14 = nan(length(lon),length(lat),length(time));
z14 = nan*ind14;
%find mean depth of isotherm of 20C

ind20 = nan(length(lon),length(lat),length(time));
z20 = nan*ind20;

for t = 1:length(time)
for j = 1:length(lon)
    for i = 1:length(lat)
        if isnan(T(j,i,:,t))
            ind14(j,i,t) = nan;
            ind20(j,i,t) = nan;
        else
    ind14(j,i,t) = min(nearest(T(j,i,:,t),14));
    z14(j,i,t) = pressure(squeeze(ind14(j,i,t)));
    ind20(j,i,t) = min(nearest(T(j,i,:,t),20));
    z20(j,i,t) = pressure(squeeze(ind20(j,i,t)));

        end
    end
end
end

%calculate the average depth

tmp14 = reshape(z14,length(lon)*length(lat)*length(time),1);
meanz14 = nan_mean(tmp14);
ind_meanz14 = nearest(pressure,meanz14);
%ind_meanz14 = 20;

tmp20 = reshape(z20,length(lon)*length(lat)*length(time),1);
meanz20 = nan_mean(tmp20);
ind_meanz20 = nearest(pressure,meanz20);

%OHC is sum(rho*cp*(dT(z)/dt)*(z14-z20)) 

OHC = nan*ind14;
meanT = nan_mean(T,4);

for t = 1:length(time)
for j = 1:length(lon)
    for i = 1:length(lat)
        if isnan(T(j,i,:,t))
            OHC(j,i,t) = nan;
        else
    OHC(j,i,t) = sum(squeeze(squeeze(rho_z(j,i,ind_meanz14:-1:ind_meanz20,t))).*squeeze(squeeze(cp(j,i,ind_meanz14:-1:ind_meanz20)))...
        .*squeeze(squeeze(T(j,i,ind_meanz14:-1:ind_meanz20,t))).*(-diff(pressure(ind_meanz14+1:-1:ind_meanz20))));%minus sign is to revert sign of integration
        end
    end
end
end
% for t = 1:length(time)
% for j = 1:length(lon)
%     for i = 1:length(lat)
%         if isnan(T(j,i,:,t))
%             OHC(j,i,t) = nan;
%         else
%     OHC(j,i,t) = sum(squeeze(squeeze(rho_z(j,i,ind_meanz14:-1:ind_meanz20,t))).*squeeze(squeeze(cp(j,i,ind_meanz14:-1:ind_meanz20)))...
%         .*squeeze(squeeze(T(j,i,ind_meanz14:-1:ind_meanz20,t)-meanT(j,i,ind_meanz14:-1:ind_meanz20))).*(-diff(pressure(ind_meanz14+1:-1:ind_meanz20))));%minus sign is to revert sign of integration
%         end
%     end
% end
% end

%do average OHC along PX06 line

lon176 = nearest(lon,176.5);
ts_OHC_ARGO = nan_mean(squeeze(OHC(lon176,:,:)),1);
ts_OHC_ARGO = ts_OHC_ARGO - mean(ts_OHC_ARGO);
%deseason 

[yy, mm, tmp, tmp, tmp, tmp]=datevec(time);

for ii=1:12
  ind=find(mm==ii);
  OHC_ave_ARGO(ii) = mean(ts_OHC_ARGO(ind));
  OHC_std_ARGO(ii) = std(ts_OHC_ARGO(ind));
  OHC_ste_ARGO(ii) = OHC_std_ARGO(ii)./sqrt(11);
  
end;
%deseason OHC
for ii=1:length(time)
  [yy,mm,tmp,tmp,tmp,tmp]=datevec(time(ii));
  ts_OHC_ARGO_ds(ii) = ts_OHC_ARGO(ii) - OHC_ave_ARGO(mm);
end;

save OHC_ARGO OHC ts_OHC_ARGO ts_OHC_ARGO_ds OHC_ave_ARGO OHC_std_ARGO OHC_ste_ARGO lon lat time pressure cp rho_z 

%load OHC from HRX and compare with Argo
load('TS_inventory_STMW_RG_Argo_new.mat', 'TS_thermostad_Argo')
load('../EAUC_XBT/PX06/thermostad_PXO6.mat','OHC_TS','time_thermostad_TS');
ts_OHC_HRX = OHC_TS;
figure(4)
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 20 8]);
[h1] = plot(time,10^(-8)*cosfilt(ts_OHC_ARGO,3,1),'color',[0.2 0.8 0.8],'LineWidth',2);hold on
[h1b] = plot(time_thermostad_TS,10^(-8)*(ts_OHC_HRX-mean(ts_OHC_HRX)),'color',[0 0 0],'LineWidth',2);hold on
plot(time_thermostad_TS,0,'color',[0 0 0],'LineWidth',2);hold on
ax1 = gca;
set(ax1,'XColor','k','YColor','k','XTickLabel',1986:2:2014,...
           'XTick',[datenum([1986:2:2014],1,1)],'FontSize',12)
xlim([725000 736000]);
ylim([-15 20]);
ylabel('OHC [10^8 J m^{-2}]','FontSize',12); 
legend([h1,h1b],'RG OI','HRX','location','north','orientation','horizontal');
%print -depsc E:\PhD\PhD_writings\PhD_Thesis_final\figures\timeseries_OHC_ARGO_HRX.eps

%calculate correlation between STMW and OHC
x = cosfilt(TS_thermostad_Argo,3,1);
y = cosfilt(ts_OHC_ARGO,3,1);
[rho] = corrcoef(x(5:end-4),y(5:end-4));

%deseason OHC HRX 
[yy, mm, tmp, tmp, tmp, tmp]=datevec(time_thermostad_TS);

for ii=1:12
  ind=find(mm==ii);
  OHC_HRX_av(ii) = mean(ts_OHC_HRX(ind));
  OHC_HRX_std(ii) = std(ts_OHC_HRX(ind));
  OHC_HRX_ste(ii) = OHC_HRX_std(ii)./sqrt(length(ind));
  
end;
for ii=1:length(time_thermostad_TS)
  [yy,mm,tmp,tmp,tmp,tmp]=datevec(time_thermostad_TS(ii));
  OHC_HRX_ds(ii) = ts_OHC_HRX(ii) - OHC_HRX_av(mm);
end;

load ../lsha_updated/OHC_LSHA.mat
OHC_LSHA = -OHC_LSHA;
ts_OHC_LSHA = mean(OHC_LSHA,1);
figure(44)
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 20 8]);
[h111b] = plot(time_thermostad_TS,10^(-8)*(ts_OHC_HRX-mean(ts_OHC_HRX)),'color',[0 0 0],'LineWidth',2);hold on
[h1111b] = plot(timeh,10^(-8)*(ts_OHC_LSHA-mean(ts_OHC_LSHA)),'color',[1 0 1],'LineWidth',2);hold on
plot(time_thermostad_TS,0,'color',[0 0 0],'LineWidth',2);hold on
ax1 = gca;
set(ax1,'XColor','k','YColor','k','XTickLabel',1986:2:2014,...
           'XTick',[datenum([1986:2:2014],1,1)],'FontSize',12)
xlim([725000 736000]);
ylim([-15 20]);
ylabel('OHC [10^8 J m^{-2}]','FontSize',12); 
legend([h111b,h1111b],'HRX','LSHA','location','north','orientation','horizontal');
%print -depsc E:\PhD\PhD_writings\PhD_Thesis_final\figures\timeseries_OHC_ARGO_HRX.eps

%Compare with surface heat fluxes from NCEP

load ../NCEP_data/NCEP_total_hflx.mat
%make a control volume 

ii = min(nearest(latt,-29));
ff = min(nearest(latt,-34));
iii = min(nearest(lont,173));
fff = min(nearest(lont,180));
tind = nearest(time_hflx,datenum(2004,01,01));
tfin = nearest(time_hflx,datenum(2014,12,01));
time2 = time_hflx(tind:tfin);
total_htflx_NZ = -total_hflx(iii:fff,ii:ff,tind:tfin);%negative sign is for the ocean
lon = lont(iii:fff);
lat = latt(ii:ff);

area_htflx = reshape(total_htflx_NZ,length(lon)*length(lat),length(time2),1);
mean_htflx = nan_mean(area_htflx,1);

load ../JRA-55/JRA-55_netHF.mat
%make a control volume 

ii = max(nearest(lattj,-29.5));
ff = min(nearest(lattj,-34));
iii = min(nearest(lontj,173));
fff = min(nearest(lontj,180));
tind = nearest(timej,datenum(2004,01,15));
tfin = nearest(timej,datenum(2014,12,01));
time3 = timej(tind:tfin);
total_htflx_NZ_JRA = netHF(iii:fff,ff:ii,tind:tfin);%negative sign is for the ocean
lon = lontj(iii:fff);
lat = lattj(ff:ii);

area_htflx_JRA = reshape(total_htflx_NZ_JRA,length(lon)*length(lat),length(time3),1);
mean_htflx_JRA = nan_mean(area_htflx_JRA,1);

figure(5)
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 20 8]);
%[h1] = plot(time,10^(-8)*cosfilt(ts_OHC,7,1),'color',[0.3 0.3 0.3],'LineWidth',2);hold on
[h1] = plot(time3,cosfilt(mean_htflx_JRA,7,1),'color',[0.8 0.3 0.8],'LineWidth',2);hold on

plot(time,0,'color',[0 0 0],'LineWidth',2);hold on
ylim([-150 100]);
ax1 = gca;
set(ax1,'XColor','w','YColor','k','XTickLabel',2004:1:2014,...
           'XTick',[datenum([2004:1:2014],1,1)],'FontSize',12)
xlim([731900 735600]);
ylabel('Surface Heat Flux [W m^{-2}]','FontSize',12); 
set(gca,'Xtick',[]);
ax2 = axes('Position',get(ax1,'Position'),...
           'XTickLabel',2004:1:2014,...
           'XTick',[datenum([2004:1:2014],1,1)],...
           'YAxisLocation','right',...
           'Color','none',...
           'XColor','k','YColor','k','FontSize',12);
hold(ax2,'on')
[h2] = plot(time,cosfilt(TS_thermostad_Argo,7,1),'color',[0.2 0.8 0.8],'LineWidth',3);
xlim([731900 735600]);
ylim([0 250]);

legend([h1,h2],'Surface Heat flux','STMW','location','north','orientation','horizontal');
ylabel('STMW inventory [km^2]','FontSize',12);  

x = cosfilt(TS_thermostad_Argo,7,1);
y = cosfilt(mean_htflx_JRA,7,1);
[rho] = corrcoef(x(5:end-13),y(5:end-4));

%print -depsc E:\PhD\PhD_writings\PhD_Thesis_final\figures\HTFLX_STMW.eps

%deseason ARGO data
[yy, mm, tmp, tmp, tmp, tmp]=datevec(time);


for ii=1:12
  ind=find(mm==ii);
  thermostadav_ARGO(ii) = mean(TS_thermostad_Argo(ind));
  thermostadstd_ARGO(ii) = std(TS_thermostad_Argo(ind));
  thermostadste_ARGO(ii) = thermostadstd_ARGO(ii)./sqrt(11);
end;

for ii=1:length(time)
  [yy,mm,tmp,tmp,tmp,tmp]=datevec(time(ii));
  TS_thermostad_Argo_ds(ii) = TS_thermostad_Argo(ii) - thermostadav_ARGO(mm);
end; 

%deseason mean heat flux and STMW
[yy, mm, tmp, tmp, tmp, tmp]=datevec(time2);

for ii=1:12
  ind=find(mm==ii);
  htflx_ave(ii) = mean(mean_htflx(ind));
  htflx_std(ii) = std(mean_htflx(ind));
  htflx_ste(ii) = htflx_std(ii)./sqrt(length(ind));
  
end;
%deseason surface heat flux 
for ii=1:length(time2)
  [yy,mm,tmp,tmp,tmp,tmp]=datevec(time2(ii));
  ts_htflx_ds(ii) = mean_htflx(ii) - htflx_ave(mm);
end;

%deseason mean heat flux and STMW JRA-55
[yy, mm, tmp, tmp, tmp, tmp]=datevec(time3);

for ii=1:12
  ind=find(mm==ii);
  htflx_ave_JRA(ii) = mean(mean_htflx_JRA(ind));
  htflx_std_JRA(ii) = std(mean_htflx_JRA(ind));
  htflx_ste_JRA(ii) = htflx_std_JRA(ii)./sqrt(length(ind));
  
end;
%deseason surface heat flux JRA-55
for ii=1:length(time3)
  [yy,mm,tmp,tmp,tmp,tmp]=datevec(time3(ii));
  ts_htflx_JRA_ds(ii) = mean_htflx_JRA(ii) - htflx_ave_JRA(mm);
end;

figure(6)
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 20 8]);
%[h1] = plot(time,10^(-8)*cosfilt(ts_OHC,7,1),'color',[0.3 0.3 0.3],'LineWidth',2);hold on
[h1] = plot(time3,cosfilt(ts_htflx_JRA_ds,7,1),'color',[0.8 0.3 0.8],'LineWidth',2);hold on
%[h1b] = plot(time3,cosfilt(ts_htflx_JRA_ds,7,1),'color',[0 0 0.8],'LineWidth',2);hold on

plot(time,0,'color',[0 0 0],'LineWidth',2);hold on
ylim([-25 25]);
ax1 = gca;
set(ax1,'XColor','w','YColor','k','XTickLabel',2004:1:2014,...
           'XTick',[datenum([2004:1:2014],1,1)],'FontSize',12)
xlim([731900 735600]);
ylabel('Surface Heat Flux [W m^{-2}]','FontSize',12); 
set(gca,'Xtick',[]);
ax2 = axes('Position',get(ax1,'Position'),...
           'XTickLabel',2004:1:2014,...
           'XTick',[datenum([2004:1:2014],1,1)],...
           'YAxisLocation','right',...
           'Color','none',...
           'XColor','k','YColor','k','FontSize',12);
hold(ax2,'on')
[h2] = plot(time,cosfilt(TS_thermostad_Argo_ds,7,1),'color',[0.2 0.8 0.8],'LineWidth',3);
xlim([731900 735600]);
ylim([-50 50]);

legend([h1,h2],'Surface Heat flux','STMW','location','north','orientation','horizontal');
ylabel('STMW inventory [km^2]','FontSize',12);  

x = cosfilt(TS_thermostad_Argo,7,1);
y = cosfilt(ts_htflx_JRA_ds,7,1);
[rho] = corrcoef(x(5:end-13),y(5:end-4));

%% Scatter of OHC in February vs MLD in August
load('ARGO_thickness_STMW_new.mat','lona','mld');

lats = min(nearest(lata,-34));
latn = min(nearest(lata,-28));
lon177 = min(nearest(lona,176.5));


mld_EAUC = squeeze(mld(lon177,lats:latn,:));
mld_EAUC_ave = nan_mean(mld_EAUC,1);

[yy, mm, tmp, tmp, tmp, tmp]=datevec(time);
% ind_w = find(mm>=8 & mm<=9);
% ind_s = find(mm>=1 & mm<=2);
ind_w = find(mm==8);
ind_s = find(mm==2);

figure;
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0 0 22 8]);
h1 = scatter(mld_EAUC_ave(ind_w),10^(-8)*ts_OHC_ARGO(ind_s),50);hold on
%h1 = scatter(mld_EAUC_ave,10^(-8)*ts_OHC,50);hold on

xlabel('Later winter MLD [m]');
ylabel('Summer OHC [10^8 J m^{-2}]');grid;
y1=10^(-8)*ts_OHC_ARGO(ind_s)';
x=mld_EAUC_ave(ind_w)';
hold on
box('on');
% %Try adding variance ellipse for fitting the best line that minimise
% %variance in both x and y
% 
% [M]=[x,y1];
% 
% %remove NaN by creating new M with non NaNs
% idx=find(~isnan(sum(M.')));
% u=x(idx);
% v=y1(idx);
% 
% %Demean
% m_u=mean(u);
% m_v=mean(v);
% nomean_u=u-m_u;
% nomean_v=v-m_v;
% [F]=[nomean_u,nomean_v];
% %Calculate covariance matrix R
% [R]=F'*F;
% %calculate eigenvalues and eigenvectors, eigenvalues are major, minor
% %ellipse axix
% [C,L]=eig(R);
% %plotting error ellipse
% plot(u,v,'gx');hold on
% draw_variance_ellipse(m_u,m_v,u,v,0.3);
% axis([60 160 -5 5]);

%print -depsc E:\PhD\PhD_writings\PhD_Thesis_final\figures\scatter_MLDwinter_OHCsummer_ellipse.eps

%compare with N2
load('monthly_ARGO_STMW_new.mat', 'N2','lona','lata','time','p_ave')

lon177 = min(nearest(lona,177));
ii = min(nearest(lata,-34));
ff = min(nearest(lata,-28));
pressure = squeeze(p_ave(lon177,ff,:,1));
iii =  ind_meanz14;
%iii = min(nearest(pressure,200));
pressure(iii)
N2_tmp1 = squeeze(nan_mean(squeeze(N2(lon177,ii:ff,1:iii,:)),1));

N2_tmp2 = nan_mean(N2_tmp1,1);

%deseason ARGO data
[yy, mm, tmp, tmp, tmp, tmp]=datevec(time);


for ii=1:12
  ind=find(mm==ii);
  N2av_ARGO(ii) = mean(N2_tmp2(ind));
  N2std_ARGO(ii) = std(N2_tmp2(ind));
  N2ste_ARGO(ii) = N2std_ARGO(ii)./sqrt(11);
end;

for ii=1:length(time)
  [yy,mm,tmp,tmp,tmp,tmp]=datevec(time(ii));
  N2_Argo_ds(ii) = N2_tmp2(ii) - N2av_ARGO(mm);
end; 

[yy, mm, tmp, tmp, tmp, tmp]=datevec(time);
% ind_w = find(mm>=8 & mm<=9);
% ind_s = find(mm>=2 & mm<=3);
ind_w = find(mm==8);
ind_s = find(mm==2);

figure;
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0 0 22 8]);
%h1 = scatter(sqrt(N2_tmp2(ind_s(1:end-2))),mld_EAUC_ave(ind_w(3:end)));hold on
h2 = scatter(10^(3)*sqrt(N2_tmp2(ind_s)),mld_EAUC_ave(ind_w));hold on

%h1 = scatter(10^(3)*sqrt(N2_tmp2),mld_EAUC_ave);hold on

ylabel('Later winter MLD [m]');
xlabel('N [10^{-3} s^{-1}]');grid;
y1=mld_EAUC_ave(ind_w)';
x=10^(3)*sqrt(N2_tmp2(ind_s))';
hold on
box('on');
% %Try adding variance ellipse for fitting the best line that minimise
% %variance in both x and y
% 
% [M]=[x,y1];
% 
% %remove NaN by creating new M with non NaNs
% idx=find(~isnan(sum(M.')));
% u=x(idx);
% v=y1(idx);
% 
% %Demean
% m_u=mean(u);
% m_v=mean(v);
% nomean_u=u-m_u;
% nomean_v=v-m_v;
% [F]=[nomean_u,nomean_v];
% %Calculate covariance matrix R
% [R]=F'*F;
% %calculate eigenvalues and eigenvectors, eigenvalues are major, minor
% %ellipse axix
% [C,L]=eig(R);
% %plotting error ellipse
% plot(u,v,'gx');hold on
% %draw_variance_ellipse(m_u,m_v,u,v,0.5);
% error_ellipse(R,[m_u m_v]);
% %axis([60 160 -5 5]);

%print -depsc E:\PhD\PhD_writings\PhD_Thesis_final\figures\scatter_MLDwinter_Nsummer.eps
