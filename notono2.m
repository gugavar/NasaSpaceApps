%% Load Data and latitude/longtigude vectors
home;
clear;
clear all;
NO2 = open('OMI2019.mat');
Cdata = open('CData.mat');
NO2=NO2.OMI2019;
Cdata= Cdata.Cdata;
lat = [-89.875:0.250:89.875];
long= [-179.875:0.250:179.875];

%% WORLD
figure(1)
latlim = [-89.875 89.875];
lonlim = [-179.875 179.875];
worldmap(latlim,lonlim) % create worldmap
load coastlines
hold on
pcolorm(lat,long,NO2) % plot NO2 data
hold on
geoshow(coastlat,coastlon,'Color','k') %show coastlines
colorbar
colormap jet
title('Time Averaged (2019) Map of NO_{2} Tropospheric Column (1/cm^2)')
set(gca,'fontsize', 18);
%% USA
figure(2)
title('Time Averaged Map of NO_{2} Tropospheric Column (1/cm^2)')
latlim = [25 50];
lonlim = [-130 -65];
worldmap(latlim,lonlim) %create worldmap
%geolimits(latlim,lonlim)
load coastlines
pcolorm(lat,long,NO2) %plot NO2 data
hold on
geoshow(coastlat,coastlon,'Color','k') %show coastlines
colorbar
colormap jet
caxis([0 10.*10^15])
states = shaperead('usastatehi',...
       'UseGeoCoords', true, 'BoundingBox', [lonlim', latlim']);
geoshow(states,'FaceColor',[1,1,1],'facealpha',0);
set(gca,'fontsize', 18);
%% Coronavirus Data
ind= find(Cdata.Case~=0);
ind=ind(2:end);
cases=Cdata.Case(ind); % total cases
deaths = Cdata.Death(ind); % total deaths
deathP = deaths./cases; 
deathP=deathP.*100; % mortality per # cases
deathP(deathP>100)=0;
latD=Cdata.Lat(ind);
lonD=Cdata.Long(ind);

%% Correlation Calculations
% round latD and lonD to lar and long values
latDr = interp1(lat,lat,latD,'nearest');
lonDr = interp1(long,long,lonD,'nearest');

NO2_vect = zeros(length(deathP),1);
for i=1:length(latDr)
   k = find(lat==latDr(i));
   j = find(long==lonDr(i));
   %NO2_vect(i)=NO2(k,j);
   NO2_vect(i)=(NO2(k-1,j)+NO2(k+1,j)+NO2(k,j)+NO2(k-1,j-1)+NO2(k,j-1)+NO2(k+1,j-1)+...
       NO2(k-1,j+1)+NO2(k,j+1)+NO2(k+1,j+1))./9;
end

DeathCors = [latDr lonDr cases deaths deathP NO2_vect];
%% Bining
b1 = find(NO2_vect<10^15);
b2 = find(NO2_vect>10^15 & NO2_vect<2.*10^15);
b3 = find(NO2_vect>2.*10^15 & NO2_vect<3.*10^15);
b4 = find(NO2_vect>3.*10^15);

deathPB =[zeros(1,4)];
deathPB(1)=(sum(deaths(b1))./sum(cases(b1))).*100;
deathPB(2)=(sum(deaths(b2))./sum(cases(b2))).*100;
deathPB(3)=(sum(deaths(b3))./sum(cases(b3))).*100;
deathPB(4)=(sum(deaths(b4))./sum(cases(b4))).*100;

figure(3)
histogram('Categories',{'<1','1-2','2-3','>3'},'BinCounts',deathPB) 
xlabel('*10^(^1^5^)^ NO2 molecules per cm^2')
ylabel('Death percentage (%)')
grid on
title('US Coronavirus death percentage vs. tropospheric NO_{2} level')
set(gca,'fontsize', 18);
%%
figure(4)
deathsTot=[sum(deaths(b1)) sum(deaths(b2)) sum(deaths(b3)) sum(deaths(b4))];
histogram('Categories',{'<1','1-2','2-3','>3'},'BinCounts',deathsTot) 
xlabel('*10^(^1^5^)^ NO2 molecules per cm^2')
ylabel('Total Deaths')
grid on
title('US Coronavirus total deaths vs. tropospheric NO_{2} level')
set(gca,'fontsize', 18);



