

%% Load OMI and Italy's Coronavirus Data 
clear all
bbita=open('bbita.mat');
bbita=bbita.bbita;
NO2 = open('OMI2019.mat');
NO2=NO2.OMI2019;

lat = [-89.875:0.250:89.875];
long= [-179.875:0.250:179.875];


%% Italy Map
figure(1)
lonlim = [5.875 18.875];
latlim = [34.875 47.875];
worldmap(latlim,lonlim)
load coastlines
hold on
pcolorm(lat,long,NO2) %plot NO2 data
hold on
colorbar
colormap jet
caxis([0 10.*10^15])
geoshow(coastlat,coastlon,'Color','k') %show coastlines
title('Time Averaged Map of NO_{2} Tropospheric Column (1/cm^2)')


%% Italy Region Bounding boxes
%bounding boxes
lonB = bbita(:,1);
latB = bbita(:,2);
% round latD and lonD to lar and long values
latDr = interp1(lat,lat,latB,'nearest');
lonDr = interp1(long,long,lonB,'nearest');
i=0;
boundingbox=zeros(19,4);
for i=1:19
   boundingbox(i,1)=latDr(2*i-1);
   boundingbox(i,2)=lonDr(2*i-1);
   boundingbox(i,3)=latDr(2*i);
   boundingbox(i,4)=lonDr(2*i);
end

%% Compute Average NO2 Level Over Italy's Regions
totalRegion=zeros(length(boundingbox),1);
for i=1:length(boundingbox)
    for j=boundingbox(i,1):0.250:boundingbox(i,3)
        for k=boundingbox(i,2):0.250:boundingbox(i,4)
               m = find(lat==j);
               n = find(long==k);
               totalRegion(i)=totalRegion(i)+NO2(m,n);
        end
    end
    totalRegion(i)=totalRegion(i)./(length(boundingbox(i,2):0.250:boundingbox(i,4)).*length(boundingbox(i,1):0.250:boundingbox(i,3)));
end

%% COVID DATA
dcita=open('dcita.mat');
dcita=dcita.dcita;
deathPita=(dcita(:,1)./dcita(:,2)).*100;
deaths=dcita(:,1);
cases=dcita(:,2);
figure(2)
plot(totalRegion,deathPita,'*')
xlabel('NO2 leven per cm^2')
ylabel('Death percentage in Italy''s Region')
grid on
lsline
%% Bining
b1 = find(totalRegion<10^15);
b2 = find(totalRegion>10^15 & totalRegion<2.*10^15);
b3 = find(totalRegion>2.*10^15 & totalRegion<3.*10^15);
b4 = find(totalRegion>3.*10^15);

deathPB =[zeros(1,4)];
deathPB(1)=(sum(deaths(b1))./sum(cases(b1))).*100;
deathPB(2)=(sum(deaths(b2))./sum(cases(b2))).*100;
deathPB(3)=(sum(deaths(b3))./sum(cases(b3))).*100;
deathPB(4)=(sum(deaths(b4))./sum(cases(b4))).*100;

figure(3)
histogram('Categories',{'<1','1-2','2-3','>3'},'BinCounts',deathPB) 
xlabel('*10^(^1^5^)^ NO_{2} molecules per cm^2')
ylabel('Death percentage (%)')
grid on
title('Italy Coronavirus death percentage vs. tropospheric NO_{2} level (Regional Data)')
set(gca,'fontsize', 18);
