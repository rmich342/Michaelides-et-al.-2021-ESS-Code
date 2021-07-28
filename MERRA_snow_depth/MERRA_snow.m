

close all
clear all

%Read in MERRA-2 data
for i=1:12
    ind=i;
    if i<10
        name=strcat('MERRA2_400.tavgU_2d_lnd_Nx.20190',string(ind),'.nc4.nc4');
    else
        name=strcat('MERRA2_400.tavgU_2d_lnd_Nx.2019',string(ind),'.nc4.nc4');
    end
    lon=ncread(name,'lon');
    lat=ncread(name,'lat');
    time=ncread(name,'time');
    tsurf=ncread(name,'TSURF');
    tsnow=ncread(name,'TPSNOW');
    snowdepth=ncread(name,'SNODP');
    snowfrac=ncread(name,'FRSNO');
    
    lons{i,:}=lon;
    lats{i,:}=lat;
    times{i,:}=time;
    tsurfs{i,:}=tsurf;
    tsnows{i,:}=tsnow;
    snowdepths{i,:}=snowdepth;
    snowfracs{i,:}=snowfrac;
    
    avg_depth(:,:,i)=mean(snowdepths{i},3);
    avg_frac(:,:,i)=mean(snowfracs{i},3);
end


%bounding box
%((-151.4259 68.5275,-148.56 68.5275,-148.56 69.5649,-151.4259 69.5649,-151.4259 68.5275))
lat_min=68.5275;
lat_max=69.5649;
lon_min=-148.56;
lon_max=-151.4259;

[a lon_min_ind]=min(abs(lon-lon_min));
[a lon_max_ind]=min(abs(lon-lon_max));
[a lat_min_ind]=min(abs(lat-lat_min));
[a lat_max_ind]=min(abs(lat-lat_max));

months={'January','February','March','April','May','June','July','August','September','October','November','December'};



%Plot Average Snow Depth and Snow Cover Fraction for each Month
figure('units','normalized','outerposition',[0 0 1 .8])
for i=1:12
    subplot(1,2,1),imagesc(lon,lat,avg_depth(:,:,i).')
    hold on
    rectangle('Position',[lon_min lat_min abs(lon_max-lon_min) abs(lat_max-lat_min)],'LineWidth',3,...
    'EdgeColor',[1 0 0]);
    colorbar
    caxis([0 1])
    title(strcat(string(months(i)),' Snow Depth (m)'))
    xlabel('Longitude')
    ylabel('Latitude')
    set(gca,'YDir','normal')
    set(gca,'FontSize', 20)
    subplot(1,2,2),imagesc(lon,lat,avg_frac(:,:,i).')
    hold on
    rectangle('Position',[lon_min lat_min abs(lon_max-lon_min) abs(lat_max-lat_min)],'LineWidth',3,...
    'EdgeColor',[1 0 0]);
    colorbar
    %caxis([0 1])
    title(strcat(string(months(i)),' Snow Fraction (-)'))
    xlabel('Longitude')
    ylabel('Latitude')
    set(gca,'YDir','normal')
    set(gca, 'FontSize', 20)
    pause();
end




subplot(1,2,1),imagesc(lon,lat,(avg_depth(:,:,10)-avg_depth(:,:,4)).')
hold on
rectangle('Position',[lon_min lat_min abs(lon_max-lon_min) abs(lat_max-lat_min)],'LineWidth',3,...
'EdgeColor',[1 0 0]);
colorbar
caxis([-1 1])
title(strcat(string(months(i)),' Snow Depth (m)'))
xlabel('Longitude')
ylabel('Latitude')
set(gca,'YDir','normal')
set(gca,'FontSize', 20)
subplot(1,2,2),imagesc(lon,lat,(avg_frac(:,:,10)-avg_frac(:,:,4)).')
hold on
rectangle('Position',[lon_min lat_min abs(lon_max-lon_min) abs(lat_max-lat_min)],'LineWidth',3,...
'EdgeColor',[1 0 0]);
colorbar
%caxis([0 1])
title(strcat(string(months(i)),' Snow Fraction (-)'))
xlabel('Longitude')
ylabel('Latitude')
set(gca,'YDir','normal')
set(gca, 'FontSize', 20)
pause();



