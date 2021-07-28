clear all
close all
%% post processing
along_track = csvread('repeat_dh_1280_gt2l_2m_dem.csv', 1, 4); %along-track repeats
lon=along_track(1:end,1);
lat=along_track(1:end,2);
rgt = along_track(1:end, 22);
beam = along_track(1:end, 23); %beam 

plot_title = 'Track 1280 gt2l:  6/21 to 9/19'; %title for the plots
track_name = '1280_gt2l'; %for nameing files

along_track_s = sortrows(along_track, 2); %sort by latitude

%extract variables from subsetted data
lon = along_track_s(1:end,1);
lat = along_track_s(1:end,2);
dh=along_track_s(1:end,9) *100; %convert to cm
sigma_dh = along_track_s(1:end,10) *100; % dh uncertatinty, converted to cm
h1 = along_track_s(1:end,5); %Height measured by IS2 on the first pass
h2 = along_track_s(1:end,7); %Height measured by IS2 on the second pass
elevation_profile = along_track_s(1:end, 24); %height sampled from Arctic DEM on first pass
t1=along_track_s(1:end,11); %time of first IS2 pass
t2=along_track_s(1:end,12); %time of 2nd IS2 pass
dtd=along_track_s(1:end,13); %time difference between two passes
%get average along-track slope estimated from IS2  (cm/cm)
slope_is2_1 = along_track_s(1:end, 14);
slope_is2_2 = along_track_s(1:end, 15);
slope_is2 = .5 * (slope_is2_1 + slope_is2_2);
%get average across-track slope estimated from IS2 (cm/cm)
slope_across_1 = along_track_s(1:end, 16);
slope_across_2 = along_track_s(1:end, 17);
slope_across = .5*(slope_across_1+slope_across_2);
%convert times to days since 1/1/2019
doy1=floor((t1-2019).*365);
doy2=floor((t2-2019).*365);


%remove large outliers
index  =((abs(dh) < 1000) & (abs(h1) <10000));
lon = lon(index);
lat = lat(index);
dh = dh(index);
h1 = h1(index);
h2 = h2(index);
sigma_dh = sigma_dh(index);
elevation_profile = elevation_profile(index);
slope_is2 = slope_is2(index);
slope_across = slope_across(index);
t1 = t1(index);
t2 = t2(index);
dtd = dtd(index);
doy1 = doy1(index);
doy2 = doy2(index);

%mask anomolous across-track slopes
slope_across(slope_across > 100) = nan;

%set up arrays for moving average
along_track = abs(lldistkm([lat(1) lon(1)], [lat lon]))*1000; %along-track distance in meters
max_distance = abs(lldistkm([lat(end) lon(end)], [lat(1) lon(1)]))*1000; %get max along-track distance
%make new arrays for data spaced every 20m along-track
dh_spaced = [];
slope_is2_spaced = [];
sigma_dh_spaced = [];
lat_spaced = [];
lon_spaced = [];
slope_across_spaced = [];
%start counters
d_min = 0;
i = 1;
while d_min < max(along_track) %loop through along-track distance
    %define range of 20m bin
    d_min = (i-1) *20;
    d_max = i*20;
    %find value(s) that fall in range(if they were perfectly spaced, there would only be one, but sometimes there are 2 due to < 20m spacing)
    dh_point = dh((along_track >= d_min) & (along_track <= d_max)); 
    indices = find((along_track >= d_min) & (along_track <= d_max)); %get indices for points in 20-m bin
    %nan-fill for intervals with no data
    if length(dh_point) ==0
       dh_spaced = [dh_spaced NaN]; 
       sigma_dh_spaced = [sigma_dh_spaced NaN];
       slope_is2_spaced = [slope_is2_spaced NaN];
       lat_spaced = [lat_spaced NaN];
       lon_spaced = [lon_spaced NaN];
       slope_across_spaced = [slope_across_spaced NaN];
    end
    if length(dh_point) > 0 
        index = indices(1); %get index of first point in bin in case there are two
        %populated spaced arrays
        dh_spaced = [dh_spaced dh_point(1)];
        sigma_dh_spaced = [sigma_dh_spaced sigma_dh(index)];
        slope_is2_spaced = [slope_is2_spaced slope_is2(index)];
        lat_spaced = [lat_spaced lat(index)];
        lon_spaced = [lon_spaced lon(index)];
        slope_across_spaced = [slope_across_spaced slope_across(index)];
    end
    if length(dh_point)==2 %get index of second point in bin if there are two
        index = indices(2);
        dh_spaced = [dh_spaced dh_point(2)];
        sigma_dh_spaced = [sigma_dh_spaced sigma_dh(index)];
        slope_is2_spaced = [slope_is2_spaced slope_is2(index)];
        lat_spaced = [lat_spaced lat(index)];
        lon_spaced = [lon_spaced lon(index)];
        slope_across_spaced = [slope_across_spaced slope_across(index)];
    end
    if length(dh_point)>2 %I don't think this ever happens but just in case
       disp('more than 2 elements found')
       disp(dh_point)
    end
    i = i + 1;
end

%run 100-point (2 km) boxcar
dh_av_2000 = movmean(dh_spaced, 100, 'omitnan');
slope_is2_av = movmean(slope_is2_spaced, 100, 'omitnan');
slope_across_av = movmean(slope_across_spaced, 100, 'omitnan');

%propogate the dh uncertainty through the moving mean:
sigma_av = nan(1, length(sigma_dh_spaced)); %initialize array

for i=1:length(sigma_dh_spaced)
    if i <= 50 %case where the center of the boxcar is < 50 points away from the start of the dataset
        sigmas = sigma_dh_spaced(1:i+49);
    elseif(i <= length(sigma_dh_spaced) - 49) && (i > 50) %case where we have full 100-point interval in the boxcar
        sigmas  = sigma_dh_spaced(i-50:i+49);
    elseif(i > length(sigma_dh_spaced) - 49) %case where the center of the boxcar is > 50 points away from the start of the dataset
        sigmas = sigma_dh_spaced(i:end);
    end
    %propogate uncertainty in the interval
    sigmas = sigmas(~isnan(sigmas));
    sigma_av(i) = sqrt(sum(sigmas.^2))/(length(sigmas));
end

%interoplate lat/lon of points where there was no data before averaging
for i=1:length(lat_spaced)
    indices = linspace(1, length(lat_spaced), length(lat_spaced));
    %linear fit of coordinates
    p_lat = polyfit(indices(~isnan(lat_spaced)), lat_spaced(~isnan(lat_spaced)), 1);
    p_lon = polyfit(indices(~isnan(lon_spaced)), lon_spaced(~isnan(lon_spaced)), 1);
    if isnan(lat_spaced(i)) %evaluate at points with no lon/lat value
        lat_spaced(i) = polyval(p_lat, i);
        lon_spaced(i) = polyval(p_lon, i);
    end
end

%sample SBAS time series at averaged is2 points

%get dimensions for full InSAR domain
X_FIRST=-152.09892273;
Y_FIRST=69.72755432;
X_STEP=0.00268968078;
Y_STEP=-0.00089656026;

nr=1587; naz=1515;
%set up arrays for InSAR cooridinates
lon_vec=zeros(nr,1);
lon_vec(1,1)=X_FIRST;
for i=2:nr
    lon_vec(i)=lon_vec(i-1)+X_STEP;
end
lat_vec=zeros(naz,1);
lat_vec(1,1)=Y_FIRST;
for i=2:naz
    lat_vec(i)=lat_vec(i-1)+Y_STEP;
end

%load in InSAR data
load('SBAS_times.mat');%geo_doy_2019 (times used in the SBAS time series)
load('SBAS_time_series_roger.mat'); %time_series (SBAS time series in cm)


%put IS2 values into array with same xy dimensions as the InSAR data. The
%third dimension allows for storage of multiple IS2 points in a sinlge
%lat/lon grid cell (assuming there's a max of 7 per cell)

%set up component vectors to grid IS2 values, 

count_vec=zeros(nr,naz,2); %tracks number of data points in the current lat/lon bin
dh_vec=zeros(nr,naz,7); %sampled SBAS dh values
l_vec =zeros(nr, naz, 7); %sampled latitude






for i=1:length(dh_av_2000)
    %get IS2 lon, lat, and dh
    lon_pluck=lon_spaced(i);
    lat_pluck=lat_spaced(i);
    dh_pluck=dh_av_2000(i);
    t1_pluck = doy1(1);
    t2_pluck = doy2(2);

   
    %find and get index for closest grid cell to IS2 point 
    [a lon_ind]=min(abs(lon_vec-lon_pluck));
    [b lat_ind]=min(abs(lat_vec-lat_pluck));
    
    
    count_vec(lon_ind,lat_ind,1)=count_vec(lon_ind,lat_ind,1)+1; %increment count of number of IS2 points in the current grid cell
    val=count_vec(lon_ind,lat_ind,1); %retieve current count
    dh_vec(lon_ind,lat_ind,val)=dh_pluck; %store IS2 dh value for this cell. If there are multiple data points in a single lat/lon bin, each point is stored in a new layer(val)
    l_vec(lon_ind,lat_ind,val) = lat_pluck; %store latitude of IS2 point
    
end
%convert empty entries to nan
dh_vec(dh_vec==0)=nan;
l_vec(l_vec == 0) = nan;

%initialize variables
SBAS_dh = []; %for storing dh from SBAS time series
lat_f = []; %latitude of sampled points

%get dates of the two IS2 passes
date_1 = doy1(1);
date_2 = doy2(1);
%find closest date in sbas time series
[m, date_1_i] = min(abs(geo_doy_2019 - date_1));
[m, date_2_i] = min(abs(geo_doy_2019 - date_2));
%loop through grid cells
for i =1:size(dh_vec, 1)
    for j = 1:size(dh_vec, 2)
        for k = 1:size(dh_vec,3)
            if(~isnan(dh_vec(i, j,k))) %only sample pixels with valid IS2 measurement 
                lat_f = [lat_f l_vec(i, j, k)]; %latitude of IS2 point               
                SBAS_dh = [SBAS_dh (time_series(i, j, date_2_i) - time_series(i, j, date_1_i))]; %elevation change at select grid cell between date 1 and 2 according to SBAS time series                
            end
        end
    end
end

%re-sort by latitude
[out, idx] = sort(lat_f);
 
%Apply snow correction

%Read in MERRA-2 data
for i=1:12 %loop through months
    ind=i;
    if i<10
        name=strcat('MERRA2_400.tavgU_2d_lnd_Nx.20190',string(ind),'.nc4.nc4');
    else
        name=strcat('MERRA2_400.tavgU_2d_lnd_Nx.2019',string(ind),'.nc4.nc4');
    end
    lonp=ncread(name,'lon'); 
    latp=ncread(name,'lat');
    time=ncread(name,'time');
    tsurf=ncread(name,'TSURF'); 
    tsnow=ncread(name,'TPSNOW'); 
    snowdepth=ncread(name,'SNODP'); %snow depth
    snowfrac=ncread(name,'FRSNO'); %snow cover fraction
    
    lons{i,:}=lonp;
    lats{i,:}=latp;
    times{i,:}=time;
    tsurfs{i,:}=tsurf;
    tsnows{i,:}=tsnow;
    snowdepths{i,:}=snowdepth;
    snowfracs{i,:}=snowfrac;
    %get mean snow depth and snow cover fraction
    avg_depth(:,:,i)=mean(snowdepths{i},3);
    avg_frac(:,:,i)=mean(snowfracs{i},3);
end


%bounding box
%((-151.4259 68.5275,-148.56 68.5275,-148.56 69.5649,-151.4259 69.5649,-151.4259 68.5275))
lat_min=68.5275;
lat_max=69.5649;
lon_min=-148.56;
lon_max=-151.4259;

% [a lon_min_ind]=min(abs(lon-lon_min));
% [a lon_max_ind]=min(abs(lon-lon_max));
% [a lat_min_ind]=min(abs(lat-lat_min));
% [a lat_max_ind]=min(abs(lat-lat_max));

months={'January','February','March','April','May','June','July','August','September','October','November','December'};

%coordinates of Merra-2 grid
lons1=lons{1};
lats1=lats{1};

delta_snow_transect=zeros(1,length(dh_av_2000));


%Sample Merra-2 at IS2 points
for i=1:length(lat_spaced)
    %get coordinates of IS2 point
    lon_pluck=lon_spaced(i);
    lat_pluck=lat_spaced(i);
    
    %get months of IS2 passes
    month1 = ceil(date_1./31);
    month2 = ceil(date_2./31);
   
    %get indices for grid 
    [a lon_ind]=min(abs(lons1-lon_pluck));
    [b lat_ind]=min(abs(lats1-lat_pluck));
    
    %get snow depth for the given time/position
    snow1=avg_depth(lon_ind,lat_ind,month1);
    snow2=avg_depth(lon_ind,lat_ind,month2);
    %get chamge in snow depth
    delta_snow_transect(1,i)=snow2-snow1;
    
    
end
%subtract snow depths from IS2
is2_corrected = dh_av_2000 - 100.*delta_snow_transect;


%% Plot output
figure('Renderer', 'painters', 'Position', [5 15 1200 700])
subplot(4, 1, 1) %plot along-track elevation profile
scatter(lat, elevation_profile,10, 'filled')
hold on
scatter(lat, h1, 2, 'k', 'filled')
hold on
xlim([min(lat) max(lat)])
ylim([0 600])
legend('Arctic DEM','ICESat-2', 'location', 'southwest', 'fontsize', 16, 'orientation', 'horizontal')
legend('boxoff')
ax = gca; 
ax.YAxis.FontSize = 16;
ax.XAxis.FontSize = 16;
ylabel('elevation (m)', 'fontsize', 20)

subplot(4, 1, 2) %plot along and across-track slopes
scatter(lat_spaced, slope_is2_av, 5, 'filled')
hold on
ax = gca; 
ax.YAxis.FontSize = 16;
ax.XAxis.FontSize = 16;
xlim([min(lat) max(lat)]) 
ylim([-.10 .10])
scatter(lat_spaced, slope_across_av,5, 'filled')
ylabel('slope (cm/cm)')
legend('along-track', 'across-track', 'fontsize', 16, 'orientation', 'horizontal', 'location', 'southeast')
legend('boxoff')


subplot(4, 1, 3:4) %plot SBAS and IS2 dh together
b = scatter(lat_f(idx), SBAS_dh(idx),5, 'red', 'filled');
hold on
ylabel('$dh (cm)$')
e = plot([min(lat) max(lat)], [0 0], 'k');
hold on
f = scatter(lat_spaced, is2_corrected, 5, 'b', 'filled');
legend([b f], {'SBAS', 'snow-corrected IS2'}, 'location', 'northwest', 'fontsize', 16, 'orientation', 'horizontal')
legend('boxoff')
hold on
ylim([-50 50])
xlim([min(lat) max(lat)])
ax = gca; 
ax.YAxis.FontSize = 16;
ax.XAxis.FontSize = 16;
sgtitle(plot_title, 'fontsize', 18,'fontsize', 20);
xlabel('Latitude', 'interpreter','latex','fontsize', 20)
ylabel('dh (cm)','fontsize', 20)
%% save processed is2/SBAS datasets as mat files

save(['track_' track_name '_smoothed_is2_minus_snow.mat'], 'is2_corrected');
save(['track_' track_name '_sampled_SBAS.mat'], 'SBAS_dh');
save(['track_' track_name '_latitude.mat'], 'lat_spaced');
save(['track_' track_name '_SBAS_latitude.mat'], 'lat_f');
save(['track_' track_name '_slope_along.mat'], 'is2_slope_av');
save(['track_' track_name '_slope_across.mat'], 'slope_across_av');





