%% Load all of the InSAR files you are considering
clear all; close all; clc;

nr=1587; naz=1515;% image size:  nr = # of x (range) pixels; 
                             % naz = # of y (azimuth) pixels
n=74; % number of slcs 

%Just read in coherent inteferograms from the 2019 year with a temporal
%baseline of 48 days or less
cells=importdata('intlist48_2019');
N=length(cells);
lambda=5.6; %wavelength

unw_phase=zeros(nr,naz,N);
amp=zeros(nr,naz,N);
coh=zeros(nr,naz,N);
unis=zeros(nr,naz,N);
date_pair=cell(2,N);
doy_pair=cell(2,N);

%Read in the unwrapped phase (unw), coherence (coh), amplitude (amp) and
%unimodally-corrected unwrapped phase (uni)
for i=1:N
    strint=cells{i};
    strunw=strrep(strint,'.int','.unw');
    stramp=strrep(strint,'.int','.amp');
    strcc=strrep(strint,'.int','.cc'); 
    struni=strrep(strint,'.int','.unimodal.unw');
% correlations
    filename_c=sprintf('%s',strcc); 
    fid=fopen(filename_c);
    dat=fread(fid,[2*nr,inf],'float','ieee-le');
    temp=dat((nr+1):end,:);
    coh(:,:,i)=temp;
    fclose(fid);
% unwrapped phase
    filename=sprintf('%s',strunw);  
    fid=fopen(filename);
    dat=fread(fid,[2*nr,inf],'float','ieee-le');
    temp=dat(nr+1:end,:);
    unw_phase(:,:,i)=temp;
    fclose(fid);
% Amplitude
    filename=sprintf('%s',stramp);  
    fid=fopen(filename);
    dat=fread(fid,[2*nr,inf],'float','ieee-le');
    temp=dat(1:2:2*nr-1,:)+dat(2:2:2*nr,:);
    amp(:,:,i)=temp;
    fclose(fid);
% unimodal unwrapped phase
    filename=sprintf('%s',struni);  
    fid=fopen(filename);
    dat=fread(fid,[2*nr,inf],'float','ieee-le');
    temp=dat(1:nr,:);
    uni(:,:,i)=temp;
    fclose(fid);
    
% date information
split1=strsplit(strint,'_');
strint2=split1{2};
split2=strsplit(strint2,'.');
d1=split1{1};
d2=split2{1};

date1=strcat(d1(5:6),'/',d1(7:8),'/',d1(1:4));
date2=strcat(d2(5:6),'/',d2(7:8),'/',d2(1:4));

date1_vec=datetime(date1,'InputFormat','MM/dd/yyyy');
date2_vec=datetime(date2,'InputFormat','MM/dd/yyyy');

doy1=day(date1_vec,'dayofyear');
doy2=day(date2_vec,'dayofyear');

date_pair{1,i}=date1;
date_pair{2,i}=date2;
doy_pair{1,i}=doy1;
doy_pair{2,i}=doy2;
end


% Get an average coherence file 
avecc = mean(coh,3);
cor_mask=avecc;
cor_mask(cor_mask<0.15)=nan;
cor_mask(cor_mask>0.15)=1;

% Get an average amplitude file 
amp_mean = mean(amp,3);
amp_mask=amp_mean;
amp_mask(amp_mask<0.1)=nan;
amp_mask(amp_mask>0.1)=1;

%generated combined amplitude and coherence mask
mask=cor_mask+amp_mask;
mask(mask<2)=nan;
mask(mask==2)=1;

%% Read in Temperature data, Daymet ADDT Info

data=csvread('14235_lat_69.0_lon_-150.0_2020-08-31_204431.csv',8);
year=data(:,1);
yday=data(:,2);
dayls=data(:,3);
prcpmm=data(:,4);
sradWm2=data(:,5);
swekgm2=data(:,6);
tmaxdegc=data(:,7);
tmindegc=data(:,8);
vpPa=data(:,9);

tavg=(tmindegc+tmaxdegc)./2;    %average daily temperature
thaw_count=zeros(365,1);        %number of thaw days
freeze_count=zeros(365,1);      %number of freeze days
ADDT=zeros(365,1);
ADDF1=zeros(365,1);
ADDF2=zeros(365,1);

thaw_count(1,1)=0;
freeze_count(1,1)=1;
ADDT(1,1)=0;
ADDF1(1,1)=0-tavg(1,1);
ADDF2(1,1)=0-tavg(1,1);

%Go through 2019 temperature record and calculate Accumulated degree days
%of thaw (ADDT) and accumulated degree days of freezing (ADDF)
for i=2:365
    if tavg(i,1)>0
        thaw_count(i,1)=thaw_count(i-1,1)+1;
        ADDT(i,1)=ADDT(i-1,1)+tavg(i,1);
        ADDF1(i,1)=ADDF1(i-1,1);
        ADDF2(i,1)=ADDF1(i-1,1);
    end
    if tavg(i,1)<0
        freeze_count(i,1)=freeze_count(i-1,1)+1;
        ADDT(i,1)=ADDT(i-1,1);
        ADDF1(i,1)=ADDF1(i-1,1)+(0-tavg(i,1));
        ADDF2(i,1)=ADDF2(i-1,1)+(0-tavg(i,1));
    end
    if tavg(i,1)==0
        ADDT(i,1)=ADDT(i-1,1);
        ADDF1(i,1)=ADDF1(i-1,1);
        ADDF2(i,1)=ADDF1(i-1,1);
        freeze_count(i,1)=freeze_count(i-1,1);
        thaw_count(i,1)=thaw_count(i-1,1);
    end
        
end

%Zero out ADDF2 to just consider autumn refreeze
ADDF2=ADDF1-1714;
ADDF2(ADDF2<0)=0;


%generate temperature/subsidence proxy curves
temp_curve=sqrt(ADDT)-4*sqrt(ADDF2);
temp_curve(temp_curve<0)=0;

N_temp_curve=temp_curve./max(temp_curve);       %Temperature curve with thaw and freeze
NADDT=sqrt(ADDT./max(ADDT));                    %temperature curve just with thaw



%visualize the temperature/subsidence curves
figure
subplot(2,1,1)
plot(NADDT,'LineWidth',2,'Color','r')
hold on
plot(N_temp_curve,'LineWidth',2,'Color','b')
xlabel('Day of Year')
xlim([1 365]);
ylabel({'Normalized Accumulated','Degree Days of Thaw'})
title('Temperature Curves')
legend('Pure Thaw','Thaw and Freeze')
set(gca,'FontSize',20);
hold on
subplot(2,1,2)
plot(1-NADDT,'LineWidth',2,'Color','r')
hold on
plot(1-N_temp_curve,'LineWidth',2,'Color','b')
xlabel('Day of Year')
xlim([1 365]);
ylabel({'Normalized Proxy','Seasonal Subsidence'})
title('Subsidence Curves')
legend('Pure Thaw','Thaw and Freeze')
set(gca,'FontSize',20);


%visualize the date pairs of all interferograms
figure
plot(cell2mat(doy_pair(1,:)),'LineWidth',2)
hold on
plot(cell2mat(doy_pair(2,:)),'LineWidth',2)
hold on
plot(cell2mat(doy_pair(2,:))-cell2mat(doy_pair(1,:)),'LineWidth',2)
legend('DOY 1','DOY 2','Temporal Baseline')
title('Day of Year Pairs, InSAR')
xlabel('Interferogram Number')
ylabel('Day of Year')
set(gca,'FontSize',20)


%% Now we will do the atmospheric correction
%remove topographically-correlated atmospheric noise from interferograms
% Depending on how large the scene is, you may want to base the correction
% solely on a few areas with the topographic relief changes quite a bit.
% But if the scene is small, you can indeed use the whole scene.
% This bit also requires a calibration pixel (or set of pixels) to ensure
% that all of the scenes are set to the same "datum"

% Load the dem
nr0=28568;
naz0=9094;
fid=fopen('elevation.dem','r');
dem0=fread(fid,[nr0,naz0],'int16'); % x length first, y length second
fclose(fid);

dem=imresize(dem0,[nr,naz]);
masked_unw_phase = unw_phase.*mask;
masked_unw_phase1 = uni.*mask;

% Pick some pixels for calibration
pixels = [1140 1070]; % [range azimuth]
[sz2,~] = size(pixels);

phase = zeros(nr,naz,N);
phase1 = zeros(nr,naz,N);
corrections=zeros(nr,naz,N);
corrections1=zeros(nr,naz,N);
phase_rshp = [];
dem_rshp = [];
for int = 1:N
    block_phase = masked_unw_phase(:,:,int);
    indx = isnan(block_phase);
    line = polyfit(dem(~indx),block_phase(~indx),1);
    correction = (line(1)*dem + line(2)); 
    corrections(:,:,int)=correction;
    %phase(:,:,int) = masked_unw_phase(:,:,int);
    phase(:,:,int) = masked_unw_phase(:,:,int) - correction;
    
    block_phase = masked_unw_phase1(:,:,int);
    indx = isnan(block_phase);
    line = polyfit(dem(~indx),block_phase(~indx),1);
    correction = (line(1)*dem + line(2)); 
    corrections1(:,:,int)=correction;
    %phase1(:,:,int) = masked_unw_phase1(:,:,int);
    phase1(:,:,int) = masked_unw_phase1(:,:,int) - correction;
    
end


%% Histogram Phase Calibration
%Calibrate the scene-wide phase based on the 5th percentile of measured
%deformation
% Pick some pixels for calibration
phase_cal = zeros(nr,naz,N);
phase_cal1 = zeros(nr,naz,N);
for int = 1:N
    pluck=phase(:,:,int);
    val = prctile(pluck(:),5);
    phase_cal(:,:,int) = pluck - val;
    
    
    pluck=phase1(:,:,int);
    val = prctile(pluck(:),5);
    phase_cal1(:,:,int) = pluck - val;
end

%% Pluck out Scenes you want to use
%Only use the interferograms corresponding to the following indices
%below is for intlist48_2019
scene_inds=[21,24,25,27,30,41,42,45,49,53,54,57,64,77];


phase_good=phase_cal(:,:,scene_inds);
phase1_good=phase_cal1(:,:,scene_inds);
doy_pair_good=doy_pair(:,scene_inds);


%% ADDT-modified SBAS (Liu et al. 2012)
%calculate the best-fitting seasonal subsidence estimate for each pixel
%using a modified form of the SBAS algorithm.


%Generate design matrix using NADDT
B=zeros(N,1);
for i=1:N
    %B(i,1)=doy_pair{2,i}-doy_pair{1,i};
    B(i,1)=NADDT(doy_pair{2,i})-NADDT(doy_pair{1,i});
end

B1=B(scene_inds,1);
B1=B1./norm(B1);
Bi=pinv(B);
B1i=pinv(B1);


phase_good=phase(:,:,scene_inds);
phase1_good=phase1(:,:,scene_inds);

%ref_phase=phase_good(1051,493);
%ref_phase1=phase1_good(1051,493);


seasonal=zeros(nr,naz);
seasonal1=zeros(nr,naz);

for i=1:nr
    for j=1:naz
        %pluck=squeeze(phase_good(i,j,:)-ref_phase);
        pluck=squeeze(phase_good(i,j,:));
        sol=squeeze(B1i*pluck);
        seasonal(i,j)=sol;
        
        %pluck=squeeze(phase1_good(i,j,:)-ref_phase1);
        pluck=squeeze(phase1_good(i,j,:));
        sol=squeeze(B1i*pluck);
        seasonal1(i,j)=sol;
    end
end
        

seasonal_cm = (seasonal.*(lambda./(4*pi)));
seasonal1_cm = (seasonal1.*(lambda./(4*pi)));

seasonal_cm_ref=seasonal_cm-prctile(seasonal_cm(:),5);
seasonal1_cm_ref=seasonal1_cm-prctile(seasonal1_cm(:),5);


figure
subplot(2,2,1),imagesc(seasonal_cm_ref.')
title({'2019 ADDT-modified SBAS','Deformation','(cm)'})
caxis([0 5])
colorbar
axis off
set(gca,'FontSize',20)
subplot(2,2,2),histogram(seasonal_cm_ref.')
title({'2019 ADDT-modified','Deformation','(cm)'})
set(gca,'FontSize',20)
subplot(2,2,3),imagesc(seasonal1_cm_ref.')
title({'2019 ADDT-modified','Deformation','(cm)','unimodal correction'})
caxis([0 5])
colorbar
axis off
set(gca,'FontSize',20)
subplot(2,2,4),histogram(seasonal1_cm_ref.')
title({'2019 ADDT-modified','Deformation','(cm)','unimodal correction'})
set(gca,'FontSize',20)



figure
subplot(1,2,1),imagesc(seasonal1_cm_ref.')
title({'2019 Seasonal','Subsidence (cm)'})
caxis([0 5])
colorbar
axis off
set(gca,'FontSize',20)
subplot(1,2,2),histogram(seasonal1_cm_ref.')
title({'2019 Seasonal','Subsidence (cm)'})
set(gca,'FontSize',20)



%% Calculate residuals from ADDT-modified SBAS (Liu et al. 2012)
%create forward model interferograms

data_hat=zeros(nr,naz,length(B1));

for i=1:length(B1)
    data_hat(:,:,i)=B1(i).*seasonal1_cm_ref;
end

data=phase1_good.*(lambda./(4*pi));

rmse=sqrt(nanmean((data-data_hat).^2,3));
%(y - yhat)    % Errors
%(y - yhat).^2   % Squared Error
%mean((y - yhat).^2)   % Mean Squared Error
%RMSE = sqrt(mean((y - yhat).^2));  % Root Mean Squared Error

residuals=data-data_hat;
sub=seasonal1_cm_ref;
seas_mask=seasonal1_cm_ref;

seas_mask(rmse>abs(seas_mask))=nan;


figure
subplot(2,2,1),imagesc(sub.')
title({'2019 ADDT-Modified SBAS','Seasonal Subsidence (cm)'})
caxis([0 5])
colorbar
axis off
set(gca,'FontSize',20)
subplot(2,2,2),histogram(sub.')
title({'2019 ADDT-Modified SBAS','Seasonal Subsidence (cm)'})
xlim([-5 10])
set(gca,'FontSize',20)
subplot(2,2,3),imagesc(rmse.')
title({'2019 ADDT-Modified SBAS','RMSE (cm)'})
caxis([0 5])
colorbar
axis off
set(gca,'FontSize',20)
subplot(2,2,4),histogram(rmse)
title({'2019 ADDT-Modified SBAS','RMSE (cm)'})
xlim([0 5])
set(gca,'FontSize',20)


%% Read in SBAS parameters and make them consistent with the other data
load Bperp.out
load Tm.out
load deltime.out
load timedeltas.out

Bperp0=Bperp;
Tm0=Tm;
deltime0=deltime;
timedeltas0=timedeltas;

Bperp=zeros(size(Bperp));
Tm=zeros(size(Tm));
deltime=zeros(size(deltime));


cells2=importdata('sbas_list');
N2=length(cells2);

sbas_dates0=string(cells2.textdata);
sbas_data0=string(cells2.data);

%Generate matrixes of temporal baselines, perpendicular baselimes, and
%delta times
for i=1:N
    % date information
    strint=cells{i};   
    split1=strsplit(strint,'_');
    strint2=split1{2};
    split2=strsplit(strint2,'.');
    d1=split1{1};
    d2=split2{1};
    
    index1=find(contains(sbas_dates0(:,1),d1));
    index2=find(contains(sbas_dates0(:,2),d2));
    index=intersect(index1,index2);
    
    Bperp(i,:)=Bperp0(index,:);
    Tm(i,:)=Tm0(index,:);
    deltime(i,:)=deltime0(index,:);
end


cellsn=importdata('geolist');
n=length(cellsn);

%get dayofyear information for each SLC
for i=1:n
    strint=cellsn{i};
    d1=strint(21:28);
    date1=strcat(d1(5:6),'/',d1(7:8),'/',d1(1:4));
    date1_v=datetime(date1,'InputFormat','MM/dd/yyyy');
    doy1_v=day(date1_v,'dayofyear');
    
    geo_date{i}=date1_v;
    geo_doy{i}=doy1_v;
end
