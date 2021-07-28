close all
clear all
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



%xovers=csvread('xovers_two_seasons_all_slopes_f.csv',1,0);
%dataset='2019/2020 Thaw Season, All Slopes'

%xovers=csvread('xovers_two_seasons_filtered_slope_f.csv',1,0);
%dataset='2019/2020 Thaw Season, |Slopes|< 5 \circ'

%xovers=csvread('xovers_full_period_all_slopes_f.csv',1,0);
%dataset='Full Period, All Slopes'

%xovers=csvread('xovers_full_period_filtered_slope_f.csv',1,0);
%ataset='Full Period, |Slopes|< 5 \circ'


xovers=csvread('xovers_full_period_low_slopes_new_filtered.csv',1,0);
dataset='Full Period, |Slopes|< 5 \circ'

lon_up=xovers(1:end,1);
lat_up=xovers(1:end,2);
dh_up=xovers(1:end,3)*100;
t1_up=xovers(1:end,7);
t2_up=xovers(1:end,8);
dtd_up=xovers(1:end,10);
doy1_up=floor((t1_up-2019).*365);
doy2_up=floor((t2_up-2019).*365);









Btemp=365;      %Temporal baseline threshold
%claculate statistics
average = mean(dh_up);
two_sigma=2*std(dh_up);        


%filter crossovers by temporal baseline and filter out values outside the
%two sigma range
dh_2019_up=[];
doy1_2019_up=[];
doy2_2019_up=[];
lon_2019_up=[];
lat_2019_up=[];

dh_2020_up=[];
doy1_2020_up=[];
doy2_2020_up=[];
lon_2020_up=[];
lat_2020_up=[];

for i=1:length(dh_up)
    if (dh_up(i))<=(average + two_sigma) && (dh_up(i)) >= (average - two_sigma) && doy2_up(i)-doy1_up(i)<Btemp
        %2019 year
        if doy1_up(i)>0 && doy1_up(i)<365
            if doy2_up(i)>0 && doy2_up(i)<365
                dh_2019_up=[dh_2019_up,dh_up(i)];
                doy1_2019_up=[doy1_2019_up,doy1_up(i)];
                doy2_2019_up=[doy2_2019_up,doy2_up(i)];
                lon_2019_up=[lon_2019_up,lon_up(i)];
                lat_2019_up=[lat_2019_up,lat_up(i)];
            end
        end
        %2020 year
        if doy1_up(i)>365 && doy1_up(i)<730
            if doy2_up(i)>365 && doy2_up(i)<730
                dh_2020_up=[dh_2020_up,dh_up(i)];
                doy1_2020_up=[doy1_2020_up,doy1_up(i)];
                doy2_2020_up=[doy2_2020_up,doy2_up(i)];
                lon_2020_up=[lon_2020_up,lon_up(i)];
                lat_2020_up=[lat_2020_up,lat_up(i)];
            end
        end
    end
end



 
%generate their associated temperature curve values
is2_addt_up_2019=sqrt(NADDT(doy2_2019_up))-sqrt(NADDT(doy1_2019_up));
is2_Ntemp_up_2019=(N_temp_curve(doy2_2019_up))-(N_temp_curve(doy1_2019_up));
is2_addt_up_2020=sqrt(NADDT(doy2_2020_up-365))-sqrt(NADDT(doy1_2020_up-365));
is2_Ntemp_up_2020=(N_temp_curve(doy2_2020_up-365))-(N_temp_curve(doy1_2020_up-365));






%remove all values between winter-winter pairs (i.e. change in thaw=0, both
%scenes are totally frozen)
dh_2019a_up=[];
is2_Ntempa_up_2019=[];
Btemp_up_2019=[];
for i=1:length(dh_2019_up)
    if ~is2_Ntemp_up_2019(i)==0
        is2_Ntempa_up_2019=[is2_Ntempa_up_2019,is2_Ntemp_up_2019(i)];
        dh_2019a_up=[dh_2019a_up,dh_2019_up(i)];
        Btemp_up_2019=[Btemp_up_2019,doy2_2019_up(i)-doy1_2019_up(i)];
    end
end
dh_2020a_up=[];
is2_Ntempa_up_2020=[];
Btemp_up_2020=[];
for i=1:length(dh_2020_up)
    if ~is2_Ntemp_up_2020(i)==0
        is2_Ntempa_up_2020=[is2_Ntempa_up_2020,is2_Ntemp_up_2020(i)];
        dh_2020a_up=[dh_2020a_up,dh_2020_up(i)];
        Btemp_up_2020=[Btemp_up_2020,doy2_2020_up(i)-doy1_2020_up(i)];
    end
end



%bin along the y-axis (i.e. reduce the statistical spread).
for i=1:length(is2_Ntempa_up_2019)
    val=is2_Ntempa_up_2019(i);
    [v inds]=find(is2_Ntempa_up_2019==val);
    mean_up_2019(i)=mean(dh_2019a_up(inds));
end
for i=1:length(is2_Ntempa_up_2020)
    val=is2_Ntempa_up_2020(i);
    [v inds]=find(is2_Ntempa_up_2020==val);
    mean_up_2020(i)=mean(dh_2020a_up(inds));
end




Nbin=51;
lambda=.1;
pluck_val=85;
cmap=parula(9);
cmap2=[1,1,1;cmap];




%linear regression 0.6 max thaw indices, no y-axis binning
x_up=is2_Ntempa_up_2019(:);
n=length(x_up);
d=find(abs(is2_Ntempa_up_2019)<prctile(x_up,pluck_val));
d_is=is2_Ntempa_up_2019(d);
d_dh=dh_2019a_up(d);
d_mean=mean_up_2019(d);
%linear regression
x_up=d_is(:);
y_up=d_dh(:);
%y_up=d_mean(:);
[p_up,S]=polyfit(x_up,y_up,1);
[f_up,sigma] = polyval(p_up,x_up,S);
se0=sqrt(sum((f_up-y_up).^2)/(n-2));
n=length(x_up);
[ b sigma2_x x_est y_est stats] = deming(x_up,y_up,lambda);
se1=stats.s_e;
rho=corrcoef(x_up,y_up);
x_up_2019=x_up;
y_up_2019=y_up;
f_up_2019=f_up;
x_est_2019=x_est;
y_est_2019=y_est;
se1_2019=se1;
rho_2019=rho;
se_2019=2.4041;
%linear regression 0.6 max thaw indices, no y-axis binning
x_up=is2_Ntempa_up_2020(:);
n=length(x_up);
d=find(abs(is2_Ntempa_up_2020)<prctile(x_up,pluck_val));
d_is=is2_Ntempa_up_2020(d);
d_dh=dh_2020a_up(d);
d_mean=mean_up_2020(d);
%linear regression
x_up=d_is(:);
y_up=d_dh(:);
%y_up=d_mean(:);
[p_up,S]=polyfit(x_up,y_up,1);
[f_up,sigma] = polyval(p_up,x_up,S);
se0=sqrt(sum((f_up-y_up).^2)/(n-2));
n=length(x_up);
[ b sigma2_x x_est y_est stats] = deming(x_up,y_up,lambda);
se1=stats.s_e;
rho=corrcoef(x_up,y_up);
x_up_2020=x_up;
y_up_2020=y_up;
f_up_2020=f_up;
x_est_2020=x_est;
y_est_2020=y_est;
se1_2020=se1;
rho_2020=rho;
se_2020=2.3646;



figure;
subplot(2,1,1)
%dx=max(abs(x_up)); %set range of xovers
sc2=hist3([x_up_2019(:),y_up_2019(:)],'Edges',[{-1:(2/(Nbin-1)):1}, {-50:(100/(Nbin-1)):50}],'EdgeColor','w','LineStyle','none'); yy1=linspace(-50, 50, Nbin); xx1=linspace(-1, 1, Nbin);
pcolor(yy1,xx1,sc2);
hcb=colorbar;
colormap(cmap2)
%colormap('hot')
grid off
hold on
plot(f_up_2019,x_up_2019,'m-','LineWidth',2)
hold on
plot(y_est_2019,x_est_2019,'k-','LineWidth',2)
ylim([-0.5 0.5])
dim = [0.2 0.5 0.3 0.3];
string1=strcat('\color{magenta} linear regression: \rho = ',string(rho_2019(1,2)),' standard error = ',string(se_2019));
string2=strcat('\color{black} Deming regression: \rho = ',string(rho_2019(1,2)),' standard error = ',string(se1_2019));
str = {string1,string2};
annotation('textbox',dim,'String',str,'FitBoxToText','on','BackgroundColor','w','FontSize',14);
%,'m--',f_up,x_up-2*sigma,'m--','LineWidth',5)
set(get(hcb,'label'),'String','Number of Crossovers')
xlabel('dh (cm)')
ylabel('Change in Normalized Accumulated Degree Days')
%title({dataset,'Case 2: Throw out |NADD|>0.6, No Smoothing'})
set(gca,'FontSize',14)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
subplot(2,1,2)
%dx=max(abs(x_up)); %set range of xovers
sc2=hist3([x_up_2020(:),y_up_2020(:)],'Edges',[{-1:(2/(Nbin-1)):1}, {-50:(100/(Nbin-1)):50}],'EdgeColor','w','LineStyle','none'); yy1=linspace(-50, 50, Nbin); xx1=linspace(-1, 1, Nbin);
pcolor(yy1,xx1,sc2);
hcb=colorbar;
colormap(cmap2)
%colormap('hot')
grid off
hold on
plot(f_up_2020,x_up_2020,'m-','LineWidth',2)
hold on
plot(y_est_2020,x_est_2020,'k-','LineWidth',2)
ylim([-0.5 0.5])
dim = [0.2 0.5 0.3 0.3];
string1=strcat('\color{magenta} linear regression: \rho = ',string(rho_2020(1,2)),' standard error = ',string(se_2020));
string2=strcat('\color{black} Deming regression: \rho = ',string(rho_2020(1,2)),' standard error = ',string(se1_2020));
str = {string1,string2};
annotation('textbox',dim,'String',str,'FitBoxToText','on','BackgroundColor','w','FontSize',14);
%,'m--',f_up,x_up-2*sigma,'m--','LineWidth',5)
set(get(hcb,'label'),'String','Number of Crossovers')
xlabel('dh (cm)')
ylabel('Change in Normalized Accumulated Degree Days')
%title({dataset,'Case 2: Throw out |NADD|>0.6, No Smoothing'})
set(gca,'FontSize',14)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);







%linear regression test statistics
fitlm(x_up_2019,y_up_2019)
fitlm(x_up_2020,y_up_2020)

