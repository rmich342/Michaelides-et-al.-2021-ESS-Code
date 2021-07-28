pro read_ak_xovers, xover, on=on, off=off, year=year, sigmacut=sigmacut


;;PARAMETERS
snow_off_date = 2019 + 157/365D
snow_back_date = 2019 + 275/365D

;;READ DATA
;:This is for the dh/dt figures for each year
filename = '../data/xovers_two_seasons_low_slopes_new_filtered.csv'
;;This is for the dh model figures
;filename = '../data/xovers_full_period_low_slopes_new_filtered.csv'
temp = read_csv(filename)
ndata = n_elements(temp.field01)
;filenames = file_search('../data/*.csv')
;ifile = 3
;temp = read_csv(filenames(ifile))

;;MAKE AND POPULATE DATA STRUCTURE
xover = replicate({ $
	lon:0D, $
	lat:0D, $
	dh:0D, $
	sigma_dh:0D, $
	h1_interp:0D, $
	h2_interp:0D, $
	t1:0D, $
	julday1:0D, $
	t2:0D, $
	julday2:0D, $
	dt_year:0D, $
	dt_days:0D, $
	sigma_h1:0D, $
	sigma_h2:0D, $
	dhdx_1:0D, $
	dhdx_2:0D, $
	n_photons_1:0D, $
	n_photons_2:0D, $
	dhdt:0D, $
	dh_model:0D}, ndata)
	
xover.lon = temp.field01
xover.lat = temp.field02
xover.dh = temp.field03
xover.sigma_dh = temp.field04
xover.h1_interp = temp.field05
xover.h2_interp = temp.field06
xover.t1 = temp.field07
xover.t2 = temp.field08
xover.dt_year = temp.field09
xover.dt_days = temp.field10
xover.sigma_h1 = temp.field11
xover.sigma_h2 = temp.field12
xover.dhdx_1 = temp.field13
xover.dhdx_2 = temp.field14
xover.n_photons_1 = temp.field15
xover.n_photons_2 = temp.field16
xover.dhdt = xover.dh/xover.dt_days*100D 	;dh/dt in cm/day

;;CONVERT DATES
year_t1 = floor(xover.t1)
days_t1 = xover.t1 - year_t1
year_t2 = floor(xover.t2)
days_t2 = xover.t2 - year_t2
for i = 0, ndata - 1 do begin
    daysinyear, year_t1(i), monthdays, ndaysinyear
    doy = days_t1(i)*ndaysinyear
    doy2ymd, year_t1(i), doy, month, day
    xover(i).julday1 = julday(month, day, year_t1(i))
    doy = days_t2(i)*ndaysinyear
    doy2ymd, year_t2(i), doy, month, day
    xover(i).julday2 = julday(month, day, year_t2(i))
endfor

end

