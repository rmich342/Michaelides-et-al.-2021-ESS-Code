pro plot_ak_xovers, xover, on=on, off=off, year=year, sigmacut=sigmacut


;;PARAMETERS
snow_off_date = 2019 + 157/365D
snow_back_date = 2019 + 275/365D
modis_2019 = [julday(6,6,2019), julday(9,17,2019)]
modis_2020 = [julday(6,11,2020), julday(9,18,2020)]
merra2_2019 = [julday(6,1,2019), julday(8,31,2019)]
merra2_2020 = [julday(6,1,2020), julday(8,31,2020)]
daymet_2019 = [julday(4,30,2019), julday(10,17,2019)]
daymet_2020 = [julday(4,30,2020), julday(10,17,2020)]

;;DERIVED PARAMETERS
xover_plot = xover
if keyword_set(off) then begin
	print, 'SNOW OFF XOVERS ONLY'
	iuse = where(xover_plot.t1 ge snow_off_date and xover_plot.t2 le snow_back_date)
	xover_plot = xover_plot(iuse)
endif
if keyword_set(on) then begin
	print, 'SNOW ON XOVERS ONLY'
	iuse = where(xover_plot.t2 le snow_off_date)
	xover_plot = xover_plot(iuse)
endif
if keyword_set(year) then begin
	print, 'Plotting only Year = ' + int2str(year)
	iuse = where(xover_plot.t1 ge year and xover_plot.t2 le year + 1)
	xover_plot = xover_plot(iuse)
endif
if keyword_set(sigmacut) then begin
	iuse = where(xover_plot.sigma_dh le sigmacut AND abs(xover_plot.dhdx_1) le sigmacut AND abs(xover_plot.dhdx_2) le sigmacut)
	;iuse = where(xover_plot.sigma_dh le sigmacut)
	;iuse = where(abs(xover_plot.dhdx_1) le sigmacut and abs(xover_plot.dhdx_2) le sigmacut)
	;iuse = where(xover_plot.sigma_h1 le sigmacut and xover_plot.sigma_h2 le sigmacut)
	xover_plot = xover_plot(iuse)
endif
isort = reverse(sort(xover_plot.dt_days))
dt_days = xover_plot(isort).dt_days
t1 = xover_plot(isort).t1
t2 = xover_plot(isort).t2
julday1 = xover_plot(isort).julday1
julday2 = xover_plot(isort).julday2
dh = xover_plot(isort).dh
dhdt = xover_plot(isort).dhdt
n = n_elements(dt_days)
print, n

;;SETUP PLOTTING
set_plot, 'ps'
!p.multi = [0, 1, 1]
!p.font = 0
device, decomposed = 0
if keyword_set(year) then begin
	tstart = year
	tend = year + 1
	tstart = julday(1,1,year)
	tend = julday(1,1,year + 1)
	xtitletext = int2str(year, 4)
endif else begin
	t1 = t1 - floor(t1)
	t2 = t2 - floor(t2)
	tstart = 0
	tend = 1
	xtitletext = 'fractional year'
endelse

;;PLOT CROSSOVER LENGTHS
xcm = 5.8
ycm = 5.0
device, filename = 'dt.eps', /color, bits = 8, /encapsulated, xsize = xcm, ysize = ycm, $
	/helvetica
loadct, 0, /silent
plot,dt_days(sort(dt_days)), psym=symcat(15), /yst,/xst,symsize=0.5,chars=1,chart=2,yr=[0,250],$
	ytit='crossover dt (days)',xtit='crossover index' 
oplot, dt_days(sort(dt_days))
device, /close_file

;;PLOT DH/DT
xcm = 9
ycm = 8
device, filename = int2str(year,4) + '_dhdt_time.eps', /color, bits = 8, /encapsulated, xsize = xcm, ysize = ycm, $
	/Helvetica, font_size = 8
loadct, 0, /silent
plot, julday1, dhdt, /nodata, $
	xr=[tstart, tend], xminor=3, xticku = 'month', $
	yr=[-2,1], /yst, $
	ytit = 'crossover dh/dt (cm/day)', $
	xtit = xtitletext
loadct,0,/silent
print,'Min/Max Timespan:', min(dt_days), max(dt_days)
for j = 0, n - 1 do $
	oplot, [julday1(j), julday2(j)],[dhdt(j), dhdt(j)], thick = 1, color = (220 - dt_days(j))*200D/220D
loadct,0,/silent
oplot, [tstart, tend], [0, 0], color = 0, linestyle = 1

;;ANNOTATE
;loadct, 3, /silent
;icolor = 161
loadct, 0, /silent
icolor = 0
shift = -0.011
if keyword_set(year) then begin
	case year of
		2019: begin
			modis = modis_2019
			merra2 = merra2_2019
			daymet = daymet_2019
			oplot, [daymet(0),daymet(1)], [-1.85, -1.85]+shift, linestyle = 0, thick = 2, color = icolor
			oplot, [daymet(0),daymet(0)], [-1.85, -1.78]+shift, linestyle = 0, thick = 2, color = icolor
			oplot, [daymet(1),daymet(1)], [-1.85, -1.78]+shift, linestyle = 0, thick = 2, color = icolor
			xyouts, total(daymet)/2D, -1.81+shift, 'DAYMET', color = icolor, chars = 0.9, charthick=2, align=0.5
			oplot, [modis(0), modis(1)], [-1.6, -1.6]+shift, linestyle = 0, thick = 2, color = icolor
			oplot, [modis(0), modis(0)], [-1.6, -1.53]+shift, linestyle = 0, thick = 2, color = icolor
			oplot, [modis(1), modis(1)], [-1.6, -1.53]+shift, linestyle = 0, thick = 2, color = icolor
			xyouts, total(modis)/2D, -1.56+shift, 'MODIS', color = icolor, chars = 0.9, charthick=2, align=0.5
			oplot, [merra2(0), merra2(1)], [-1.35, -1.35]+shift, linestyle = 0, thick = 2, color = icolor
			oplot, [merra2(0), merra2(0)], [-1.35, -1.28]+shift, linestyle = 0, thick = 2, color = icolor
			oplot, [merra2(1), merra2(1)], [-1.35, -1.28]+shift, linestyle = 0, thick = 2, color = icolor
			xyouts, total(merra2)/2D, -1.31+shift, 'MERRA-2', color = icolor, chars = 0.9, charthick=2, align=0.5
		end
		2020: begin
			modis = modis_2020
			merra2 = merra2_2020
			daymet = daymet_2020
			oplot, [daymet(0),daymet(1)], [-1.85, -1.85]+shift, linestyle = 0, thick = 2, color = icolor
			oplot, [daymet(0),daymet(0)], [-1.85, -1.78]+shift, linestyle = 0, thick = 2, color = icolor
			oplot, [daymet(1),daymet(1)], [-1.85, -1.78]+shift, linestyle = 0, thick = 2, color = icolor
			xyouts, total(daymet)/2D, -1.81+shift, 'DAYMET', color = icolor, chars = 0.9, charthick=2, align=0.5
			oplot, [merra2(0), merra2(1)], [-1.35, -1.35]+shift, linestyle = 0, thick = 2, color = icolor
			oplot, [merra2(0), merra2(0)], [-1.35, -1.28]+shift, linestyle = 0, thick = 2, color = icolor
			oplot, [merra2(1), merra2(1)], [-1.35, -1.28]+shift, linestyle = 0, thick = 2, color = icolor
			xyouts, total(merra2)/2D, -1.31+shift, 'MERRA-2', color = icolor, chars = 0.9, charthick=2, align=0.5
			oplot, [modis(0), modis(1)], [-1.6, -1.6]+shift, linestyle = 0, thick = 2, color = icolor
			oplot, [modis(0), modis(0)], [-1.6, -1.53]+shift, linestyle = 0, thick = 2, color = icolor
			oplot, [modis(1), modis(1)], [-1.6, -1.53]+shift, linestyle = 0, thick = 2, color = icolor
			xyouts, total(modis)/2D, -1.56+shift, 'MODIS', color = icolor, chars = 0.9, charthick=2, align=0.5
		end
	endcase	
endif
xyouts, julday(11,26,year), 0.775, 'b.)', chart=3, align=0, chars = 1.25

device, /close_file

;;PLOT 2
if keyword_set(detail) then begin
	xcm = 12D
	ycm = 12D
	device, filename = 'dhdt_time_detail.eps', /color, bits = 8, /encapsulated, xsize = xcm, ysize = ycm
	loadct, 0, /silent
	i = 200
	iplot = where(dt_days ge i - 25 and dt_days le i, nplot)
	plot, t1, dhdt, chars=1.0, xminor=2, xr=[tstart, tend], /nodata, yr=[-2,1], /yst, chart=2, $
		ytit = 'crossover dh/dt (cm/day)', xtit = xtitletext, tit = 'Subset: ' + int2str(i - 25) + '-' + int2str(i) + ' days'
	if nplot gt 0 then begin
		t1_plot = t1(iplot)
		t2_plot = t2(iplot)
		dhdt_plot = dhdt(iplot)
		dh_plot = dh(iplot)
		for j = 0, nplot - 1 do $
			oplot, [t1_plot(j), t2_plot(j)],[0, dhdt_plot(j)], color = i*250D/225D, thick=2
	endif
	loadct,0,/silent
	oplot, [tstart, tend], [0, 0], color = 0, linestyle = 2
	device, /close_file
endif

;;RESET PLOT OPTIONS
set_plot, 'x'
!P.multi = 0

end

