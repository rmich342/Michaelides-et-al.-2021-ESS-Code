pro plot_ak_dhdt_model, xovers, output, year=year

;;PARAMETERS
l_size = [0.1]
l_smooth = 40; [29]*1; [1.4]*20
n_run = 100
fraction = 0.75

;;INITIAL RUN WITH NO SMOOTHING
process_ak_xovers, xovers, tvector, model, n_xover, model_size, model_smooth, resid_norm, $
	/readdata, lambda_size = l_size(0), lambda_smooth = 0.0, year = year
n_model = n_elements(model)
n_xovers = n_elements(xovers)

;;OUTPUT STRUCTURE
output = replicate({ $
	n_xover:0, $
	julday:fltarr(n_model), $
	model:fltarr(n_model), $
	l_size:0.0, $
	l_smooth:0.0, $
	model_size:0.0, $
	model_smooth:0.0, $
	resid_norm:0.0}, n_run + 2)

;;VALUES OF NO SMOOTHING RUN
output(0).n_xover = n_xover
output(0).julday = tvector
output(0).model = model
output(0).l_size = l_size(0)
output(0).l_smooth = 0.0
output(0).model_size = model_size
output(0).model_smooth = model_smooth
output(0).resid_norm = resid_norm

;;RUN WITH CHOSEN SMOOTHING
process_ak_xovers, xovers, tvector, model, n_xover, model_size, model_smooth, resid_norm, $
	lambda_size = l_size(0), lambda_smooth = l_smooth(0), year = year
output(1).n_xover = n_xover
output(1).julday = tvector
output(1).model = model
output(1).l_size = l_size(0)
output(1).l_smooth = 0.0
output(1).model_size = model_size
output(1).model_smooth = model_smooth
output(1).resid_norm = resid_norm

;;SETUP PLOTTING
set_plot, 'ps'
!p.font = 0
xcm = 9
ycm = 8
yearname = ''
if n_elements(year) eq 1 then $
	yearname = int2str(year) else $
	for i = 0, n_elements(year) - 1 do yearname = yearname + '+' + int2str(year,4)
device, filename = 'ak_xover_' + yearname + '.eps', /color, bits = 8, /encapsulated, xsize = xcm, ysize = ycm, $
	/helvetica, font_size = 8
device, decomposed = 0
loadct, 39, /silent
!p.multi = 0

;;DO PROCESSING
n_use = floor(fraction*n_xovers)
seed = 20.1343
;xovers_montecarlo = xovers
for i = 2, n_run do begin
	;;DO MONTE CARLO ERROR PROPAGATION
	;for j = 0, n_xovers - 1 do xovers_montecarlo(j).dh = xovers(j).dh + randomn(seed, 1)*xovers(j).sigma_dh
	;;DO JACKKNIFE VARIANCE CALCULATION
	xovers_montecarlo = RANDOM_SAMPLE(seed, xovers, n_use)	
	process_ak_xovers, xovers_montecarlo, tvector, model, n_xover, model_size, model_smooth, resid_norm, $
		lambda_size = l_size(0), lambda_smooth = l_smooth(0), year = year
	output(i).n_xover = n_xover
	output(i).julday = tvector
	output(i).model = model
	output(i).model_size = model_size
	output(i).model_smooth = model_smooth
	output(i).resid_norm = resid_norm
endfor
model_sigma = stddev(output(2:n_run - 1).model, dimension = 2)

;;PLOT
xrange = [julday(1,0,year(0)), julday(1,1,year(0)+1)]
plot, output(0).julday, output(0).model*100D, $
	;position = [0.07, 0.15, 0.33, 0.95], /normal, $
	xr=xrange, /xst, xticku = 'month', $
	yr=[-0.4,0.4]*100D, /yst, $
	xtit=yearname, ytit = 'dh (cm)', chart = 2 
	;, tit='All Crossovers, 2019-2020 (N = 2201)'
oplot, xrange, [0, 0], linestyle = 1, thick = 2
oplot, output(1).julday, 100D*output(1).model, color = 250, thick = 3
oplot, output(1).julday, 100D*(output(1).model + model_sigma), color = 210, thick = 2
oplot, output(1).julday, 100D*(output(1).model - model_sigma), color = 210, thick = 2
oplot, julday(1, 0, year(0)) + floor((xovers.t1 - floor(xovers.t1))*365D), fltarr(n_xovers), psym = symcat(15), symsize = 0.5, color=60
oplot, julday(1, 0, year(0)) + floor((xovers.t2 - floor(xovers.t2))*365D), fltarr(n_xovers), psym = symcat(15), symsize = 0.5, color=60
xyouts, julday(11,26,year), 34, 'b.)', chart=3, align=0, chars = 1.25

;;RESET PLOT OPTIONS
device, /close_file
set_plot, 'x'

;;WRITE DATA
data = [transpose(output(1).julday), transpose(output(1).model), transpose(model_sigma)]
formats = ["f10.2", "f10.4", "f10.4"]
write_generic_data_column, data, formats, outfile='ak_xover_model_' + int2str(year, 4) + '.txt'

end

