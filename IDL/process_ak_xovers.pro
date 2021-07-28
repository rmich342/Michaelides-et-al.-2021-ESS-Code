pro process_ak_xovers, xover, tvector, x, ndata, model_size, model_smooth, resid_norm, $
	on=on, off=off, year=year, sigmacut=sigmacut, weights=weights, $ 
	lambda_size=lambda_size, lambda_smooth=lambda_smooth
;;
;; year is a required keyword specifying one or more target years in form [year1, year2, ...]
;;


;;PARAMETERS
snow_off_date = 2019 + 157/365D
snow_back_date = 2019 + 275/365D

;;SUBSET DATA
if keyword_set(off) then begin
	print, 'SNOW OFF XOVERS ONLY'
	iuse = where(xover.t1 ge snow_off_date and xover.t2 le snow_back_date)
	xover = xover(iuse)
endif
if keyword_set(on) then begin
	print, 'SNOW ON XOVERS ONLY'
	iuse = where(xover.t2 le snow_off_date)
	xover = xover(iuse)
endif
if n_elements(year) eq 1 then begin
	iuse = where(xover.t1 ge year and xover.t2 le year + 1)
	xover = xover(iuse)
endif

if keyword_set(sigmacut) then begin
	iuse = where(xover.sigma_dh le sigmacut AND abs(xover.dhdx_1) le sigmacut AND abs(xover.dhdx_2) le sigmacut)
	;iuse = where(xover.sigma_dh le sigmacut)
	;iuse = where(abs(xover.dhdx_1) le sigmacut and abs(xover.dhdx_2) le sigmacut)
	;iuse = where(xover.sigma_h1 le sigmacut and xover.sigma_h2 le sigmacut)
	xover = xover(iuse)
endif
dt_days = xover.dt_days
t1 = xover.t1
t2 = xover.t2
julday1 = xover.julday1
julday2 = xover.julday2
dh = xover.dh
dhdt = xover.dhdt
sigma_dh = xover.sigma_dh
ndata = n_elements(dt_days)

;;RECAST MULTI-YEAR DATA ONTO SAME ANNUAL TIMELINE
t1 = floor((t1 - floor(t1))*365D)
t2 = floor((t2 - floor(t2))*366D)
min_t = min(t1)
max_t = max(t2)
nmodel = max_t - min_t + 1									;number of elements in time series vector
if n_elements(year) eq 1 then $
	tvector = julday(1,1,year) + findgen(nmodel) + min_t else $
	tvector = julday(1,1,year(0)) + findgen(nmodel) + min_t

;;INVERSE PROBLEM SETUP
;;Want to solve Gx = d, where d is our vector of raw crossover height differences, x
;;is the target surface change time series, and G is the design matrix relating the two. 
;;G includes a smoothing matrix appended at the bottom.  The solution is:
;;  	Gx = d 
;;      or 
;;  	G+ = VS+U' (from SVD)
;;  	x = VS+U'd
;;
;;    	M'Mx = M'd (from Cholesky)
;;    	or M'WM = M'Wd (weighted)
;;    	where W is diagonal matrix of 1/variance

;;PARAMETERS OF INVERSION
if ~keyword_set(lambda_size) then lambda_size = 0.01
if ~keyword_set(lambda_smooth) then lambda_smooth = 0.00
weights = 0

;;WEIGHTS
if keyword_set(weights) then begin
	inverse_variance_dh = 1D/(sigma_dh)^2D
	W = diag_matrix(inverse_variance_dh)
endif

;;CREATE OBSERVATION DATA VECTOR (data = d)
data = dh
if keyword_set(weights) then begin
	data = matrix_multiply(data, W)				;W##data
endif

;;CREATE JACOBIAN MATRIX (jmatrix) BY ROW FROM JACOBIAN INPUT STRUCTURE
jmatrix = fltarr(nmodel, ndata)
for i = 0, ndata - 1 do begin
	jmatrix(t1(i) - min_t, i) = -1D
	jmatrix(t2(i) - min_t, i) = 1
endfor
if keyword_set(weights) then begin
	jmatrix = matrix_multiply(jmatrix, W)		;W##jmatrix
endif

;;NORMALIZE USING MINIMUM MODEL SIZE
if lambda_size gt 0.0 then begin
	ones = fltarr(nmodel) + 1.0*lambda_size
	sizecoda = diag_matrix(ones)
	;sizecoda(0) = 100*lambda_size				;this pins model to zero at its first value
	jmatrix = [[jmatrix], [sizecoda]]
	data = [reform(data), fltarr(nmodel)]		;for size normalization
endif	

;;NORMALIZE USING MODEL SMOOTHNESS
if lambda_smooth gt 0.0 then begin
	smoothcoda = fltarr(nmodel, nmodel)
	smoothcoda(0, 0) = 0.0
	for i = 1, nmodel - 1 do begin
		smoothcoda(i, i) = 1.0*lambda_smooth
		smoothcoda(i - 1, i) = -1.0*lambda_smooth	
	endfor
	jmatrix = [[jmatrix], [smoothcoda]]
	data = [reform(data), fltarr(nmodel)]		;for size normalization
endif

;;MAKE MATRICES FOR INVERSION
m = transpose(jmatrix)##jmatrix
rhs = transpose(jmatrix)##data

;;INVERT USING SVD OF JACOBIAN 
la_svd, M, sv, U, V
nsv = n_elements(sv)							;number of singular values 
svi = dblarr(nsv, nsv)							;make matrix for inverse singular values
svi(indgen(nsv)*(nsv + 1)) = 1/sv				;set diagonal elements of svi to be 1/sv          
x = V ## svi ## transpose(U) ## rhs				;solution vector b
model_size = norm(x)
model_smooth = norm(x(1:nmodel - 1) - x(0:nmodel - 2))

;;GET MODELED DH
for i = 0, ndata - 1 do begin
	i1 = where(tvector eq t1(i))
	i2 = where(tvector eq t2(i))
	xover(i).dh_model = x(i2) - x(i1)
endfor
resid_norm = norm(xover.dh - xover.dh_model)

;;PRINT
if keyword_set(verbose) then begin
	print, 'Plotting only Year = ' + int2str(year)
	print, 'Number of xovers: ', ndata
	print, 'Model size: ', model_size
	print, 'Model smooth: ', model_smooth
	print, 'Model residual: ', resid_norm
endif

end

