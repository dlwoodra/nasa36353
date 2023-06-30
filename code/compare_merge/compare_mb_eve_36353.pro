pro compare_mb_eve_36353, wave, irradiance, relEerr

window,0,xs=10,ys=10
wdelete
!p.color=0 & !p.background='ffffff'x & !p.charsize=1.5

@config36353

workingdir=file_dirname(routine_filepath()) ; in code/compare_merge
datadir = workingdir+'/../../data/'

restore,'rocket36'+numberstr+'_megsb_irr.sav' ; use spectra

; recall solar motion is detected at image #50 forward
;print,'do not use images after #38 at ',spectra[50].time
;; Tom said
;;  Apogee = 293.3 km at T+279.15 sec
;;  Good ADCS data from T+123 to T+360 sec
;;   T+133 sec is first time for good 10-sec integration
;;   The pitch/yaw after T+360 is very good, but roll axis is not good
;;   Recommend cal data not use anything after T+360 sec
; Tom said there is a roll that occurs near the end of the flight for SAM/ESP
;
;solar = lindgen(21)+14L ; indices of solar measurements
;solar = lindgen(30)+32L ; indices of solar measurements

; 6/30/23 restrict time range to good spectra only
;solar = solaridx_a[where(spectra[solaridx_a].time gt 150 and spectra[solaridx_a].time lt 400)]
apogeesec = float(apogeesecstr)
solar = solaridx_b[where(spectra[solaridx_b].time gt apogeesec-41 and $
                         spectra[solaridx_b].time lt apogeesec+41)]
fill_value = 1e-8
upperlimit = 1e-2

for i=0,n_elements(solar)-1 do begin
  bad=where(spectra[solar[i]].sp lt 1e-8 or spectra[solar[i]].sp gt upperlimit,n_bad)
  if n_bad gt 0 then spectra[solar[i]].sp[bad] = fill_value
endfor

;stop
;spectra.sp[3660:*]=1e-8 ; remove large spike beyond 106.2 nm
;spectra_cps.sp[3660:*] = .01 & spectra_cps.sp=spectra_cps.sp>.01

new = remove_atm_absorption_36353( spectra[solar].time, spectra[solar[0]].w, $
                             spectra[solar].sp )

;  0.96864051 ; June 18 (169) at 19:00 UTC
;au = 0.9839 ; 2010123
; au now comes from config36353
new /= au

;stop
; replace data with atm corrected values
spectra[solar].sp = transpose(new)

solarlo=solar[n_elements(solar)/2 - 5]
solarhi=solarlo+10
; avg just 20-30 (20 is T+202 sec, using 22 as T+222sec to be safer)
rsp=mean(spectra[solarlo:solarhi].sp,dim=2) > (sens3700 * 0.1) ; 0.1 DN
rwave=spectra[0].w

wave=rwave
irradiance=rsp

plot,rwave,rsp,/ylog,yr=[1e-8,1e-2],xr=[30,110],xs=1,ps=10,tit=strtrim(theyd,2)+' 36.'+numberstr,xtit='Wavelength (nm)',ytit='W/m^2/nm'
for i=-1,5 do begin
   sf = 10.^i
   tmp=sens3700*sf
   oplot,rwave,tmp,co='fe0000'x
   xyouts,80,tmp[2900],strtrim(string(sf,form='(f10.1)'),2)+'*S',/data,co='fe0000'x
endfor
stop

; L2b
;sftp evescl1:/evenetapp/testing/data/level2b/2013/294/EVS_L2B_2013294_007_01.fit.gz .

datafile = file_search(datadir+'/EVS_L2B_'+strtrim(theyd,2)+'*fit.gz',count=count)

d=eve_read_whole_fits(datafile[0])
; avg just 15 minutes
;evesp=mean(d.spectrum[1140:1155].irradiance,dim=2) ; filter?
; find low and high time range from 18:32-18:45
junk=min(abs(d.spectrum.sod - 18L*60*60 + 32L*60), lo)
junk=min(abs(d.spectrum.sod - (18L*60*60 + 45L*60)), hi)
lo=lo[0]
hi=hi[0]

;evesp=mean((smooth(d.spectrum.irradiance,[1,6],/edge_trun))[*,1140:1155],dim=2)
evesp=mean((smooth(d.spectrum.irradiance,[1,6],/edge_trun))[*,lo:hi],dim=2)
evewave=d.spectrummeta.wavelength
eveyd=d.spectrum[0].yyyydoy

;plot,rwave,rsp,/ylog,yr=[1e-7,1e-2],xr=[30,105],xs=1,ps=10,tit=+strtrim(theyd,2)+' 36.258',xtit='Wavelength (nm)',ytit='W/m^2/nm'
;oplot,evewave,evesp,co='fe0000'x,ps=10
;xyouts,31,6e-3,/data,'EVE L2B 15-min avg V6',color='fe0000'x
;xyouts,31,4e-3,/data,'Rocket 36.258'

file='rkt_36'+numberstr+'_megsb_full_compare.png'
width=800 & height=800
p = plot( rwave, rsp, $
  dimensions=[width,height],/stairstep,xstyle=1,xrange=[33,105],$
  /ylog,yrange=[1e-7,1e-2],$
  xtitle='Wavelength (nm)',ytitle='Irradiance (W/m^2/nm)',$
  title=strtrim(theyd,2)+' NASA 36.'+numberstr+' MEGS-B')
gd=where(evewave gt 33 and evewave lt 106.5)
p = plot( evewave[gd],evesp[gd], /overplot, /stairstep, color='blue' )
p = text(/norm,.15,.75,'EVE L2B 15-min avg '+file_basename(datafile[0]),color='blue')
p = text(/norm,.15,.8,'36.'+numberstr)
p.save,file,width=width,height=height
stop
p.close

width=800 & height=600
!p.multi=[0,1,2]
for i=30,100,5 do begin
   xlo=float(i) & xhi = xlo+8.
   inside=where(rwave gt xlo and rwave lt xhi)
   eveinside=where(evewave gt xlo and evewave lt xhi)
   yrange=[min(rsp[inside])>1e-7,max(rsp[inside])>max(evesp[eveinside])]
   plot,rwave,rsp,/ylog,yr=yrange,xr=[xlo,xhi],xs=1,ps=10,tit=+strtrim(theyd,2)+' 36.'+numberstr,$
        xtit='Wavelength (nm)',ytit='Irradiance (W/m^2/nm)'
   oplot,evewave,evesp,co='fe0000'x,ps=10
   xyouts,i+.2,.4*10.^!y.crange[1],/data,'EVE L2B 15-min avg',color='fe0000'x
   xyouts,i+.2,.1*10.^!y.crange[1],/data,'Rocket 36.'+numberstr

   ratio = rsp[inside]/evesp[eveinside] ; same wavelength sampling
   plot,rwave[inside],ratio,xr=[xlo,xhi],xs=1,yr=[0,2],ys=1,$
     title='Rocket/EVE',ytit='Ratio R/E',xtit='Wavelength (nm)',ps=10
   oplot,!x.crange,[1,1],lines=1,co='aa00'x
   ;stop
   
   ;
   ; show same stuff with 1 angstrom smooth
   ;
   smk = 5
   smrsp = smooth(rsp,smk,/edge_trunc)
   smevesp = smooth(evesp,smk,/edge_trunc)
   plot,rwave,smrsp,/ylog,yr=yrange,xr=[xlo,xhi],xs=1,ps=10,tit='0.1nm smooth '+strtrim(theyd,2)+' 36.'+numberstr,$
        xtit='Wavelength (nm)',ytit='Smoothed Irradiance (W/m^2/nm)'
   oplot,evewave,smevesp,co='fe0000'x,ps=10
   xyouts,i+.2,.4*10.^!y.crange[1],/data,'EVE L2B 15-min avg',color='fe0000'x
   xyouts,i+.2,.1*10.^!y.crange[1],/data,'Rocket 36.'+numberstr

   smratio = smrsp[inside]/smevesp[eveinside] ; same wavelength sampling
   plot,rwave[inside],smratio,xr=[xlo,xhi],xs=1,yr=[0,2],ys=1,$
     title='Rocket/EVE',ytit='Smoothed Ratio R/E',xtit='Wavelength (nm)',ps=10
   oplot,!x.crange,[1,1],lines=1,co='aa00'x
   ;stop
   
endfor


; errors
; assume senserr3700 is in units of sens3700, so ratio is fractional
; uncertainty
; since E = cps*S, then the uncertainty in E should be calculated as
; varE/E^2 = varcps / cps^2 + varS / S^2 + 2 covarcpsS / E
; assume covariance is negligible between S and cps (probably OK)
; but there is definitely correlation between bins because of the
; finite resolution

; variance in S
varS = senserr3700^2

; what is the variance in cps?
; assume it is stddev of irradiance measurements that contribute to rsp?
;solarlo = 42
;solarhi = 48
cps = mean(spectra_cps[solarlo:solarhi].sp,dim=2,/nan)>1e-4   ; averaged irradiance
varcps = variance(spectra_cps[solarlo:solarhi].sp,dim=2) ; > 0.1

varE = rsp^2 * (varcps/(cps^2) + vars/(sens3700^2))
sigmaE = sqrt(varE)

relEerr = sigmaE/rsp <1.

!p.multi=[0,1,2]
!p.color=0
plot,spectra[0].w, rsp,/ylog,ps=10,tit='36.'+numberstr+' MEGS-B Irradiance',ytit='W/(m^2 nm)',xr=[33,105],xs=1,xtit='Wavelength (nm)'
plot,spectra[0].w, relEerr,/ylog,tit='36.'+numberstr+' MEGS-B Relative Combined Uncertainty !7r!3=sqrt([!7r!3!Bcps!A2!n/cps!A2!N + !7r!3!BS!A2!N/S!A2!N])',ytit='Fractional',yr=[0.01,1],xr=[33,105],xs=1,xtit='Wavelength (nm)'
stop

plot,spectra[0].w, senserr3700/sens3700,/ylog,ps=10,tit='36.'+numberstr+' MEGS-B Relative Sensitivity Uncertainty',ytit='Fractional',yr=[0.001,1],xr=[33,105],xs=1,xtit='Wavelength (nm)'
plot,spectra[0].w, sqrt(varcps)/cps,/ylog,tit='36.'+numberstr+' MEGS-B Relative Precision',ytit='Fractional',yr=[0.005,1],ys=1,xr=[33,105],xs=1,xtit='Wavelength (nm)'

!p.multi=0
plot,spectra[0].w, relEerr,/ylog,tit='36.'+numberstr+' MEGS-B Relative Uncertainties',ytit='Fractional',yr=[0.001,1],xr=[33,105],xs=1,xtit='Wavelength (nm)',thick=2
oplot,spectra[0].w, senserr3700/sens3700,ps=10,co='fe'x
oplot,spectra[0].w, sqrt(varcps)/cps,ps=10,co='fe5050'x
xyouts,35,.002,/data,'Total Uncertainty'
xyouts,35,.004,/data,'Sensitivity Uncertainty',co='fe'x
xyouts,35,.008,/data,'cps Uncertainty',co='fe5050'x

prec=sqrt(varcps)/cps

savfile = datadir+'mb_36'+numberstr+'_irr_at_1au.sav'
save, file = savfile, wave, irradiance, relEerr, /compress
stop
return
end
