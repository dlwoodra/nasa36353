function s12join, s1, s2, idx=idx

if size(idx,/type) eq 0 then idx=1424L

s = [s1[0:idx-1], s2[idx:*]]

return,s
end


;+
; This procedure makes a bunch of plots and returns the MEGS-A merged
; spectrum irradiance (rsp) and relative uncertainty (relEerr)
; 
; :Params:
;    wave: out, optional, type=dblarr
;      wavelength in nm
;   irradiance: out, optional, type=dblarr
;      irradiance 2d array of time and wavelength
;   relEerr: out, optional, type-dblarr
;      uncertainty based on stddev and sensitivity only, 2d array same
;      as irradiance
;
;-
pro compare_ma_eve_36353, wave, irradiance, relEerr

window,0,xs=10,ys=10
wdelete
!p.color=0 & !p.background='ffffff'x & !p.charsize=1.5

@config36353
;numberstr and theyd

workingdir=file_dirname(routine_filepath()) ; in code/compare_merge
datadir = workingdir+'/../../data/'

; not in data dir
restore,'rocket36'+numberstr+'_megsa_irr.sav' ; use spectra

; recall solar motion is detected at image #38 forward
;print,'do not use images after #50 at ',spectra[50].time
;; Tom said
;;  Apogee = 293.3 km at T+279.15 sec
;;  Good ADCS data from T+123 to T+360 sec
;;   T+133 sec is first time for good 10-sec integration
;;   The pitch/yaw after T+360 is very good, but roll axis is not good
;;   Recommend cal data not use anything after T+360 sec
;
;solar = lindgen(21)+14L ; indices of solar measurements (0-25)
;solar = lindgen(30)+32L
;solar = solaridx_a

; 6/30/23 restrict time range to good spectra only
;solar = solaridx_a[where(spectra[solaridx_a].time gt 150 and spectra[solaridx_a].time lt 400)]
apogeesec = float(apogeesecstr)
solar = solaridx_a[where(spectra[solaridx_a].time gt apogeesec-41 and $
                         spectra[solaridx_a].time lt apogeesec+41)]

; replace unrealistic values in solar spectra with a fill value
fill_value = 1e-8
upperlimit = 2e-2 ; guess at solar max for 30.4
;for i=3544,3599 do spectra[solar].sp2[i] = fill_value ;long wavelength slit 2
;for i=0,200 do spectra.sp2[i] = 1e-8 ; short wavelength slit 2

spectra[solar].sp2 = spectra[solar].sp2 > fill_value
spectra[solar].sp1 = spectra[solar].sp1 > fill_value

; remove bad high values
for i=0,n_elements(solar)-1 do begin
  ; slit 2
  bad=where(spectra[solar[i]].sp2 gt upperlimit)
  spectra[solar[i]].sp2[bad] = fill_value

  ; slit1
  bad =where(spectra[solar[i]].sp1 gt upperlimit)
  spectra[solar[i]].sp1[bad] = fill_value
endfor

new1 = remove_atm_absorption_36353( spectra[solar].time, spectra[solar[0]].w1, $
                              spectra[solar].sp1 )
new2 = remove_atm_absorption_36353( spectra[solar].time, spectra[solar[0]].w2, $
                              spectra[solar].sp2 )

;print,'ERROR: without the RADAR data no atmospheric absorption can be done properly'
;stop

;au = 0.96864051 ; June 18 (169) at 19:00 UTC
;theyd = 2013294
;au = 1.00934 ; Oct 21, 2013
;au = 0.9839
;au = 0.98589 ; from lisird lasp_vsop87_1au_correction_PT1M 9/9/21 17:30 UTC
; call this interactively to get 1-AU factor from LISIRD
; s=get_lisird_data(dataset='lasp_vsop87_1au_correction_PT1M',mintime='2021-09-09T17:30:00',maxtime='2021-09-09T17:35:00',/jd)
; earth is far from sun, so need to increase irradiance

; au is now defined in config36353.pro
new1 /= au
new2 /= au

;; apply ma1deg and ma2deg MEGS-A degradation
;restore,'rkt36336_degradation.sav' ; ma1deg, ma2deg, mbdeg
;; need to interpolate from 2048 (reverse wavelength order) to the 3600
;; output bins
;stop

;stop
; replace data with atm corrected values
spectra[solar].sp1 = transpose(new1)
spectra[solar].sp2 = transpose(new2)

; choose images to average
; use timeseries figure of 36.8 from read_ma_36290
;idx=[22,23,24,25,26,27,28,29,30,31,32,33]
;solaridx=[42,43,44,45,46,47,48]
solaridx = solaridx_a
; TODO: use both solaridx_a and solaridx_b to find the best common indices

rsp1=mean(spectra[solaridx].sp1,dim=2) > (a1sens3600*.1) ; .1 DN
rsp2=mean(spectra[solaridx].sp2,dim=2)> (a2sens3600*.1) ; .1 DN

rwave1=spectra[0].w1
rwave2=spectra[0].w2 ; same as w1
; remove bad spike in rsp2 beyond 38.42nm
bad=where(rwave2 gt 38.42,n_bad)
if n_bad gt 0 then rsp2[bad]=0.d

wave=rwave1 ; since rwave1 = rwave2

; test sharpening slit 1 to match slit 2
ksharp=[-1.d,-1,6,-1,-1] & ksharp/=total(ksharp)
testrsp1 = convol(rsp1,ksharp)

; test smoothing slit 2 to match slit 1
kdull=[1,2,3,4,5,4,3,2,1] & kdull/=total(kdull)
testrsp2 = convol(rsp2,kdull)

; report 17.1 nm line comparisons
gd17=where(rwave1 ge 17 and rwave1 lt 17.2)
print,'17.1nm slit 1=',total(rsp1[gd17])
print,'17.1nm slit 2=',total(rsp2[gd17])
print,'17.1nm slit 1/slit2=',total(rsp1[gd17]) / total(rsp2[gd17])

gd174=where(rwave1 ge 17.35 and rwave1 lt 17.55)
print,'17.4nm slit 1/slit2=',total(rsp1[gd174]) / total(rsp2[gd174])

gd177=where(rwave1 ge 17.65 and rwave1 lt 17.85)
print,'17.7nm slit 1/slit2=',total(rsp1[gd177]) / total(rsp2[gd177])

;s12idx=1424 ; bin to switch from s1 to s2 ; 17.24 nm
;s12idx=1390 ; bin to switch from s1 to s2 ; 16.90 nm
s12idx=1193 ; bin to switch s1 to s2 ; 14.93 nm

rsp=s12join(rsp1,rsp2, idx=s12idx) ;[rsp1[0:s12idx-1],rsp2[s12idx:*]]
irradiance=rsp

r1wrange=[5,18]
r2wrange=[16,39]

a12sens3600=s12join(a1sens3600, a2sens3600, idx=s12idx) ;[a1sens3600[0:s12idx-1],a2sens3600[s12idx:*]]
a12senserr3600=s12join(a1senserr3600, a2senserr3600, idx=s12idx) ;[a1senserr3600[0:s12idx-1],a2senserr3600[s12idx:*]]

plot,rwave1[0:s12idx],rsp1[0:s12idx],/ylog,yr=[1e-8,1e-2],xr=[5,38],xs=1,ps=10,tit=strtrim(theyd,2)+' 36.'+numberstr,xtit='Wavelength (nm)',ytit='W/m^2/nm'
oplot,rwave2[s12idx:*],rsp2[s12idx:*],ps=10,co='fe0000'x
for i=-1,5 do begin
   sf = 10.^i
   tmp=a12sens3600*sf
   oplot,rwave1,tmp,co='fe'x
   xyouts,17,tmp[s12idx],strtrim(string(sf,form='(f10.1)'),2)+'*S',/data,co='fe'x
endfor
stop

; just slit 1 and slit 2 overlap
!p.multi=[0,1,2]
!p.charsize=1.5
;proposed_switchover2 = 14.7
proposed_switchover2 = 14.93 ; Tom approved this one 3/24/19 after OS issues
proposed_switchover  = 16.9 ; Tom approves this one 3/8/19
old_switchover       = 17.24
plot,rwave1,rsp1,/ylog,ps=10,yr=[1e-6,1e-3],xs=1,xr=[14,18], $
  xtit='Wavelength (nm)',ytit='Irradiance (W/m^2/nm)',tit='36.'+numberstr+' MEGS-A Full-Disk Solar Spectrum'
oplot,rwave2,rsp2,co='fe0000'x,ps=10
oplot,proposed_switchover*[1,1],[1e-6,1e-3],lines=2,co='fe'x
oplot,proposed_switchover2*[1,1],[1e-6,1e-3],lines=2,co='fe'x
oplot,old_switchover*[1,1],[1e-6,1e-3],lines=2
legend_dlw,['MEGS-A slit 1', 'MEGS-A slit 2'],textcolor=[0,'fe0000'x],charsize=2
plot,rwave1,rsp1/rsp2,ps=10,yr=[0.5,2], xr=!x.crange,xs=1,$
  xtit='Wavelength (nm)',ytit='S1 / S2',tit='36.'+numberstr+' MEGS-A Slit Ratio'
oplot,proposed_switchover*[1,1],!y.crange,lines=2,co='fe'x
oplot,proposed_switchover2*[1,1],!y.crange,lines=2,co='fe'x
oplot,old_switchover*[1,1],!y.crange,lines=2
oplot,!x.crange,[1,1],lines=1,co='aa00'x
stop

;url='http://lasp.colorado.edu/eve/data_access/evewebdata/products/level2b/' + $
;    string(theyd/1000,form='(i04)') + '/' + $
;    string(theyd mod 1000,form='(i03)') + $
;    '/EVS_L2B_'+strtrim(theyd,2)+'_006_02.fit.gz'
;localfilenamel2b='data/'+strtrim(theyd,2)+'.fit.gz'
;
;oUrl = obj_new('IDLnetURL')
;result = oURL->Get(url=url,filename=localfilenamel2b)
;obj_destroy, oUrl

; L2b
datafile = file_search(datadir+'EVS_L2B_'+strtrim(theyd,2)+'*fit.gz',count=count)
; have to use v7, no l2b v6 files available
;datafile = file_search(localfilenamel2b,count=count)

d=eve_read_whole_fits(datafile[0])

; avg just 15 minutes
;evesp=mean(d.spectrum[1140:1155].irradiance,dim=2) ; filter?
evesp=mean((smooth(d.spectrum.irradiance,[1,6],/edge_trun))[*,1140:1155],dim=2)
;stop
evewave = d.spectrummeta.wavelength

; lyman-apha composite shows 3.62 on June 18,2018
; equivalent dates in 2010 are May 16 and June 7,8 2010 and
; Jan 27,28 2011

; f10.7 was 71, and f107avg was 112(WRONG AVG)
; 70 on may 16 2010
; 72.8, 71.2 on June 7,8 2010 (Choose June 8, 2010159)
; 78.1 on both Jan 27,28, 2011

; options determined Mar 6, 2019:
;  <F10_30> is a 30-day mean, <F10_72> is 72-day mean
; YYYY DOY MM/DD  F10  <F10_30> <F10_72> lya
; 2010 152 06/01 74.8,   74.08    75.27  3.8 
; 2010 153 06/02 76.1,   74.1     75.29  3.75
; 2010 154 06/03 76.8,   74.2     75.39  3.72
; on 2018 06/18 rocket day
; F10=76.1, <F10_30>=74.13, <F10_72>=73.86, lyman-alpha=3.62
; so these 3 days are probably reasonable, but the bottom line
; is that there is no truly equivalent data for EVE.
; Lyman-alpha is not very close, and neither are the F10 means.

;; L3 2010159
;d=eve_read_whole_fits('data/EVE_L3_2010159_006_02.fit')
;evesp = d.data.sp_irradiance
;evewave = d.spectrummeta.wavelength
;eveyd=d.data.yyyydoy
; Instead just get the merged file and look at several days around here
url='https://lasp.colorado.edu/eve/data_access/evewebdata/products/merged/latest_EVE_L3_merged.ncdf'
localfilename=datadir+'/latest_EVE_L3_merged.ncdf'

oUrl = obj_new('IDLnetURL')
result = oURL->Get(url=url,filename=localfilename)
obj_destroy, oUrl

read_netcdf,localfilename,d
;eveyd=2010153
;eveyd=2013294
eveyd=2011174 ; 36.353 approx 100
;eveyd=2010139 ; matches 17.4-17.5 mean irradiance from EVE
; find this date in eve data to use instead of the one daily file
;stop
evesp = reform(d.mergeddata.sp_irradiance[*,where(d.mergeddata.yyyydoy eq eveyd)])
evewave = d.spectrummeta.wavelength ; 6-106, 5000 elements
;stop


file='rkt_36'+numberstr+'_megsa_full_compare.png'
width=800 & height=800
p = plot( rwave1, rsp, $
  dimensions=[width,height],/stairstep,xstyle=1,xrange=[6,37],$
  /ylog,yrange=[1e-6,1e-2],$
  xtitle='Wavelength (nm)',ytitle='Irradiance (W/m^2/nm)',$
  title=strtrim(theyd,2)+' NASA 36.'+numberstr+' MEGS-A')
gd=where(evewave gt 6 and evewave lt 37)
p = plot( evewave[gd],evesp[gd], /overplot, /stairstep, color='red' )
p = text(/norm,.15,.75,'EVE L3 daily avg V6 on '+strtrim(eveyd,2),color='red')
p = text(/norm,.15,.8,'36.'+numberstr)
p.save,file,width=width,height=height
;stop
p.close

width=800 & height=600
for i=5,30,5 do begin
  xlo = float(i) & xhi = xlo+8.
  inside=where(rwave1 gt xlo and rwave1 lt xhi)
  ;plot,rwave1,rsp,/ylog,yr=[1e-7>min(rsp[inside]),max(rsp[inside])],xr=[xlo,xhi],xs=1,ps=10,tit='2018169 36.'+numberstr,$
  ;      xtit='Wavelength (nm)',ytit='Irradiance (W/m^2/nm)'
  ;oplot,evewave,evesp,co='fe'x,ps=10
  ;xyouts,i+.2,6e-4,/data,'EVE L3 daily avg V6 on '+strtrim(eveyd,2),color='fe'x
  ;xyouts,i+.2,4e-4,/data,'Rocket 36.'+numberstr

  gd=where(evewave gt xlo and evewave lt xhi)
  yrange=[1e-7>(min(rsp[inside])<min(evesp[gd])),max(rsp[inside])>max(evesp[gd])]
  p = plot( rwave1, rsp, $
    dimensions=[width,height],/stairstep,xstyle=1,xrange=[xlo,xhi],$
    /ylog,yrange=yrange,$
    xtitle='Wavelength (nm)',ytitle='Irradiance (W/m^2/nm)',$
    title=strtrim(theyd,2)+' NASA 36.'+numberstr+' MEGS-A')
  p = plot( evewave[gd],evesp[gd], /overplot, /stairstep, color='red' )
  p = text(/norm,.15,.75,'EVE L3 daily avg V6 on '+strtrim(eveyd,2),color='red')
  p = text(/norm,.15,.8,'36.'+numberstr)
  file='rkt_36'+numberstr+'_megsa_compare_'+string(i,form='(i2.2)')+'.png'
  p.save,file,width=width,height=height
  ;stop
  p.close
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
varS = a12senserr3600^2

; what is the variance in cps?
; assume it is stddev of irradiance measurements that contribute to rsp?
cps1 = mean(spectra_cps[solaridx].sp1,dim=2)>1e-4 ; averaged irradiance
cps2 = mean(spectra_cps[solaridx].sp2,dim=2)>1e-4 ; averaged irradiance
cps = s12join(cps1,cps2, idx=s12idx) ;[cps1[0:s12idx-1],cps2[s12idx:*]]

varcps1 = variance(spectra_cps[solaridx].sp1,dim=2) > 0.1
varcps2 = variance(spectra_cps[solaridx].sp2,dim=2) > 0.1
varcps = s12join(varcps1, varcps2, idx=s12idx)

varE = rsp^2 * (varcps/(cps^2) + vars/(a12sens3600^2))
sigmaE = sqrt(varE)

relEerr = sigmaE/rsp <1.

!p.multi=[0,1,2]
plot,spectra[0].w1, rsp,/ylog,ps=10,tit='36.'+numberstr+' MEGS-A Irradiance',ytit='W/(m^2 nm)',xr=[5,38],xs=1,xtit='Wavelength (nm)'
plot,spectra[0].w1, relEerr,/ylog,tit='36.'+numberstr+' MEGS-A Relative Combined Uncertainty !7r!3=sqrt([!7r!3!Bcps!A2!n/cps!A2!N + !7r!3!BS!A2!N/S!A2!N])',ytit='Fractional',yr=[0.01,1],xr=[5,38],xs=1,xtit='Wavelength (nm)'
stop

plot,spectra[0].w1, a12senserr3600/a12sens3600,/ylog,ps=10,tit='36.'+numberstr+' MEGS-A Relative Sensitivity Uncertainty',ytit='Fractional',yr=[0.001,1],xr=[5,38],xs=1,xtit='Wavelength (nm)'
plot,spectra[0].w1, sqrt(varcps)/cps,/ylog,tit='36.'+numberstr+' MEGS-A Relative Precision',ytit='Fractional',yr=[0.005,1],ys=1,xr=[5,38],xs=1,xtit='Wavelength (nm)'

!p.multi=0
plot,spectra[0].w1, relEerr,/ylog,tit='36.'+numberstr+' MEGS-A Relative Uncertainties',ytit='Fractional',yr=[0.001,1],xr=[5,38],xs=1,xtit='Wavelength (nm)',thick=2
oplot,spectra[0].w1, a12senserr3600/a12sens3600,ps=10,co='fe'x
oplot,spectra[0].w1, sqrt(varcps)/cps,ps=10,co='fe5050'x
xyouts,6,.002,/data,'Total Uncertainty'
xyouts,6,.004,/data,'Sensitivity Uncertainty',co='fe'x
xyouts,6,.008,/data,'cps Uncertainty',co='fe5050'x

prec = sqrt(varcps)/cps

savfile = datadir+'ma_36'+numberstr+'_irr_at_1au.sav'
save, file=savfile, wave, irradiance, relEerr, /compress

stop
return
end
