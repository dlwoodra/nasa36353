function convert_to_02nm, indat, max_ma_wave
  n_wave = (105L - 6 )/.02          ;- 1 ; 4950
  wave = dindgen(n_wave)*.02 + 6.00 ; match MEGS-B locations
  outdat = dblarr(3,n_wave)
  outdat[0,*] = wave
  aunc = reform(indat[2,*] * indat[1,*]) ; uncertainty with units
  tmp=min(abs(wave - max_ma_wave), max_ma_wave_idx)
  for i=0L,max_ma_wave_idx do begin
     outdat[1,i] = (indat[1,(i*2)-1] + indat[1,i*2]) * 0.5d ; mean
     outdat[2,i] = sqrt(aunc[(i*2)-1]^2 + aunc[(i*2)]^2) / 2.
  endfor
  outdat[2,*] /= outdat[1,*] ; change back to relative uncertainty
                                ; copy MEGS-B portion
  tmp=min(abs(indat[0,*] - max_ma_wave),in_idx)
  outdat[1,max_ma_wave_idx:*] = indat[1,in_idx:*]
  outdat[2,max_ma_wave_idx:*] = indat[2,in_idx:*]  
  return,outdat
end


function get_nrleuv, rwave  

  workingdir=file_dirname(routine_filepath()) ; in code/compare_merge
  datadir = workingdir+'/../../data/'

  ; for reference, Harry Warren showed an NRLEUV spectrum
  ; at 120 nm, Harry gets 1e8 ph/cm^2/sec
  nrlph = read_dat(datadir+'nrleuv_ref_qs_hr_v2.dat')
  nrlhrwave = reform(nrlph[0,*]) / 10. ; angstroms to nm
  deltaw = rwave-shift(rwave,1) & deltaw[0] = deltaw[1]
  ;binratio = 1./(deltaw)
  ;binratio *= (.02/(rwave-shift(rwave,1))) ; relative to MB sampling
  ;stop
  nrlhrirr_orig = ph2watt(nrlhrwave, reform(nrlph[1,*])) ; ph/cm^2/sec to W/m^2/bin
  nrlhrirr_orig /= .005 ; PHIL (div by bin size)
  ;expect binratio to be factor of 2 for rwave at .01 nm sampling
  ;since nrl is half angstrom sampling 

  ;FWHM = (rwave[1]-rwave[0])/(nrlhrwave[1]-nrlhrwave[0]) ; nm/(pixel size)
  FWHM = 0.01/.005 ; nm/(pixel size) ; PHIL
  sigma = FWHM/2.d/sqrt(2.d*alog(2)) ; PHIL
  x = findgen(801)/10.-40. ; PHIL
  kernel = exp(-x^2/2./sigma^2) ; PHIL
  nrlhrirr=convol(reform(nrlhrirr_orig), kernel, /normalize) ; PHIL

  ;; convolve with 1 angstrom triangle kernel to match EVE
  ;kernel = [1,2,3,4,5,6,7,8,9,10,11,10,9,8,7,6,5,4,3,2,1]
  ;kernel /= total(kernel)
  ;nrlhrirr = convol(nrlhrirr_orig, kernel, /edge_truncate, /center, /nan)
  
  data = {wave:rwave, irr:rwave*0.,nrlwave:nrlhrwave,nrlirr:nrlhrirr}
  delta = rwave[1] - rwave[0]
  low = rwave - 0.5d*(delta)
  hiw = low + delta
  for i=0L,n_elements(rwave)-1 do begin
     gd=where(nrlhrwave ge low[i] and nrlhrwave le hiw[i],n_gd)
     ; integrate nrl onto rocket wavelength scale
     if n_gd gt 0 then $
        if total(nrlhrirr[gd]) gt 0 then $
           data.irr[i] = mean( nrlhrirr[gd] )
        ;if total(nrlhrirr[gd]) gt 0 then $
        ;   data.irr[i] = int_tabulated(nrlhrwave[gd], nrlhrirr[gd],/double)
  endfor
;  data.irr *= binratio

  return,data
end

function test_get_nrleuv

  @config36353
  ;numberstr = '290'
  ;version='1_0'
  ;theyd = 2013294 ; oct 21

  restore,'rocket36'+numberstr+'_megsb_irr.sav' ; use spectra
  wave = spectra[0].w
  ;!p.color=0 & !p.background='ffffff'x
  ;plot,wave,spectra[25].sp,/ylog,yr=[1e-7,2e-3],xtit='Wavelength (nm)',ytit='Irradiance (W/m^2/nm)',xr=[33,106],xs=1
  nrl=get_nrleuv(wave)
  ;oplot,nrl.wave,nrl.irr,co='fe'x
  
  width=800 & height=800
  p=plot(wave,spectra[25].sp,xtitle='Wavelength (nm)',ytitle='Irradiance (W/m^2/nm)',$
    title='NASA Rocket 36.'+numberstr+' MEGS-B Spectrum with NRLEUV v'+version,$
         dimensions=[width,height], /stairstep, yrange=[1e-7,1e-2],xstyle=1,xrange=[33,106],/ylog)
  porig=p
  p=plot(nrl.wave,nrl.irr,/stairstep,co='blue',/overplot)
  p = text(/norm,.14,.9,'NRLEUV',co='blue')

  ; NEED TO CHANGE DATES/DOWNLOAD
  ;read_netcdf,'data/see__egs_L2_2018169_012_03.ncdf',see
  read_netcdf,'data/see__egs_L2_'+strtrim(theyd,2)+'_012_02.ncdf',see
  gd=where(see.wave gt 28 and see.wave lt 107)
  p=plot(see.wave[gd],see.flux_median[gd],/stairstep,co='red',/overplot)
  p = text(/norm,.14,.92,'TIMED-SEE',co='red')

  ;restore,'data/FISM_daily_2018169_v02_01.sav'
  restore,'data/FISM_daily_'+strtrim(theyd,2)+'_v02_01.sav'
  p=plot(fism_wv,fism_pred,co='purple',/stairstep,/overplot)
  p = text(/norm,.14,.94,'FISM',co='purple')


  ;restore,'data/megsb_irradiance_2018_old_code.sav'
  ;p = plot(wavelength, irradiance,co='dark green',/stairstep,/overplot)
  ;p = text(/norm,.14,.96,'OLD Method',co='dark_green')
  stop
  
  p.save,'rkt_nrlspectrum.png',width=width,height=height

  for wlo=32,102,5 do begin
     porig.xrange = [wlo,wlo+6]
     p.save,'rkt_nrlspectrum'+strtrim(wlo,2)+'.png',width=width,height=height
  endfor  
  stop
  p.close

  rkt1a={wave:wave, sp:spectra[25].sp}
  print,'saving rkt_comparison_data.sav'
  save,file='rkt_comparison_data.sav',fism_wv,fism_pred, see, nrl, rkt1a,/compress

  stop
  return,0
end

pro compare_with_old_rocket, outdat02, max_ma_wave

  @config36353

  workingdir=file_dirname(routine_filepath()) ; in code/compare_merge
  datadir = workingdir+'/../../data/'

  ; first restore the old rocket spectrum
  restore,datadir+'rachels_36258_Sep10.sav' ; wavelength, irradiance, uncertainty
  ; probably has different sampling for MA and MB
  gd=where(wavelength ge 6.0 and wavelength lt 105.0)
  olddat = dblarr(3,n_elements(wavelength[gd]))
  olddat[0,*] = wavelength[gd]
  olddat[1,*] = irradiance[gd]
  olddat[2,*] = uncertainty[gd]
  oldrkt36258 = convert_to_02nm(olddat, max_ma_wave) ; like level 2/2b
  ; there is something funny about oldrkt36258
  ; wavelengths longer than 27 nm are off by 1 bin
  ; wavelengths shorter than 27 nm are OK
  oldrkt36258[1,1051:-1] = oldrkt36258[1,1050:-2]
  
  ;call it oldrkt36258
  
  stop
  
  ;
  ; plot two rockets
  ;
  width=800 & height=800
  xr=[6,12]
  p1=plot(outdat02[0,*],outdat02[1,*],/ylog,yr=[1e-7,1e-3],xrange=xr,/stairstep,xstyle=1,$
          xtitle='Wavelength (nm)', ytitle='W/m^2/nm',title='Rocket Spectra',$
          LAYOUT=[1,2,1],dim=[width,height],margin=[.1,.1,.05,.1],co='black')
  op = plot(oldrkt36258[0,*],oldrkt36258[1,*],color='red',/stairstep,/overplot)

  t = text(/norm,.14,.88,'36.'+numberstr+' v'+version,co='black')
  t = text(/norm,.14,.86,'36.258 old',co='red')

  ratio = outdat02[1,*] / oldrkt36258[1,*]
  smwindow=9
  smratio = smooth(outdat02[1,*],smwindow,/edge_trunc) / smooth(oldrkt36258[1,*],smwindow,/edge_trunc)

  p2=plot(outdat02[0,*],ratio,xrange=xr,xstyle=1,/stairstep,/current,$
          xtitle='Wavelength (nm)',ytitle='Ratio New/Old',yr=[0.5,1.5],$
          LAYOUT=[1,2,2],margin=[.1,.1,.05,.1])
  for count=7,13 do op2=plot([0,200],[1,1]*count/10.,linestyle='dot',/overplot)
  op2=plot([0,200],[1,1],co='green',linestyle='dash',/overplot)
  ; overplot a 0.1 nm smooth
  opr=plot(outdat02[0,*],smratio,co='red',/overplot)
  p1.xrange = [5,105]
  p2.xrange = [5,105]
  p1.save,'mb_rkt_oldnew_comp.png',width=width,height=height
  stop

  ; plot the same thing with different wavelength ranges
  for i=5,100,5 do begin
     xr=[i,i+6]
     p1.xrange = xr
     p2.xrange = xr
     p1.save,'mb_rkt_oldnew_comp_'+string(i,form='(i03)')+'.png',width=width,height=height
  endfor
  stop
  p1.close
  
  return
end


pro compare_ma_mb_36353

@config36353
  
;;version='3_8'
;version='1_0'
;numberstr = '290'
;theyd = 2013294


workingdir=file_dirname(routine_filepath()) ; in code/compare_merge
datadir = workingdir+'/../../data/'

; read WHI solar ref
whi=read_dat(datadir+'ref_solar_irradiance_whi-2008_ver2.dat')
gd=where(whi[0,*] lt 200)       ; trim off all the longer wavelengths
whi=whi[*,gd]
whiwave=reform(whi[0,*])
whiirr =reform(whi[3,*])

;compare_ma_eve_36353, ma_wave, ma_irr, ma_err
restore,datadir+'ma_36'+numberstr+'_irr_at_1au.sav'
ma_wave = temporary(wave)
ma_irr = temporary(irradiance)
ma_err = temporary(relEerr)

;compare_mb_eve_36353, mb_wave, mb_irr, mb_err
restore,datadir+'mb_36'+numberstr+'_irr_at_1au.sav'
mb_wave = temporary(wave)
mb_irr = temporary(irradiance)
mb_err = temporary(relEerr)

;***
;***
;***
;; artificially bump up the uncertainties by 10% everywhere
;; until we get the 2017 calibration applied
;ma_err = (ma_err+0.1)<1.
;mb_err = (mb_err+0.1)<1.
;***
;***
;***

; fill ma spectrum below 6 nm
min_ma_wave = 6.0
max_ma_wave = 37. ; THIS CONTROLS MA/MB CROSSOVER
;max_ma_wave = 33.33 ; THIS CONTROLS MA/MB CROSSOVER
min_mb_wave=33.
max_mb_wave=105.0

fill_value = 1e-8

bad=where(ma_wave lt min_ma_wave)
ma_irr[bad] = fill_value
ma_err[bad] = 10.

; fill mb spectrum below 33 and above 105
bad=where(mb_wave lt min_mb_wave or mb_wave gt max_mb_wave)
mb_irr[bad] = fill_value
mb_err[bad] = 10.

; merge into one spectrum
maidx=where(ma_wave lt max_ma_wave)
mbidx=where(mb_wave ge max_ma_wave)

rwave = [ ma_wave[maidx], mb_wave[mbidx] ]
rirr = [ ma_irr[maidx], mb_irr[mbidx] ]
rerr = [ ma_err[maidx], mb_err[mbidx] ]

;plot,rwave,rirr,/ylog,yr=[1e-7,1e-2],xs=1
width=800 & height=800
p=plot(rwave,rirr,xtitle='Wavelength (nm)',ytitle='Irradiance (W/m^2/nm)',$
  title='NASA Rocket 36.'+numberstr+' MEGS-A, MEGS-B Spectrum v'+version,$
  dimensions=[width,height], /stairstep, yrange=[1e-7,1e-2],xstyle=1,xrange=[5,106],/ylog)
p.save,'rkt_fullspectrum.png',width=width,height=height
stop
p.close

outdat = transpose([[rwave],[rirr],[rerr]])
; trim off invalid ends
gd=where(outdat[0,*] gt min_ma_wave and outdat[0,*] lt max_mb_wave)
outdat=outdat[*,gd]

; version 1 - first cut, almost no flight corrections

comments = ['Date generated '+systime(),'',$
            '  Solar Spectral Irradiance from ',$
            '  NASA SDO EVE Rocket 36.'+numberstr+' Results',$
            '    MEGS-A and B Combined Spectrum',$
            '    '+humandatestr,$
            '------------------------------------',$
            '    Don Woodraska, Brian Templeman, Phil Chamberlin,',$
            '    Tom Woods, Frank Eparvier, Andrew Jones','',$
   ' This spectrum spans the EUV from approximately 6 to 105 nm.','', $
   '*****', $
   '***** WARNING: THIS IS PRELIMINARY DATA - DO NOT PUBLISH', $
   '*****', $
   '***** WARNING: This data is not ready for publication and may change.', $
   '*****', $
   '***** WARNING: Uncertainties are underestimated!', $
   '*****', $
   '', $
   ' Applied calibration from SURF 2017',$
   ' Applied NRL-MSIS00E atmospheric correction, about 3% at 30.4 nm',$
   '', $
   ' MEGS-A slit 1 wavelength range is 0.01 nm sampling from 6-14.93 nm', $
   ' MEGS-A slit 2 wavelength range is 0.01 nm sampling from 14.93-'+strtrim(max_ma_wave,2)+' nm', $
   ' MEGS-B wavelength range is 0.02 nm sampling from '+strtrim(max_ma_wave,2)+'-105 nm', $
   ' MEGS-B used 2017 380MeV scaled to 2010 ratio of 140/380MeV', $
   ' MEGS-A1 used 2017 285MeV adjusted for 3_7', $
   ' MEGS-A2 used 2017 140MeV', $
   '', $
   '*** pointing was not updated for this flight ***', $
   '', $
   'Version '+strtrim(version,2), $
   '', $
   ' Day of year='+strtrim(theyd mod 1000L,2), $
   ' Lyman-alpha composite v4=6.48e-3 W/m^2','',$
   ' F10='+strtrim(ft7,2)+', <F10>81day='+strtrim(ft7a,2), '', $
   ' Applied 1-AU adjustment '+strtrim(au,2),'',$
   ' Solar zenith angle (deg) '+strtrim(sza,2),'',$
   ' Apogee(km)='+apogeekmstr,' at T+'+apogeesecstr+' seconds', $
   'MEGS-A Slit 1 to slit 2 crossover occurs at 14.93nm', $
   'MEGS-A Slit 2 to MEGS-B crossover occurs at '+strtrim(max_ma_wave,2)+'nm', $
;      '', $
;      'The flight profile from 36.336 was used since it is not available for this flight.', $
;      '', $
;      'The apogee of 284 km occurred at T+273 seconds.',$
;   'Uncertainty is from SURF combined with in-flight scatter from 248-318 ', $
;   'seconds after launch','',$
   ''   ]

rktdatfile='rkt_36'+numberstr+'_irradiance_v'+strtrim(version,2)+'.dat'

write_dat, outdat, filename=rktdatfile,$
  comments=comments, $
  format='(f7.3,e13.6,f7.3)',$
  col='wavelength (nm), irradiance (W/m^2/nm), relative uncertainty (fraction)'
print,'wrote ',rktdatfile

; now write .02 nm sampling file
rktdatfile='rkt_36'+numberstr+'_irradiance_v'+strtrim(version,2)+'_02nm.dat'
outdat02 = convert_to_02nm(outdat,max_ma_wave) ; like level 2/2b
write_dat, outdat02, filename=rktdatfile,$
  comments=comments, $
  format='(f7.3,e13.6,f7.3)',$
  col='wavelength (nm), irradiance (W/m^2/nm), relative uncertainty (fraction)'
print,'wrote ',rktdatfile

rkt36336 = read_dat(datadir+'/rkt_36336_irradiance_v3_9_02nm.dat')

;
; plot two rockets
;
xr=[6,12]
p1=plot(outdat02[0,*],outdat02[1,*],/ylog,yr=[1e-7,1e-3],xrange=xr,/stairstep,xstyle=1,$
        xtitle='Wavelength (nm)', ytitle='W/m^2/nm',title='Rocket Spectra',$
        LAYOUT=[1,2,1],dim=[width,height],margin=[.1,.1,.05,.1],co='black')
op = plot(rkt36336[0,*],rkt36336[1,*],color='red',/stairstep,/overplot)

t = text(/norm,.14,.88,'36.'+numberstr+' v'+version,co='black')
t = text(/norm,.14,.86,'36.336 v39',co='red')

ratio = outdat02[1,*] / rkt36336[1,*]

p2=plot(outdat02[0,*],ratio,xrange=xr,xstyle=1,/stairstep,/current,$
        xtitle='Wavelength (nm)',ytitle='Ratio',yr=[0.5,1.5],$
        LAYOUT=[1,2,2],margin=[.1,.1,.05,.1])
for count=7,13 do op2=plot([0,200],[1,1]*count/10.,linestyle='dot',/overplot)
op2=plot([0,200],[1,1],co='green',linestyle='dash',/overplot)
   p1.xrange = [5,105]
   p2.xrange = [5,105]
p1.save,'mb_rkt_comp.png',width=width,height=height
stop

; plot the same thing with different wavelength ranges
for i=5,100,5 do begin
   xr=[i,i+6]
   p1.xrange = xr
   p2.xrange = xr
   p1.save,'mb_rkt_comp_'+string(i,form='(i03)')+'.png',width=width,height=height
endfor
stop
p1.close

; lather rinse repeat with the old rocket spectrum Rachel created Sep 2010
compare_with_old_rocket, outdat02, max_ma_wave

stop

!p.multi=[0,1,2]

;rkt1a=dblarr(n_elements(whi[0,*]))
;halfwidth=(whi[0,100]-whi[0,99])/2.d
;for i=0L,n_elements(rkt1a)-1 do begin
;   if whi[0,i] gt 5 and whi[0,i] lt 106 then begin
;      halfwidth=(whi[0,i]-whi[0,i-1])/2.d      
;      gd=where(rwave ge whi[0,i]-halfwidth and rwave lt whi[0,i]+halfwidth,n_gd)
;      if n_gd gt 0 then rkt1a[i] = mean(rirr[gd])
;  endif
;endfor
lobins=reform(whi[0,*] - (whi[0,1]-whi[0,0])/2.)
hibins=reform(whi[0,*] + (whi[0,1]-whi[0,0])/2.)
;rkt1aold=rkt1a
rkt1a = int_tabulated_array( rwave, rirr, lobins, hibins ) / (whi[0,1]-whi[0,0])
;plot,rwave,rirr,/ylog
;oplot,whi[0,*],rkt1aold,co='fe0000'x
;oplot,whi[0,*],rkt1a,co='fe'x
;stop,'compare rkt1aold and rkt1a'

nrl = get_nrleuv( rwave )
nrl1a = get_nrleuv( reform(whi[0,*]) )

xr=[6,12]
p1=plot(whi[0,*],whi[3,*],/ylog,yr=[1e-6,1e-3],xrange=xr,/stairstep,xstyle=1,$
        xtitle='Wavelength (nm)', ytitle='W/m^2/nm',title='Irradiance v'+version,$
        LAYOUT=[1,2,1],dim=[width,height],margin=[.1,.1,.05,.1],co='blue')
;op = plot(rwave,rirr,color='red',symbol='dot',/overplot)
op = plot(whi[0,*],rkt1a,color='red',/stairstep,/overplot)

;;op = plot(nrl.wave,nrl.irr,color='blue',symbol='dot',/overplot)
;op = plot(nrl1a.wave,nrl1a.irr,color='blue',/stairstep,/overplot)
;t = text(/norm,.14,.9,'NRLEUV',co='blue')
t = text(/norm,.14,.88,'36.'+numberstr+' v'+version,co='red')
t = text(/norm,.14,.86,'WHI',co='blue')
roverwhi = rkt1a/reform(whi[3,*])
;rovernrl = rkt1a/nrl1a.irr
p2=plot(whi[0,*],roverwhi,xrange=xr,xstyle=1,/stairstep,/current,$
        xtitle='Wavelength (nm)',ytitle='Ratio',yr=[0.5,1.5],$
        LAYOUT=[1,2,2],margin=[.1,.1,.05,.1],co='blue')
;op2=plot(whi[0,*],rovernrl,co='blue',/stairstep,/overplot)
for count=7,13 do op2=plot([0,200],[1,1]*count/10.,linestyle='dot',/overplot)
op2=plot([0,200],[1,1],co='green',linestyle='dash',/overplot)
;t = text(/norm,.14,.1,'36.'+numberstr+'/NRL',co='blue')
t = text(/norm,.14,.08,'36.'+numberstr+'/WHI',co='blue')
   p1.xrange = [5,105]
   p2.xrange = [5,105]
p1.save,'mb_comp.png',width=width,height=height

; plot the same thing with different wavelength ranges
for i=5,105,5 do begin
   xr=[i,i+6]
   p1.xrange = xr
   p2.xrange = xr
   p1.save,'mb_comp_'+string(i,form='(i03)')+'.png',width=width,height=height
endfor
stop
p1.close

for wlo = 6,106,5 do begin
  xr=[wlo,wlo+6]
  in = where(rwave gt xr[0] and rwave lt xr[1])
  yr = [min(rkt1a[in])>(1e-6),max(rkt1a[in])>(1e-3)]
  plot,whi[0,*],whi[3,*],/ylog,yr=[1e-6,1e-3],xr=xr,ps=10,xs=1,xtit='Wavelength (nm)',ytit='W/m^2/nm',tit='Irradiance v'+version
  oplot,rwave,rirr,co='fe'x,ps=3
  oplot,whi[0,*],rkt1a,co='fe'x,ps=10
;  oplot,nrl.wave,nrl.irr,co='fe0000'x,ps=3
;  oplot,nrl1a.wave,nrl1a.irr,co='fe0000'x,ps=10
  xyouts,xr[0]+0.1,5e-4,'WHI'
;  xyouts,xr[0]+0.5,5e-4,'NRLEUV',co='fe0000'x
  xyouts,xr[0]+0.1,2e-4,'36.'+numberstr,co='fe'x
  plot, whi[0,*], rkt1a / reform(whi[3,*]),xr=!x.crange,xs=1,ys=1,ps=10,tit='Ratio v'+version,ytit='Ratio',xtit='Wavelength (nm)',yr=[0.5,1.5]
;  oplot, whi[0,*], rkt1a / nrl1a.irr,co='fe0000'x,ps=10
  xyouts,xr[0]+.2,.6,'36.'+numberstr+' / WHI'
;  xyouts,xr[0]+.2,.7,'36.'+numberstr+' / NRLEUV',co='fe0000'x
  oplot,!x.crange,[1,1],lines=1
  stop
endfor

plot,whi[0,*],whi[3,*],/ylog,yr=[1e-7,1e-2],xr=[6,106],ps=10,xs=1,xtit='Wavelength (nm)',ytit='W/m^2/nm',tit='Irradiance v'+version
oplot,whi[0,*],rkt1a,co='fe'x,ps=10
xyouts,8.1,5e-3,'WHI'
xyouts,8,2e-3,'36.'+numberstr,co='fe'x
plot,whi[0,*],rkt1a/whi[3,*],yr=[0,2],ys=1,xs=1,ps=10,xr=!x.crange,ytit='36.'+numberstr+'/WHI v'+version,tit='Ratio',xtit='Wavelength (nm)'
oplot,!x.crange,[1,1],lines=1
stop

p=plot(rwave,rirr,xtitle='Wavelength (nm)',ytitle='Irradiance (W/m^2/nm)',$
  title='NASA Rocket 36.'+numberstr+' MEGS-A, MEGS-B Overlap v'+version,$
  dimensions=[width,height], /stairstep, yrange=[0,4e-4],xstyle=1,xrange=[33,39])

p=plot(nrl.wave,nrl.irr,/stairstep,co='blue',/overplot)
p = text(/norm,.14,.9,'NRLEUV',co='blue')

; should overplot the lines used for wavelength fitting here for validation
ma1linesfull=read_dat(datadir+'../code/wavelength/megsa1_solar_lines.dat')
ma1lines=ma1linesfull[where(ma1linesfull.cal ne 0)]
ma2linesfull=read_dat(datadir+'../code/wavelength/megsa_solar_lines.dat')
ma2lines=ma2linesfull[where(ma2linesfull.cal ne 0)]
; plot ma1lines
for i=0,n_elements(ma1lines)-1 do begin
  waveline = ma1lines[i].wave
  name = ma1lines[i].name
  junk=min(abs(rwave-waveline),thisidx)
  val = rirr[thisidx] ; rocket spectrum irradiance at theoretical line peak
  p=plot(/overplot,[1,1]*ma1lines[i].wave, [1e-8,val])
  t = text(/data, orientation=90, waveline, val, name+' '+strtrim(waveline,2))
end
; overplot ma2 lines
for i=0,n_elements(ma2lines)-1 do begin
  waveline = ma2lines[i].wave
  name = ma2lines[i].name
  junk=min(abs(rwave-waveline),thisidx)
  val = rirr[thisidx] ; rocket spectrum irradiance at theoretical line peak
  p=plot(/overplot,[1,1]*ma1lines[i].wave, [1e-8,val],color='red')
  t = text(/data, orientation=90, waveline, val, name+' '+strtrim(waveline,2),color='red')
end
stop

!p.multi=0


; L2b
d=eve_read_whole_fits(datadir+'/EVS_L2B_'+strtrim(theyd,2)+'_007_01.fit.gz')
; avg just 15 minutes
;evesp=mean(d.spectrum[1140:1155].irradiance,dim=2) ; filter?
; 6 sample smooth, then average 15 minutes
evesp=mean((smooth(d.spectrum.irradiance,[1,6],/edge_trun))[*,1140:1155],dim=2)
evewave=d.spectrummeta.wavelength



for i=33,36,3 do begin
  plot,ma_wave,ma_irr,xr=[i,i+3.5],xs=1,/ylog,yr=[1e-7,1e-3],ps=10,xtit='Wavelength (nm)',ytit='W/m^2/nm',tit='Version '+version ; MA
;  oplot,ma_wave,ma_irr*(1.+ma_err),ps=10
;  oplot,ma_wave,ma_irr*(1.-ma_err),ps=10

  oplot,mb_wave,mb_irr,co='fe0000'x,ps=10 ; MB
;  oplot,mb_wave,mb_irr*(1.+mb_err),co='fe0000'x,ps=10 ; MB
;  oplot,mb_wave,mb_irr*(1.-mb_err),co='fe0000'x,ps=10 ; MB

  oplot,evewave,evesp,co='fe'x,ps=10 ; EVE

  ;oplot,ma_wave*2.,ma_irr,lines=1

  xyouts,.12,.9,/norm,'EVE L2B 15-min avg V7',color='fe'x
  xyouts,.12,.85,/norm,'Rocket MEGS-A slit2 36.'+numberstr
  xyouts,.12,.8,/norm,'Rocket MEGS-B 36.'+numberstr,co='fe0000'x
stop
endfor

width=800 & height=600
p=plot(rwave,rirr,xtitle='Wavelength (nm)',ytitle='Irradiance (W/m^2/nm)',$
  title='NASA Rocket 36.'+numberstr+' MEGS-A, MEGS-B Overlap v'+version,$
  dimensions=[width,height], /stairstep, yrange=[1e-8,1e-3],xstyle=1,xrange=[33,39],/ylog)
p=plot(mb_wave,mb_irr,/stairstep,co='blue',/overplot)
p=plot(ma_wave,ma_irr,/stairstep,co='red',/overplot)
p = text(/norm,.14,.9,'MB',co='blue')
p = text(/norm,.2,.9,'MA',co='red')
stop

;convert to a function
integrate_ma_mb_overlap_lines, ma_wave, ma_irr, mb_wave, mb_irr, ctr_wave, ratioAoverB
;; low=[33.48d0,34.65,35.15,35.54, 36.70,37.33,37.94]
;; hiw=[33.60d0,34.85,35.35,35.65, 36.88,37.47,38.04]
;; ctr = (low + hiw)/2.d
;; trat=dblarr(n_elements(ctr))
;; for i=0L,n_elements(trat)-1 do begin
;;   gda=where(ma_wave ge low[i] and ma_wave le hiw[i])
;;   gdb=where(mb_wave ge low[i] and mb_wave le hiw[i])
;;   trat[i] = int_tabulated(ma_wave[gda],ma_irr[gda],/double) / $
;;             int_tabulated(mb_wave[gdb],mb_irr[gdb],/double)
;; endfor
plot,ctr_wave,ratioAoverB,ps=-4,xtit='Wavelength (nm)',ytit='MEGS-A/MEGS-B line ratio v'+version
stop

;p=plot(rwave,rirr,xtitle='Wavelength (nm)',ytitle='Irradiance (W/m^2/nm)',$
;  title='NASA Rocket 36.336 MEGS-A, MEGS-B Overlap v'+version,$
;  dimensions=[width,height], /stairstep, yrange=[0,4e-4],xstyle=1,xrange=[33,39])
;p=plot(mb_wave,mb_irr,/stairstep,co='blue',/overplot)
;p=plot(ma_wave,ma_irr,/stairstep,co='red',/overplot)
;p=plot(mb_wave,mb_irr*6.,/stairstep,co='purple',/overplot)
;p = text(/norm,.14,.9,'MB',co='blue')
;p = text(/norm,.2,.9,'MA',co='red')
;p = text(/norm,.07,.9,'MB*6',co='purple')
;stop

;arat=interpol(ma_irr,ma_wave,rwave) / rirr
;brat=interpol(mb_irr,mb_wave,rwave) / rirr
;
;plot,rwave,arat

stop
stop
return
end
