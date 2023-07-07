;+
; This function removes the dark from all of the MEGS-B images in the
; sequence. It uses linear interpolation from a pre-solar dark and
; post-solar dark reference image pair. It was tuned for the 36.290
; rocket flight in Oct 21 (294), 2013 from WSMR.
; Modified for 36.258 then 36.353
; 
; :Params:
;    bmegs: in, required, type=array of structures
;      This is the input data structure
; :Returns:
;    Status is 0 if good, non-zero if bad
;-
function remove_megsb_dark_36353, bmegs, output, tunedmask

  @config36353
  
; determine a pre-observation dark and a post-observation dark
; goal is to linearly interpolate dark for each solar image


; create a dark region mask that surrounds the megsb solar strip
topmask=shift(tunedmask,0,180) & topmask[*,0:180]=0 ; this wraps vertically to the bottom
botmask=shift(tunedmask,0,-180) & botmask[*,1024-180:1023]=0
darkmask=topmask or botmask
;zero-out edge pixels
darkmask[0:3,*]=0 & darkmask[2044:2047,*]=0 ; left/right VC edges
darkmask[*,0] = 0 & darkmask[*,1023]=0 ; top/bottom rows
; use dark mask to locate pixels for each row to use for scaling

; define good image indices for dark images 
; good images do not show noise stripes across the width
; of the CCD

; 36.290 - MEGS-B
; 0-2 dark
; 3,4 corrupted dark
; 5 is flatfield
; 6,7,8 have low frequency weak wavy vertical noise pattern
; 9 dark
; 10,11 partial solar
; 12-18 solar
; 19-20 low frequency wavy vertical noise pattern
; 21-42 good solar
; 43,47 dimming
; 48,49 dark
; 50-53 ff
; 54-56 dark
;preidx=[20,21,22,23,24,29] ; not sure if 4 images is enough for a 3-image median
;postidx=[68,69,73,75,76]

preidx = dark1idx_b ;[81,83,84,86,89,90] ; remove_megsb_spikes
postidx = dark2idx_b ;[130,138,139]

; determine a pre-observation dark and a post-observation dark
; goal is to linearly interpolate dark for each solar image
predark = median(float(bmegs[preidx].image),dim=3,/even) ; added /even for 36.290
postdark= median(float(bmegs[postidx].image),dim=3,/even); added /even for 36.290

plot,total(tunedmask,2),ys=1,xs=1,xtit='column',ytit='mask height',tit='36.'+numberstr+' MEGS-B'
;stop
xs=1900
ys=1024
window,xs=xs,ys=ys
tvscl,hist_equal(congrid(tunedmask,xs,ys))
;stop
tvscl,hist_equal(congrid((bmegs[110].image - predark)>0,xs,ys))
;stop
tvscl,hist_equal(congrid(((bmegs[110].image - predark)>0)*tunedmask,xs,ys))
;stop
; these are candidates for eof/pca

wdelete

; prevent spikes in postdark
; plot postdark-predark and look for "normal" scatter
; the individual pixels should not be changing much
lobad=where(postdark-predark lt -6,n_lobad)
if n_lobad gt 0 then predark[lobad]=postdark[lobad]

; prevent spikes in predark
hibad=where(postdark-predark gt 5,n_hibad)
if n_hibad gt 0 then postdark[hibad]=predark[hibad]

; try a bandpass filter

;predark=postdark ; to test, use the same pre and post value ***
; define a reference index for the predark and postdark
pretime= preidx[-1] ;last one - might be overcorrecting 0 ; mean(preidx) ; move to beginning 
posttime = postidx[0] ;max(postidx) ;mean(postidx) 
stop
for i=0,n_elements(output)-1 do begin
  ; values before the pretime index
  if i lt pretime then output[i].dark = predark else begin
     ; values after the posttime index
     if i gt posttime then output[i].dark = postdark else begin
        ; values between pretime and posttime indices
        ; normal case is i > pretime and i < posttime
        ; interpolation is just a special combination of a weighted sum
        prefrac = (i-pretime) / (posttime-pretime)
        postfrac = 1. - prefrac
        output[i].dark = (predark*prefrac) + (postdark*postfrac)
     endelse
  endelse
endfor

for i=0,n_elements(output)-1 do begin
  output[i].dark = postdark
  ; adjust dark to match row by row
  ; only allow up to +40 DN above dark and -20 DN below dark
  deltaimg = (darkmask*(bmegs[i].image-output[i].dark))<20L > (-20L)
  for row=0L,1023 do begin
     ; use darkmask pixel differences, these are always dark
     ; except for flatfields
     if total(darkmask[*,row]) lt 19 then diffarr = median(deltaimg[*,row],19) else $
        diffarr = median(deltaimg[where(darkmask[*,row] eq 1),row],19)
     ; instead could we try a histogram-approach and select the peak?
     ; that would give the most likely, not the center value (median)
     ; probably needs more samples to determine accurately
     ; sticking with median
     gd=where(diffarr gt -20 and diffarr lt 20,n_gd)
     diff = mean(diffarr[gd]) ; exclude clipped portions
     ; skip flatfields
     if abs(diff) gt 19.9 then diff=0. ; just use normal dark
;if abs(diff) gt 1 then stop
     output[i].dark[*,row] += diff
  endfor
endfor
;stop

; 36.258 special top-hat adjustment
; this rocket was warmer than later rockets showing a significant
; square wave that changes differently in the middle than at the edges
; during the flight
; use the outer horizontal strips to estimate the additional
; adjustment needed
tmp = float(bmegs.image) - output.dark
;stop
;for i=0,n_elements(output)-1 do begin
;   ; TODO: skip flatfields
;   ;botstrip = mean(median(output[i].image[*,10:100],3), dim=2)
;   botstrip = mean(median(reform(tmp[*,10:100,i]),3), dim=2)
;   bimg=rebin(botstrip,2048,1024)
;   output[i].dark[*,0:511] += bimg
;
;   topstrip = mean(median(reform(tmp[*,922:1022,i]),3), dim=2)
;   timg = rebin(topstrip,2048,1024)
;   output[i].dark[*,512:1023] += timg
;endfor

; subtract the dark everywhere
output.image = float(bmegs.image) - output.dark
stop

; seems to be off by 0.25 DN before solar obs in bin 3213, 97.26 nm
; perform a fine-tune to enforce each column with a linear trend
newoutput = detrend_mb_dark(output, preidx, postidx)
; fine tuning was causing large negative values in output.image
; replaced with flight-filtering for L0D
output = newoutput
stop

; custom fix for 64.3 nm on 36.258
if numberstr eq '258' then begin
   ; detrend_mb_dark can overcorrect, so prevent large negative values
   for i=0,n_elements(output)-1 do begin
      bad=where(output[i].image lt -10,n_bad)
      output[i].image[bad] = 0.
   endfor
   ;stop
endif


output.total_counts = total(total(output.image,1,/double),1,/double)
for i=0L,n_elements(output)-1 do $
   output[i].sp_cps = total(output[i].image * tunedmask) / total(tunedmask)
output.sp_cps *= 0.1 ; 10 sec integrations
stop

return,0
end

;+
; Use the bandpass filter function
;-
function filter_rocket_image, image, lowpass, highpass
  newimg=float(image)
  if size(lowpass,/type) eq 0 then begin
    lowpass=0. ; needs to be zero to preserve magnitudes
    highpass=0.05 ; seems to work on solar images
  endif

  for col=0L,2047 do begin
    ; top first
    newimg[col,0:511]=bandpass_filter(image[col,0:511],lowpass,highpass)
    ; bottom
    newimg[col,512:1023]=bandpass_filter(image[col,512:1023],lowpass,highpass)
  endfor
return,newimg
end


;+
;-
pro mb_filterimgspike, arr, output_in, heavy=heavy

output = output_in[arr]       ; only use specified subset 9/8/20 for 36.258 2010 rocket

imgdim=size(output[0].image,/dim)

; loop over only the specified subset of images in output

enable_heavy_filter = 0 ; only heavily filter dark

; we cannot reliably determine when to apply heavy filtering, so the user
; has to specify it as a keyword
if keyword_set(heavy) then begin
   ; this is dark corrected dark data
   enable_heavy_filter=1
endif

; this part removes waves in the dark images and solar images

if enable_heavy_filter eq 0 then begin
   for i=0,n_elements(output)-1 do begin
      ; looks like each tap has a different noise threshold
      
      ; 0:511 (slit 2 with sam)
      next = (i+1) mod (n_elements(output)) ; wrap to compare last to first
      ; the value 10. is based on stddev of the signal spread
      ; dark has to already be removed from the image
      ; a rough estimate of dark is the VC column data
      noDarkSlit2 = output[i].image ; - median(output[i].image[2044:2047,0:511])
      nextNoDarkSlit2 = output[next].image ; - median(output[next].image[2044:2047,0:511])
      ; subtract approx dark from VC if VC mean is far above 0 DN
      ; TL is 0-3 that are shifted over to 2044-2047 in assembly
      if median(noDarkSlit2[2044:2047,0:511]) gt 500 then begin
         noDarkSlit2     -=     median(noDarkSlit2[2044:2047,0:511]) ; only approx dark
         nextNoDarkSlit2 -= median(nextNoDarkSlit2[2044:2047,0:511])
      endif
      bad=where( noDarkSlit2 gt (nextNoDarkSlit2 + sqrt(nextNoDarkSlit2>1.)*10.),n_bad)
      bot=output[i].image
      bot[bad]=output[next].image[bad] ; replace bad pixels using next image
      botspike=fltarr(imgdim)
      botspike[bad] = output[i].image[bad]

      ; 512:1023 (slit 1)
      noDarkSlit1 = output[i].image
      nextNoDarkSlit1 = output[next].image
      ; subtract approx dark from VC if VC mean is far above 0 DN
      ; BR is 2044-2047 that are shifted to 0-3 in assembly
      if median(noDarkSlit1[0:3,512:1023]) gt 500 then begin
         noDarkSlit1     -=     median(noDarkSlit1[0:3,512:1023])
         nextNoDarkSlit1 -= median(nextNoDarkSlit1[0:3,512:1023])
      endif
      ;bad=where(output[i].image gt output[next].image+sqrt(output[next].image>1.)*10.,n_bad)
      bad=where( noDarkSlit1 gt (nextNoDarkSlit1 + sqrt(nextNoDarkSlit1>1.)*10.),n_bad)
      top=output[i].image
      top[bad]=output[next].image[bad]
      topspike=fltarr(imgdim)
      topspike[bad] = output[i].image[bad]
   
      ; here combine the proper parts for each side
      output[i].image[*,0:511] = bot[*,0:511]
      output[i].spikes[*,0:511] = botspike[*,0:511]
      output[i].image[*,512:*] = top[*,512:*]
      output[i].spikes[*,512:*] = topspike[*,512:*]
   endfor
endif

; dark heavy filter added for 36.258 2010 rocket
; this cannot filter consecutive strikes at the same pixel very well
tmpout = output
if enable_heavy_filter ne 0 then begin
   for i=0L,n_elements(output)-1 do begin
      next = (i+1) mod (n_elements(output))           ; wrap to compare last to first
      prev = (i-1)                                     ; minus 1 is OK since IDL wraps automatically
      tmp = (output[next].image + output[prev].image)*0.5 ; mean
      ;bad=where(output[i].image gt tmp+5., n_bad) ; clip to only 5 DN above mean to allow some normal noise

      ; superior filtering uses a running 5-image median, temporal mean, and vertial 3 pixel median
      ; it might remove solar signal in solar images
      lo=(i-2) > 0
      hi=(lo+5) < (n_elements(output)-1L)
      lo = hi-5 ; re-adjust lo value if hi would be beyond the max
      ; can't do this next line with solar images
      ; HEAVILY RESTRICT ref value to within 5 DN of calculated estimates
      ref = median(output[lo:hi].image,dim=3) < (tmp+5.) < (median(output[i].image,3)+5.) ; combine spatial median
      ; HEAVILY RESTRICT image difference from ref to be less than 5 DN
      bad=where(abs(output[i].image - ref) gt 5., n_bad) ; clip to only 5 DN above mean to allow some normal noise
      if n_bad gt 0 then begin
         tmpout[i].spikes[bad] = output[i].image[bad] ; copy bad pixels to spikes image
         ;tmpout[i].image[bad] = tmp[bad]              ; replace bad pixels
         tmpout[i].image[bad] = ref[bad]               ; replace bad pixels
      endif
   endfor
   
   output = tmpout
endif

;; plot results
;window,1,xs=1024,ys=512
;window,2,xs=1024,ys=512
;for i=0,n_elements(output)-1 do begin
;   ; correct input index is arr[i] for input
;   wset,0
;   plot,yr=[-20,20], output_in[arr[i]].image,tit='filterimgspike #'+strtrim(arr[i],2)
;   oplot,output[i].image,co='fe'x
;   ;stop
;   ; for solar, look at images
;   wset,1
;   tvscl,hist_equal(congrid(output_in[arr[i]].image,1024,512))
;   ;stop
;   wset,2
;   tvscl,hist_equal(congrid(output[i].image,1024,512))
;   ;stop
;endfor
;wset,0

output_in[arr] = output ; copy back to input argument

return
end


;+
;-
function remove_megsb_spikes_36353, input, output

@config36353

imgdim=size(input[0].image,/dim)

output.image = input.image

; first cut is to remove jumps up pairwise, comparing to next image
; if current image is larger than next image by sqrt(DN) then replace with next
; only work on "good" solar data
;lo=11 ; 13 and 14 have RF noise
;hi=48 ; MA stops on 37 (34,35 show RF noise)
;identify types

; 36.290 - MEGS-B
; 0-2 dark
; 3,4 corrupted dark
; 5 is flatfield
; 6,7,8 have low frequency weak wavy vertical noise pattern
; 9 dark
; 10,11 partial solar
; 12-18 solar
; 19-20 low frequency wavy vertical noise pattern
; 21-42 good solar
; 43,47 dimming
; 48,49 dark
; 50-53 ff
; 54-56 dark

; 36.258
;dark1idx=[0,1,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,26,27,28,29]
;ffidx=[2,25,70] ; check
;;solaridx=[11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47]
;solaridx=[30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67]
;dark2idx=[68,69,71,72,73,74,75,76]

; 36.353
dark1idx = dark1idx_b ;[81,83,84,86,89,90] ; remove_megsb_spikes
ffidx = ffidx_b ;[133,134]
solaridx = solaridx_b ;[94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129]
dark2idx = dark2idx_b ;[130,138,139]


mb_filterimgspike, dark1idx, output, /heavy
mb_filterimgspike, ffidx, output ;, /heavy ; need at least 3 to run /heavy
mb_filterimgspike, solaridx, output
mb_filterimgspike, dark2idx, output, /heavy

return,output
end


pro make_mb_spectra, input, tunedmask, spectra, sensitivity_in=sensitivity_in, sens_sp=sens_sp, waveimg=waveimg, sensitivity_error_in=sensitivity_error_in, sens_err_sp=sens_err_sp

  @config36353
  
workingdir = file_dirname(routine_filepath()) ; in code
datapath = workingdir+'/../data/'

restore, datapath + 'rkt36'+numberstr+'_megsb_full_wave.sav'

wimg = waveimg

; create a mean image to get initial values
; choose max pixel from 97.7 at [1813,285] and a "dark" pixel as [1813,400]
tmp = reform(input.image[1813,285] - (input.image[1813,400]>1)) ; 97.7 signal at each time
sm = smooth(median(tmp,3),7,/edge_trunc) ; a time smoothed value of 97.7
tmp = max(sm,maxidx)

mb = mean(input[maxidx-5:maxidx+5].image,dim=3)
mb *= tunedmask

; make a reference spectrum
megsb_image_to_spectrum, mb, spflux, spwave, imgmask=tunedmask, imgwave=wimg

spectra_rec = {time:0.d, w:spwave, sp:spflux, err:spflux*0.}
spectra = replicate(spectra_rec, n_elements(input))
spectra.time=input.time

for i=0,n_elements(input)-1 do begin
   ;img=abs(input[i].image)>1e-9
   ; allow negative values, just not near zero values
   img=input[i].image
   bad=where(abs(img) lt 1e-12,n_bad)
   if n_bad gt 0 then img[bad]=1e-12
   megsb_image_to_spectrum, img, spwave, spflux, imgmask=tunedmask, imgwgd=tunedmask, imgwave=wimg, imgprec=sqrt(img), spprecout=sperr

   spectra[i].w = spwave
   spectra[i].sp = spflux
   ;spectra[i].err = sqrt(total(sperr^2,2,/double))
   spectra[i].err = sperr
   
   if max(spflux) gt 10 then $
     plot,spwave,spflux,/ylog,yr=[.1,1e6],tit='Image #'+strtrim(i,2) else $
     plot,spwave,spflux,/ylog,yr=[1e-8,1e-2],tit='Image #'+strtrim(i,2)
   wait,.01
   oplot,spwave,spflux          ; make it draw the screen
   heap_gc
   ;if i eq 20 then stop

   if i eq maxidx then begin
     p = image( hist_equal(congrid(input[i].image>0,1024,512)), $
      image_dimensions=[1024,512], dimension=[1024,512], rgb_table=1 )
     p.save,'megsb_'+strtrim(i,2)+'_img.png',width=1024,height=512
     p.close
   endif

endfor

if size(sensitivity_in,/type) ne 0 then megsb_image_to_spectrum,sensitivity_in*tunedmask,spwave,sens_sp,imgmask=tunedmask,imgwgd=tunedmask,imgwave=wimg

if size(sensitivity_error_in,/type) ne 0 then megsb_image_to_spectrum,sensitivity_error_in*tunedmask,spwave,sens_err_sp,imgmask=tunedmask,imgwgd=tunedmask,imgwave=wimg
;stop
return
end

pro clean_sens, sensitivity, sensitivity_error
  
  ; remove infinity, and nan from brian's sensitivity image and error image

  bad=where(finite(sensitivity) eq 0 or sensitivity ne sensitivity $
            ;or sensitivity lt 1e-8 or sensitivity gt 1e-3 $
            , n_bad)
  if n_bad gt 0 then begin
     stop
     sensitivity[bad]=1e-13
     ;  sensitivity_error[bad] = 1.
  endif

  ; new in v38
  ;if size(sensitivity_error,/type) eq 0 then sensitivity_error=sensitivity*0. + 1.0 ; 10a0%
  if size(sensitivity_error,/type) eq 0 then sensitivity_error=sensitivity*0.1 ; 10% in absolute units

  
; replace the unfilled parts with a linear extrapolation vertically
; the tunedmask will clip off parts where there is no sun
  colsum=total(sensitivity,2)   ; sum vertically

  min_valid_value = 1e-8
  
  for i=0,2047 do begin
     if colsum[i] lt 1e-6 and i lt 4 then sensitivity[i,*]=sensitivity[4,*] ; short edge
     if colsum[i] lt 1e-6 and i gt 2000 then sensitivity[i,*]=sensitivity[i-1,*] ; long edge
     ; fit a line to the "good" part
     gd=where( sensitivity[i,*] gt min_valid_value, n_good,comp=comp) ; lowest good value is around 4e-8
     if n_good gt 2 then begin
        fit = ladfit(gd,reform(sensitivity[i,gd]),/double)
        feval = fit[0] + fit[1]*dindgen(1024)
        feval = feval > min_valid_value ; min valid value to prevent values being too low
        org = reform(sensitivity[i,*])
        ; just use the middle part for stdev/med to avoid spiky edges
        range = gd[n_good*.25 : n_good*.75]
        std = stddev(org[range]-feval[range])
        med = median(org[range])
        thresh = std / med ; one-sigma
        plot, sensitivity[i,*], tit='col#'+strtrim(i,2)+' thresh='+strtrim(thresh,2)
        rel_diff = abs( reform(sensitivity[i,*]) - feval ) / feval
        ; replace top/bottom edges that show drop offs
        ; top edge
        idx = 1023
        while rel_diff[idx] gt thresh do begin
           sensitivity[i,idx] = feval[idx]
           idx--
        endwhile
        ; bottom edge
        idx = 0
        while rel_diff[idx] gt thresh do begin
           sensitivity[i,idx] = feval[idx]
           idx++
        endwhile
        oplot, sensitivity[i,*], co='fe'x
;        if i mod 100 eq 0 then stop
     endif else stop ; this should not happen, there should always be good data

     ; just replace sensitivity_error, there is nothing useful in it yet
     ; use difference from fit as a noise estimate
;     sensitivity_error[i,*] = (abs(org-feval)/(feval)) + 0.0001 ; at least .01%
;     sensitivity_error[i,*] *= sensitivity[i,*]*.2 ; give it units, limit to 20%
     
  endfor
  stop
  return
end

;+
; Need to tune for each rocket to try to correct bad images that are misassembled
;-
function fix_mb_corrupted_image, bmegs

; if special filtering is not known then just return
;return,bmegs

; replace corrupted images in 36.353 (old way uses a median)
; (old way) bmegs[95].image  = median(bmegs[[94,96]].image,dim=3,/even)
tmp = bmegs[95].image
; top and bottom rows are OK, but middle needs to be fixed
; have to fix top and bottom separately
rtmp = rotate(tmp,2) ; 180 deg rotation
; now need to shift to move VC to the right place
cn = 344
rtop = rtmp[*,0:511]
rtop = shift(rtop,cn) ; guess

rbot = rtmp[*,512:1023]
rbot = shift(rbot,-1*cn)
; examine then reassemble

new=tmp ; keep lowest 70 and highest 70 rows
rn = 70
new[*,rn:511] = rtop[*,rn:511]
new[*,512:1023-rn] = rbot[*,0:511-rn]
new[*,rn] = bmegs[94].image[*,rn] ; copy from previous image
new[*,1023-rn] = bmegs[94].image[*,1023-rn]

bmegs[95].image = new

; streak from [1896,148:511] between 100nm and 101nm in a dim place
; might be correctable with dark offset

return, bmegs
end


;+
; Custom rocket flight analysis for MEGS-B only
;
; :Examples:
;   IDL> s = read_mb_36336()
;
;-
function read_mb_36353

  @config36353

save_filtered_img = 0 ; set to 1 to write rkt36###_megsb_dark_particles_.sav
; only need to save once all the images are identified

window,0,xs=10,ys=10
wdelete
!p.color=0 & !p.background='ffffff'x & !p.charsize=1.5

; this expects a /code/ dir to contain this file and a /data/ dir at the same level as /code/ containing Tom's saesets

; change to this directory
workingdir = file_dirname(routine_filepath())
cd,workingdir

;tomsMASaveFile = workingdir+'/../data/TM2_36353_Flight_MEGS-A_adata.sav' ; contains adata
tomsMBSaveFile = workingdir+'/../data/TM2_36353_Flight_MEGS-B_bdata.sav' ; contains bdata

;BDATA           STRUCT    = -> <Anonymous> Array[140]
;IDL> help,bdata,/str
;** Structure <2639768>, 5 tags, length=4194328, data length=4194328, refs=1:
;   TIME            DOUBLE          -809.89897
;   PIXEL_ERROR     LONG                -2
;   FID_INDEX       LONG            147627
;   TOTAL_COUNTS    DOUBLE       2.8562857e+09
;   IMAGE           UINT      Array[2048, 1024]

if file_test(tomsMBSaveFile) ne 1 then begin
  print,'ERROR: cannot locate Toms MB save file '+tomsMBSaveFile
  stop
endif

print,'INFO: restoring Toms MEGS-B sav file'
restore,tomsMBSaveFile
print,'INFO: done'
; need to create bmegs array of structures to match 36.290 and 36.336
rec = {time:0.d, pixel_error:0L, fid_index:0, image:fltarr(2048,1024)}

tmp = replicate(rec,n_elements(bdata))
tmp.time = bdata.time
tmp.pixel_error = bdata.pixel_error
tmp.fid_index = bdata.fid_index
tmp.image = float(bdata.image)

; need data to be called bmegs, fix corrupted images, too
print,'INFO: calling fix_mb_corrupted_image'
bmegs = fix_mb_corrupted_image(tmp)

;stop


reltime = bmegs.time

bmegs.time = reform(reltime)
heap_gc

; look at each image, and the difference with the previous one
xsize=1920 & ysize=1024
window,0,xs=xsize,ys=ysize
;window,1,xs=800,ys=400
;window,2,xs=800,ys=400
for i=0,n_elements(bmegs)-1 do begin
   if bmegs[i].time lt 900 then continue
   wset,0
   tmpimg = float(bmegs[i].image)
   tmpimg[*,0:511] -= median(tmpimg[0:3,0:511])
   tmpimg[*,512:1023] -= median(tmpimg[2043:2047,512:1023])
   ;tvscl,hist_equal(congrid(bmegs[i].image,xsize,ysize))
   tvscl,hist_equal(congrid(tmpimg,xsize,ysize))
;   device,decomp=0
;   tvscl,hist_equal(congrid(bmegs[i].image mod 256, xsize, ysize))
   wset,1
   previdx = (i + n_elements(bmegs) - 1) mod n_elements(bmegs) 
   deltaimg = float(bmegs[i].image) - float(bmegs[previdx].image)
   !p.multi=0
   plot,deltaimg,title='index #'+strtrim(i,2)+'-#'+strtrim(previdx,2)+' T='+strtrim(bmegs[i].time,2),yr=median(deltaimg)+[-1.,1.]*6.*stddev(deltaimg),ytit='Difference (DN)'
   wset,2
   !p.multi=[0,1,2]
   plot,total(deltaimg,1)/2048.,ystyle=1,xstyle=1, $
      title='horizontal mean difference of #'+strtrim(i,2)+'-#'+strtrim(previdx,2),xtit='column',ytitle='mean DN/column'
   plot,total(deltaimg[*,0:511],2)/512.,ystyle=1,xstyle=1,$
      title='vertical mean difference of #'+strtrim(i,2)+'-#'+strtrim(previdx,2),xtit='row',ytitle='mean DN/row'
;   device,decomp=1
   oplot,total(deltaimg[*,512:1023],2)/512.,co='fe'x
   if i eq 95 then stop
   stop
end
;stop

;stop

; fix/replace corrupted images

bmegs.image[2044,0] = bmegs.image[2045,0] ; replace framestart bit in MEGS-B rocket
bmegs.image[3,1023] = bmegs.image[3,1022]     ; replace other half pair


; inventory the images
; 36.353 MEGS-B
; 
; ***
; LAUNCH, internal power change may change dark/noise
; ***
; 81 dark has diagonal ripple patterns
; 82 dark has diagonal and horizontal ripple patterns - discard
; 83,84 dark has diagonal ripple patterns
;
; 85 dark has diagonal and horizontal ripple patterns - corrupted
;
; 86 dark has diagonal ripple patterns
; 87 FF partially on
; 88 dark has diagonal and horizontal ripple patterns - discard
; 89,90 dark has diagonal ripple patterns 
; 91 dark shows both horizontal and diagonal ripple pattern - might be OK
; 92 first light smear - acquiring sun - discard
; 93 nearly full spectrum - acquiring sun, not centered, horizontal pattern - discard
; 94 sun centered full spectrum - horizontal pattern and diagonal pattern 
; 95 OK 
; 96 sun centered full spectrum - diagonal and horizontal pattern 
; 97 sun centered full spectrum - diagonal pattern
; 98 sun centered full spectrum - diagonal pattern 
; 99 sun centered full spectrum - diagonal and horizontal pattern (short side)
; 100 sun centered full spectrum - diagonal pattern 
; 101 sun centered full spectrum - diagonal and horizontal pattern
; 102,103 sun centered full spectrum - diagonal pattern 
; 104,105,106 sun centered full spectrum - diagonal pattern
; 107 sun centered full spectrum - diagonal and horizontal pattern
; 108,109 sun centered full spectrum - diagonal pattern
; 110 sun centered full spectrum - diagonal and horizontal pattern
; 111,112 sun centered full spectrum - diagonal pattern
; 113 sun centered full spectrum - diagonal and horizontal pattern
; 114 sun centered full spectrum - diagonal pattern 
; 115 sun centered full spectrum - diagonal and horizontal pattern
; 116,117 sun centered full spectrum - diagonal pattern
; 118 sun centered full spectrum - diagonal and horizontal pattern 
; 119 sun centered full spectrum - diagonal pattern
; 120 sun centered full spectrum - diagonal pattern 
; 121 sun centered full spectrum - diagonal and horizontal pattern
; 122,123 sun centered full spectrum - diagonal pattern
; 124 sun centered full spectrum - diagonal and horizontal pattern 
; 125,126 sun centered full spectrum - diagonal pattern
; 127 sun centered full spectrum - diagonal and horizontal pattern  noticable dimmer solar signal
; 128 sun centered full spectrum - diagonal pattern very dim
; 129 sun centered full spectrum - diagonal pattern  very dim nearly all solar signal is gone from short wavelengths
; 130 dark has diagonal ripple patterns 
; 131 dark shows both lower frequency horizontal and diagonal ripple pattern - discard
; 132,133,134,135 FF (132,135 is too dim)
; 136,137 shows both lower frequency horizontal and diagonal ripple pattern - discard
; 138,139 dark has diagonal ripple patterns 

;dark1idx_b=[81,83,84,86,89,90] ; remove_megsb_spikes
;ffidx_b=[133,134]
;solaridx_b=[94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129]
;;dark2idx=[130,136,137,138,139]
;dark2idx_b=[130,138,139]
;
;preidx=dark1idx_b ; for remove_megsb_dark
;postidx=dark2idx_b

; 36.258 MEGS-B
; 0 small waves in dark
; 1 corrupted dark partially reversed
; 2,25 flatfield
; 3 dark - does not match 5, sharp top-hat
; 4,6,13,15,17,19,28 larger waves, funny dark
; 5,7,8w,9,10,11w,12,14,16,18,20,21,22,23,24,29 dark dull top-hat
; (best 20,21,29)
; 26 post-flatfield bad dark
; 27 corrupted
; 30 first partial sunlight, small waves
; 31,32,33,34,35w,36,37w,38,39,40,41,42W,43,44W,45,46,47,48,49W,50,51W,52,53w,54,55,56W,57,58W,59,60w,61,62,63W,64dim,65mostly absorbed,66mostly absorbed,67only >977 sun
;  --dark lines near 50 nm and 82 nm
; 68,69 dark
; 70 flatfield
; 71 post-flatfield bad dark
; 72w,73,74w,75,76 dark



; 36.290 - MEGS-B
; 0-2 dark
; 3,4 corrupted dark
; 5 is flatfield
; 6,7,8 have low frequency weak wavy vertical noise pattern
; 9 dark
; 10,11 partial solar
; 12-18 solar
; 19-20 low frequency wavy vertical noise pattern
; 21-42 good solar
; 43,47 dimming
; 48,49 dark
; 50-53 ff
; 54-56 dark

orig=bmegs
; adjust locations by 4 pixels to fix real pixel locations
; shift columns to align lines across center (only one way will work)
; directions are opposite because readout is from opposite sides
for i=0,n_elements(bmegs)-1 do begin
   bmegs[i].image[*,0:511] = shift(bmegs[i].image[*,0:511],-4,0)
   bmegs[i].image[*,512:1023] = shift(bmegs[i].image[*,512:1023],4,0)
endfor

print,'INFO: figuring out MEGS_B tunedmask'
restore,'../data/megsb_default_mask.sav' ; larger than necessary mask
zmask = default_mask
print,'INFO: default mask'
tvscl,hist_equal(congrid(zmask * bmegs[110].image,xsize,ysize))
wait,1
;stop
; adjust for 36.290 to get close to "tuned"
topmask=shift(zmask,0,5)
botmask=shift(zmask,0,15)
zmask=topmask and botmask ; only use where both are good

; clip long wavelength edge
zmask[2043:*,*]=0
print,'INFO: adjsted default mask to prevent vertical clipping'
tvscl,hist_equal(congrid(zmask * bmegs[110].image,xsize,ysize))
wait,1
;stop

; apply mask to all images
; remove rough dark to get tuned mask
;tmpimg=median(bmegs[20:25].image,dim=3)
tmpimg=median(bmegs[100:110].image,dim=3)
tmpimg[*,0:511] = tmpimg[*,0:511] - median(tmpimg[*,0:511])
tmpimg[*,512:*] = tmpimg[*,512:*] - median(tmpimg[*,512:*])

megsb_tuned_mask, tmpimg>0, zmask, tunedmask

print,'INFO: tuned tunedmask to 36.353/36.389'
tvscl,hist_equal(congrid((tmpimg*tunedmask)>0,xsize,ysize))
;wait,1
stop

; for some reason, megsb_tuned_mask is off by ~30 pixels
;tunedmask = shift(tunedmask,0,30) ; shift it to match solar data
; it also needs to be taller vertically to prevent clipping bright lines

; need to grow the tuned mask for 36.258
;tunedmasklow=shift(tunedmask,0,27)
;tunedmaskhigh=shift(tunedmask,0,35)
;tunedmask = tunedmasklow OR tunedmaskhigh

; tunedmask is pretty good for 36.258, too

;; hand tune to 36.353
;tunedmasklow = shift(tunedmask,0,1)
;tunedmaskhigh = shift(tunedmask,0,-3)
;tunedmask = tunedmasklow OR tunedmaskhigh

;window,xs=1024,ys=512
;print,'INFO: no mask'
;tvscl,hist_equal(congrid(tmpimg>0,xsize,ysize)
;stop
;print,'INFO: tuned tunedmask to 36.353'
tvscl,hist_equal(congrid((tmpimg*tunedmask)>0,xsize,ysize))
wait,1
;stop
wdelete

; change data type of image to float
imgdim=size(bmegs[0].image,/dim)
newrec={time:0.d, pixel_error:0L, fid_index:0L, total_counts:0.d, sp_cps:0.d, $
  image:fltarr(imgdim), dark:fltarr(imgdim), spikes:fltarr(imgdim)}
result = replicate(newrec,n_elements(bmegs)) ; create output structure

; copy data that is not changing
result.time = bmegs.time
result.pixel_error = bmegs.pixel_error
result.fid_index = bmegs.fid_index

orig = result ; keep a copy for testing

print,'INFO: removing dark'
; remove dark
status=remove_megsb_dark_36353(bmegs, result, tunedmask)
nodarkimg=result

stop
; now do particle filtering
print,'INFO: particle filtering'
status=remove_megsb_spikes_36353(nodarkimg, result)
nospikes=result

; save these results for Phil and the others
if save_filtered_img eq 1 then begin
   print,'saving intermediate data for analysis by others'
   tmp=strsplit(systime(),/extract)
   ts=tmp[1]+'_'+tmp[2]+'_'+tmp[4]
   file='../data/rkt36353_megsb_dark_particles_'+ts+'.sav'
   rec={time:0., image:fltarr(2048,1024)}
   gd=where(nospikes.time gt -60,n_gd)
   mb_no_spikes = replicate(rec,n_gd)
   mb_no_spikes.time = nospikes[gd].time
   mb_no_spikes.image = nospikes[gd].image
   stop,'*** you really do not want to save this if it is not necessary***'
   save,file=file, mb_no_spikes, /compress
   mb_no_spikes=0b
   print,'INFO: read_mb_36353 - saved '+file
endif
   
heap_gc
stop

; this is where we save data for Phil and other to image analysis
; on particle filtered, dark corrected data


;; ***
;; TESTING CODE
;; TRY REMOVING PARTICLES, THEN REMOVING DARK
;status=remove_megsb_spikes_36353(bmegs, result)
;nospikes2=result
;status=remove_megsb_dark_36353(nospikes2, result, tunedmask)
;nodarkimg2=result
;; compare new nodarkimg2 and old nospikes
;stop
;; Results - changing the order increases the mean and median total
;;           irradiance in the spectrum, but allows more particles
;;           through, seems to make lines brighter - need to check this
;; ***


;from Michael
;restore,'data/SURF_Sep10_sensitivity_megs_b_fw0.sav' ; sensitivity includes gain
; from Brian 3/4/19
;restore,'data/sept_2017_sensitivity_MEGSB_fw0_380MeV.sav' ; sensitivity includes gain

;; from Brian 7/19/19 v38 ***
;sensfile='data/sensitivity_MEGSB_fw0_380MeV_2017.sav' ; scaled to A2 in overlap
; alpha, beta, filter, instrument, sensitivity[2048,1024], sensitivity_error[2048,1024]

; for 36.290 (not used)
;sensfile='data/sensitivity_MEGSB_fw0_380MeV_2013.sav' ; scaled to A2 in overlap?
print,'INFO: read_mb_36353 - calling get_36353_sensitivity'
sensitivity = get_36353_sensitivity(/megsb)

; from Brian 9/13/19 v39 rejected
;sensfile='data/sensitivity_MEGSB_second_order_2017.sav' ; 
; sensitivity[2048,1024], sensitivity_vec[3700], wave_vec[3700]
;    this file contains the 2017 MEGS B sensitivity map that is 
;    corrected for second order using the ratio of 380MeV to 140MeV from 2009 
;    order sorting data.  That ratio is then applied to the 2017 380MeV 
;    sensitivity map.
; THE Helium CONTINUUM RATIO TO WHI IS ~0.5.

; from Brian 9/15/19 v310 candidate rejected
;sensfile='data/sensitivity_MEGSB_second_order_full_2017.sav' ; 
; sensitivity[2048,1024], sensitivity_vec[3700], wave_vec[3700]
;   Brian generated a new sensitivity map applying the second order
;   correction to the entire image instead of just to wavelengths
;   above 72nm.
; THIS MADE THE Helium continuum ratio to WHI 0.4

; from Brian 9/17/19 v39 candidate
;sensfile='data/sensitivity_MEGSB_second_order_v39_2017.sav' ; 
; sensitivity[2048,1024], sensitivity_vec[3700], wave_vec[3700]
;    This version does not scale the various sensitivity maps up or
;    down to match irradiance at any particular wavelength.
; THE Helium CONTINUUM RATIO TO WHI IS ~0.5.

; from Brian 9/17/19 v39 candidate
;sensfile='data/sensitivity_MEGSB_second_order_v40_2017.sav' ; 
; sensitivity[2048,1024], sensitivity_vec[3700], wave_vec[3700]
;    This version does not scale the various sensitivity maps up or
;    down to match irradiance at any particular wavelength.
; THE Helium CONTINUUM RATIO TO WHI IS ~0.5.

; from Brian 9/27/19 v39 candidate
;sensfile='data/sensitivity_MEGSB_second_order_v39_2017.sav' ; 
; sensitivity[2048,1024]
;    This version has the MEGS-A adjustment to wavelengths
;    shorter than 42.399nm.
; THE Helium CONTINUUM RATIO TO WHI IS ~0.5.
; The Hydrogen continuum is abou 50% higher than WHI, higher than
; solar max. Rejected, trying again with same name
; seems slightly better, but not enough short wavlength changes
; try to manually fix



; from Brian 4/11/19
;sensfile='data/sensitivity_MEGSB_fw0_140MeV_2017_solar.sav' ; scaled to A2 in overlap
;sensfile='data/sensitivity_MEGSB_fw0_140MeV_2017.sav'
;sensfile='data/sensitivity_MEGSB_fw0_120MeV_2017.sav'

;restore,sensfile                ; sensitivity includes gain
; need to trim Brians sensitivity map to remove bad vertical edges
sensitivity_orig=sensitivity


; V3.9 change
print,'INFO: read_mb_36353 - restoring wavelength map'
workingdir = file_dirname(routine_filepath()) ; in code
datapath = workingdir+'/../data/'
restore,datapath+'rkt36'+numberstr+'_megsb_full_wave.sav' ; waveimg
 ; make this a function that applies correction to sensitivity
 ; to make megsb match megsa in overlap
; mbCor = EXP( 10.950772 + waveimg *  (-0.25827815 ))
; idx = WHERE(mbCor LT 1.0 OR waveimg LT 33.0)
; mbCor[idx] = 1.0
; idx = WHERE( waveimg GT 42.399143)
; mbCor[idx] = 1.0
; sensitivity *= mbCor

print,'INFO: read_mb_36353 - calling clean_sens'
 clean_sens, sensitivity, sensitivity_error
 sensitivity_error = sensitivity_error < sensitivity ; clip to 100% error
 sensitivity_error = sensitivity^2*.01/min(total(sensitivity*tunedmask,2)) ; min is 1%
; sensitivity_error = total(sensitivity*tunedmask,2)/total(tunedmask,2) ; 1 DN

 ;stop

print,'INFO: read_mb_36353 - calling make_mb_spectra using nospikes'
 nospikes.image *= 0.1d         ; divide by integration time 10 seconds
make_mb_spectra, nospikes, tunedmask, spectra_cps, $
  sensitivity_in=sensitivity, sens_sp=sens_sp, waveimg=waveimg, $
  sensitivity_error_in=sensitivity_error, sens_err_sp=sens_err_sp

; prevent bad long wavelengths
bad=where(spectra_cps[0].w gt 106.1)
for i=0,n_elements(spectra_cps)-1 do begin
   spectra_cps[i].sp[bad] = spectra_cps[i].sp[bad[0]-1] ; pin to last value
   spectra_cps[i].err[bad] = 1.
endfor

print,'INFO: read_mb_36353 - saving rkt36'+numberstr+'_mb_countspectrum.sav'
save,file='rkt36'+numberstr+'_mb_countspectrum.sav',/compress,spectra_cps

; the image process has drop outs, but this is how to do it
calimg=nospikes
for i=0,n_elements(calimg)-1 do calimg[i].image = calimg[i].image*sensitivity
make_mb_spectra, calimg, tunedmask, spectra_cal, waveimg=waveimg
; prevent bad long wavelengths
bad=where(spectra_cal[0].w gt 106.1)
for i=0,n_elements(spectra_cal)-1 do begin
   spectra_cal[i].sp[bad] = spectra_cal[i].sp[bad[0]-1] ; pin to last value
   spectra_cal[i].err[bad] = 1.
endfor

stop

; alternative approach is to interpolate the averaged sensitivity as a
; spectrum
sens2048=total(sensitivity*tunedmask,2)/total(tunedmask,2) ;2048 ; 1 DN
;senserr2048=total(sensitivity_error*tunedmask,2)/total(tunedmask,2) ;2048
senserr2048 = (total(sensitivity_error*tunedmask,2))
; clip to at least 1%, and less than 100%
;senserr2048 = senserr2048 > (sens2048*.01) < sens2048

; NOT NEEDED in V38
;; remove zeros
;;stop
;cnt=0L
;while sens2048[cnt] lt 1e-12 do cnt++
;;cnt=10L
;sens2048[0:cnt-1]=sens2048[cnt]
;senserr2048[0:cnt-1]=senserr2048[cnt]
;cnt=2047L
;while sens2048[cnt] lt 1e-12 do cnt--
;sens2048[cnt:*]=sens2048[cnt]
;senserr2048[cnt:*]=senserr2048[cnt]



w2048=total(double(waveimg)*tunedmask,2,/double)/total(tunedmask,2,/double) ; 2048
sens3700=interpol(sens2048,w2048,spectra_cps[0].w)
senserr3700=interpol(senserr2048,w2048,spectra_cps[0].w)

;spectra=spectra_cps
;for i=0,n_elements(spectra)-1 do spectra[i].sp = (spectra[i].sp>.1)*sens3700
;simplespectra=spectra
; DUE TO STRANGE EFFECTS, REPLACE SPECTRA with SPECTRA_CAL 10/2/19
spectra=spectra_cal

ctr0 = solaridx_b[n_elements(solaridx_b)/2L]
ctr = ctr0 + [-5,5]
msp = mean(spectra_cps[ctr[0]:ctr[1]].sp,dim=2) ; mean "best" spectrum
xr=[solaridx_b[0],solaridx_b[-1]]

; atmos absorption
plot,ps=-4,spectra.sp[3213]/mean(spectra[ctr[0]:ctr[1]].sp[3213]),yr=[-0.05,1.2], xr=xr, tit='36.'+numberstr+' Wavelength '+strtrim(string(spectra[0].w[3213],form='(f6.3)'),2)+' nm',ys=1,xtit='Index'
oplot,ps=-1,spectra_cps.sp[3213]/mean(spectra_cps[ctr[0]:ctr[1]].sp[3213]),co='fe'x
oplot,!x.crange,[0,0],lines=1
oplot,!x.crange,[1,1]*0.9,lines=1
oplot,!x.crange,[1,1],lines=1
oplot,[1,1]*ctr[0],!y.crange,lines=1
oplot,[1,1]*ctr[1],!y.crange,lines=1
axis,90,.5,/data,xrange=spectra[round(!x.crange)].time,xs=1,xtit='Time (sec)'
stop

!p.multi=[0,1,2]

plot,spectra_cps[ctr0].w,spectra_cps[ctr0].sp,yr=[.1,1e4],/ylog, tit='MEGS-B 36.'+numberstr,ytit='cps',xtit='Wavelength (nm)',xs=1,ps=10
plot,spectra[ctr0].w,spectra[ctr0].sp,yr=[5e-7,1e-2],/ylog,tit='36.'+numberstr,ytit='cps*S (W/m^2/nm)',xtit='Wavelength (nm)',xs=1,ps=10

save,file='rocket36'+numberstr+'_megsb_irr.sav',/compress,spectra,spectra_cps,sens2048,sens3700,senserr3700
print,'INFO: read_mb_36353 - saved rocket36'+numberstr+'_megsb_irr.sav'

stop
return,1
end
