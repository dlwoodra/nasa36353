;+
; NAME:
;  megsb_tuned_mask
;
; PURPOSE:
;  Determine the best MEGS-B mask for the solar spectrum
;
; CATEGORY:
;  L1
;
; CALLING SEQUENCE:
;  IDL> megsb_tuned_mask, imgfull, imgmask, tunedmask
;
; INPUTS:
;  imgfull is one image [2048,1024]
;  imgmask is a guessed-mask that is larger than necessary
;
; OPTIONAL INPUTS:
;  none
;
; KEYWORD PARAMETERS:
;  none
;
; OUTPUTS:
;  tunedmask is the output mask
;
; OPTIONAL OUTPUTS:
;  none
;
; COMMON BLOCKS:
;  none
;
; SIDE EFFECTS:
;
;
;
; RESTRICTIONS:
;  Image is expected to be greater than 1e-14 (threshold)
;
;
; PROCEDURE:
;
;
;
; EXAMPLE:
;
;
;
; MODIFICATION HISTORY:
;  05/13/08 DLW modified to adjust line detection thresholds
;  
; $Log: megsb_tuned_mask.pro,v $
; Revision 1.4  2009/04/07 15:30:01  evesdp
; update
;
; Revision 1.3  2008/07/01 21:51:44  evesdp
; forcing strip to be 2 pixels narrower
;
; Revision 1.2  2008/05/15 22:51:49  evesdp
; tweaked thresholds
;
;-

pro megsb_tuned_mask, imgfull_in, imgmask, tunedmask

; INPUTS:
;  imgfull_in is one image [2048,1024]
;  imgmask is a guessed-mask that is larger than necessary
; OUTPUT:
;  tunedmask is the output mask

;rt=[763.023, -0.295935, -1.78369e-5]
;bt=[899.382, -0.287142, -1.56213e-5]
;; orig
;rt=[781.796, -0.276607, -2.49433e-05]
;bt=[928.212, -0.265820, -2.52693e-05]

  ; DLW hand examination of line edges from 36.353, 6/7/23
  ; wavelengths 33.5, 58.4, 102.5 define a quadratic
  ;stop

  ; these are the values from 36.389 for reference
  ;lowerx = [141,805,1946]
  ;lowery = [748,537,147]
  ;upperx = [132,802,1950]
  ;uppery = [891,694,316]

  ; 36.353 is very close to 36.389
  lowerx = [141,805,1945]
  lowery = [746,537,148]
  upperx = [132,803,1950]
  uppery = [891,695,316]
  rt = poly_fit(lowerx, lowery, 2)
  bt = poly_fit(upperx, uppery, 2)

  print,'INFO: megsb_tuned_mask - fit rt = ',rt
  print,'INFO: megsb_tuned_mask - fit bt = ',bt


  goto, rt_bt_known


;
; testing code, for tuning/fitting
;

;offset = 0.5
offset = 0.25
imgfull = imgfull_in*imgmask>offset<1
imgfull -= offset

poly_order=2


;tune the mask to match imgfull slit image height across the detector
; try 8 bins
;n_bins=64/2L ;number of bins
n_bins=64L ;number of bins
binsize=2048./n_bins ;256
xlo=(lindgen(n_bins))*(binsize) ;bin left-edges
xhi=xlo+binsize-1L ;bin right edges
topedge=lonarr(n_bins)
topedge--
botedge=topedge
xcm=topedge*1.
;mi=((imgfull>.2)-0.2)*imgmask
threshold=(1e-14)
;mi=((imgfull>threshold)-threshold)*imgmask
mi=median(((imgfull>threshold)-threshold)*imgmask,4)
;window,1
;!p.multi=0
xw=findgen(2048)
for i=0L,n_bins-1 do begin
    ;row sum
    sum=smooth(total(mi[xlo[i]:xhi[i],*],1),5,/edge_tr)
    ;wset,1
    ;plot,sum
    ;gd=where(sum gt 0.1,n_gd)
    gd=where(sum gt threshold*.5,n_gd)
    if n_gd gt 1 then begin
        medval=median(sum[gd],/even)
        ;half power on top
        j=0L
        gds=gd[n_gd/2]
        val=min(abs(sum[gd[0]:gds] - 0.5*medval),idx)
        topedge[i] = idx + gd[0]
        val=min(abs(sum[gds:gd[n_gd-1]] - 0.5*medval),idx)
        botedge[i] = idx + gds

        ;wset,0
        ;tvscl,hist_equal(congrid(mi[xlo[i]:xhi[i],*],1024,512))
        ;stop
    endif
endfor
bad=where(botedge-topedge lt 130,n_bad,comp=gd,ncomp=n_gd)
if n_bad gt 0 then begin
    botedge[bad]=-1
    topedge[bad]=-1
endif
xctr=(xlo+xhi)*0.5

for j=0L,1 do begin
    r=poly_fit(xctr[gd],(botedge[gd]-topedge[gd]),2,yfit=yfit,yband=yband)
    q=yfit-(botedge[gd]-topedge[gd])
    q *= q                      ;squared difference

    ;sort the differences from min to max
    qs=sort(q)
    ;stop
    gd=gd[qs[sort(xctr[qs[0:n_gd*.75]])]] ;keep lower 3 quartiles
    n_gd=n_elements(gd)
    ;refit only lower 3 quartiles
    r=poly_fit(xctr[gd],(botedge[gd]-topedge[gd]),2,yfit=yfit,yband=yband)
endfor

; fudge the top edge
topedge[gd[(n_elements(gd)-6):(n_elements(gd)-1)]] -= 7 ; lower by 5 pixels

rt=poly_fit(xctr[gd],topedge[gd],poly_order)
bt=poly_fit(xctr[gd],botedge[gd],poly_order)

print,'megsb_tuned_mask - rt = ',rt
print,'megsb_tuned_mask - bt = ',bt

;
; FLIGHT CODE
;
rt_bt_known:
tmpx=findgen(2048)
topcurve = poly(tmpx,rt)
botcurve = poly(tmpx,bt)
tunedmask=bytarr(n_elements(imgmask[*,0]),n_elements(imgmask[0,*]))

;
; 06/12/08 DLW narrow the tuned mask by one pixel on top and bot
;              (top/bot definition is upside-down, so bot>top always
;
;botcurve -= 1
;topcurve += 1

;
; 03/27/10 DLW Expand the tuned mask to ensure flux is not lost at long wavelengths
;
;botcurve += 1
;topcurve -= 1


; enforce valid ranges only!
;topcurve = (topcurve > 0) < botcurve
;botcurve = (botcurve > topcurve) < 1023
x=where(topcurve lt 0,n_x)
if n_x gt 0 then stop

for i=0L,2047 do begin
    tunedmask[i,topcurve[i]:botcurve[i]]=1
endfor

return
end
