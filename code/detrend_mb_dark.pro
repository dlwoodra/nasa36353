pro imgstats, in, imgmean, imgstd, imgmin, imgmax, imgmed
  
  imgstd  = stddev(in.image,dim=3)
  imgmean = mean(in.image,dim=3)
  imgmin  = min(in.image,dim=3, max=imgmax)
  imgmed = median(float(in.image),dim=3,/even)

  return
end

function best_imgmean, input_in

  input = input_in
  imgstats, input, themean, thestd, themin, themax, themed

  ; perform aggressive particle filtering
  ; based on flight EVE dark_image_particle_filter.pro 
  maxdifference = 7 ; DN
  ; the rocket is warmer/noisier than flight EVE
  ; increased maxdifference to reflect larger noise
  adjusted_input = dark_image_particle_filter((input.image), maxdifference)
  
  postmean = mean(adjusted_input,dim=3) ; short cut
  return,postmean

  ; reassign and continue filtering
  input.image = adjusted_input
  
  ;stop
  
  ; median filter input.image to remove spikes/splotches
  for i=0,n_elements(input)-1 do begin
     prev=i-1
     next=i+1
     input[i].image = median(input_in[[prev,i,next]].image,dim=3)
  endfor
  
  ; replace values outside reasonable range
  ; the normal std is around 4-4.5 DN
  for i=0,n_elements(input)-1 do begin
     ;bad=where(input[i].image gt (themin+5.),n_bad)
     bad=where(input[i].image gt (themin+3.*thestd),n_bad)
     if n_bad gt 0 then input[i].image[bad] = (themin[bad]+5.);<themed[bad]<themean[bad] ; aim for middle
  endfor
  imgstats, input, postmean, poststd, postmin, postmax
  img = postmean

  return, img
 end

function detrend_mb_dark, input, preidx, postidx

  ; preidx and postidx are proportional to time

  output=input

  corr = fltarr(2048,1024,n_elements(output))

  r1img = best_imgmean(output[preidx]) ; more noisy
  r2img = best_imgmean(output[postidx]) ; less noisy
  
  ;r1img = mean(output[preidx].image,dim=3) ; [2048,1024]
  ;r2img = mean(output[postidx].image,dim=3)
  r1ctr = mean(preidx)
  r2ctr = mean(postidx)
  slopeimg = (r2img-r1img)/(r2ctr-r1ctr) ; slope = dy/dx
  gd=where(abs(slopeimg) lt .15,comp=bad) ; upper limit on slope
  slopeimg[bad] = 0.
  offsetimg = r1img - slopeimg*r1ctr ; b = y1 - m*x1

  ;nx = n_elements(output[0].image[*,0])
  ;ny = n_elements(output[0].image[0,*])
  for i=0,n_elements(output)-1 do begin
     fit = slopeimg*float(i) + offsetimg
     bad=where(abs(output[i].image - fit) gt 10,n_bad)
     if n_bad gt 0 then fit[bad]=0. 
     corr[*,*,i] = slopeimg*float(i) + offsetimg ;line
  endfor

  for i=0,n_elements(output)-1 do output[i].image -= corr[*,*,i]
  ;stop
  return,output
end
