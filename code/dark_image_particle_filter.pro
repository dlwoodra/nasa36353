; docformat = 'rst'

;+
; :Author:
;    Don Woodraska
;
; :Copyright:
;    2023 Regents of the University of Colorado, Sponsored by NASA
;    under contract NAS5-02140. BSD3
;-

;+
; This is an internally used function for dark_image_particle_filter
; that performs the image filtering. It relies on consecutive image
; persistence to identify and replace values in the image array that
; exceed a difference threshold.
;
; :Params:
;    imgarr: in, required, type=fltarr
;      The unchanged image array used as a reference.
;      Dimensions are identical to darkarr.
;   maxdifference: in, optional, type=float
;      The maximum difference that is allowed from image to image in
;      time. This is only used to identify spikes up, not spikes down.
;      The default value is 5.
;      Values can only be positive. Code returns error if negative.
;      Users should carefully consider the impact to the distribution
;      of making this too small. A value of 2-3 sigma that works well
;      for EVE CCDs might not be large enough for some noisier data,
;      or might be too large for less noisy data.
;
; :Returns:
;   This function returns a new copy of imgarr with the values
;   clipped.
;
; :Example:
;   As an example, we create a 6 image array.
;   Image #1-4 have pixel 1,1 initially set to 200. The
;   dark_image_particle_filter code replaces image 1:4 pixel 1,1 with
;   the mean.
;   A rough test could be something simple as in ::
;      IDL> imgarr=findgen([6,2,2]) + 100 & imgarr[1:4,1,1]=200
;      IDL> print,imgarr
;      100.000      101.000      102.000      103.000      104.000      105.000
;      106.000      107.000      108.000      109.000      110.000      111.000
;      112.000      113.000      114.000      115.000      116.000      117.000
;      118.000      200.000      200.000      200.000      200.000      123.000
;      IDL> print,dark_image_particle_filter(imgarr,mean=mn)
;      100.000      101.000      102.000      103.000      104.000      105.000
;      106.000      107.000      108.000      109.000      110.000      111.000
;      112.000      113.000      114.000      115.000      116.000      117.000
;      118.000      120.500      120.500      120.500      120.500      123.000
;      IDL> print,mn
;      102.500      108.500
;      114.500      120.500
;
;   In this case, the mean of pixel 1,1 is (118+123)/2 = 120.5 since all
;   the rest were replaced.
;
;-
function dipf_clip_images, imgarr, maxdifference, use60=use60

   darkarr = float(imgarr)

   if maxdifference le 0. then begin
      print,'ERROR: dark_image_particle_filter-dipf_clip_images - maxdifference is not greater than zero, '+strtrim(maxdifference,2)
      stop ; fatal
   endif
   
   ; TEMPORAL PARTICLE FILTERING ONLY
   ; This should work well for dark data where values change little from one
   ; image to the next image.

   ; filter particle strikes larger then maxdifference from previous/next images
   ; only do 3 times (can remove up to 5 consecutive hits)

   nimg = n_elements(imgarr[*,0,0]) ; number of images
   if keyword_set(use60) then begin
      for unusedloopcounter=0, ((nimg-1) < 2) do $
         darkarr = darkarr $
         < shift(darkarr+maxdifference,1,0,0) < shift(darkarr+maxdifference,-1,0,0) $
         < shift(darkarr+maxdifference,2,0,0) < shift(darkarr+maxdifference,-2,0,0) $
         < shift(darkarr+maxdifference,nimg/2,0,0) < shift(darkarr+maxdifference,-1*nimg/2,0,0)
   endif else begin
      for unusedloopcounter=0, ((nimg-1) < 2) do $
         darkarr = darkarr $
         < shift(darkarr+maxdifference,1,0,0) < shift(darkarr+maxdifference,-1,0,0) $
         < shift(darkarr+maxdifference,2,0,0) < shift(darkarr+maxdifference,-2,0,0)
   endelse


   ; if memory allows, this could be improved by adding additional terms
   ; for +/-2, +/-3, etc.

   ; recall 1d to 2d indexing
   ; x2didx = 1didx mod n_elements(x)
   ; y2didx = 1didx / n_elements(x)
   ; so 1d pixel is at 2d position [x2didx, y2didx]
   
   ; repeat to remove consecutive spikes at the same pixel in neighboring images
   ; each iteration cleans more up to 3 iterations
   ; after 3 no changes occur
   ; To test how it creates a triangle, try this test case
   ;IDL> d=fltarr(10)+100. & d[2:6]=200 
   ; set 5 values in a row to something large
   ;IDL> print,d
   ;      100.000      100.000      200.000      200.000      200.000      200.000
   ;      200.000      100.000      100.000      100.000
   ; run first iteration and show new values
   ;IDL> d = (d<shift(d+5,1))<shift(d+5,-1) & print,d
   ;      100.000      100.000      105.000      200.000      200.000      200.000
   ;      105.000      100.000      100.000      100.000
   ; run second iteration and show new values
   ;IDL> d = (d<shift(d+5,1))<shift(d+5,-1) & print,d
   ;      100.000      100.000      105.000      110.000      200.000      110.000
   ;      105.000      100.000      100.000      100.000
   ; run third iteration and show new values
   ;IDL> d = (d<shift(d+5,1))<shift(d+5,-1) & print,d
   ;      100.000      100.000      105.000      110.000      115.000      110.000
   ;      105.000      100.000      100.000      100.000
   ; no change in 4th iteration
   ;IDL> d = (d<shift(d+5,1))<shift(d+5,-1) & print,d
   ;      100.000      100.000      105.000      110.000      115.000      110.000
   ;      105.000      100.000      100.000      100.000     
   ; a square wave becomes a triangle
   ; it does not matter what the shape is as long as it is different from the original
   return, darkarr
end


;+
; This is an internal function for use only by
; dark_image_particle_filter. It replaced pixels that exceed the
; difference with a NaN value.
;
; :Params:
;    darkarr: in, required, type=fltarr
;      The image array that will have NaNs.
;      Dimensions are [N,X,Y]
;    imgarr: in, required, type=fltarr
;      The unchanged image array used as a reference.
;      Dimensions are identical to darkarr.
;    total_replaced: out, optional, type=long
;      The number of pixels that are foudn to be different between
;      darkarr and imgarr that were replaced with NaN
;
; :Returns:
;    This function returns the 1-d array of indices into darkarr that
;    were changed.
;
;-
function dipf_replace_with_nan, darkarr, imgarr, total_replaced

   all_pixels_to_replace = where(abs(darkarr - float(imgarr)) gt .1, total_replaced) ; these are the pixels that were replaced

   darkarr[all_pixels_to_replace] = !values.f_nan
   return, all_pixels_to_replace
end


;+
; This function is only intended for images that are not supposed to
; be changing, like the dark images. This will remove solar photons or
; particle spikes/streaks from images that exceed the maxdifference.
;
; This function determines which pixels are different from the
; previous and next images by some amount. Those pixels are replaced
; by the mean of the remaining pixels and returned by this function.
;
; This function applies a temporal-only filter. No spatial filtering
; is performed.
;
; Optional keywords allow some internal calculation information to be
; returned like the total pixels replaced, number of pixels in each
; image that were replaced, and the mean of the filtered images.
;
; :Params:
;    imgarr: in, required, type=fltarr
;      The 3-d array containing N images of dimension X, and Y.
;      The dimension order is imgarr[N,X,Y]. Filtering is only 
;      applied over the first dimension. The N dimension corresponds
;      to time, and this expects consecutive images in time to have 
;      good values that are within maxdifference of each other. It is
;      recommended that N be at least 3, but memory limitations
;      dictate the maximum N. In practice, use 100 or less.
;   maxdifference: in, optional, type=float
;      The maximum difference that is allowed from image to image in
;      time. This is only used to identify spikes up, not spikes down.
;      The default value is 5.
;      Values can only be positive. Code returns error if negative.
;      Users should carefully consider the impact to the distribution
;      of making this too small. A value of 2-3 sigma that works well
;      for EVE CCDs might not be large enough for some noisier data,
;      or might be too large for less noisy data.
;
; :Keywords:
;    total_replaced: out, optional, type=long
;      The total number of pixels that were identified/replaced.
;    number_replaced_per_image: out, optional, type=lonarr
;      The number of pixels in each image that was replaced.
;      The array is the same lendgth as the first dimention of imgarr.
;    mean_of_images: out, optional, type=fltarr
;      The mean of imgarr calculated using values that are within
;      the maxdifference.
;    use60: in, optional, type=boolean
;      Set this keyword to perform additional heavy 60-second filtering
;
; :Returns:
;    This function returns an array of images identical to the argument
;    with values at or above maxdifference replaced with the mean.
;
; Modification History:
;  6/9/23 DLW Original file creation for use in EVE version 8.
; 
;-
function dark_image_particle_filter, imgarr, maxdifference, total_replaced=total_replaced, $
                                     number_replaced_per_image=number_replaced_per_image, $
                                     mean_of_images=mean_of_images, use60=use60

   if size(maxdifference,/type) eq 0 then maxdifference=5.

   if maxdifference le 0. then begin
      print,'ERROR: dark_image_particle_filter - cannot use maxdifference <= 0 ',strtrim(maxdifference,2)
      return,-1
   endif
   
   ; memory discussion
   ; EVE images are 2048x1024 each, as floats each is 8 MB
   ; a 100 image array is then 800 MB and there are 2 copies in memory
   ; more memory is needed for all the intermediate calculation as well

   ; WARNING
   ; This code can exceed physical memory on older machines if large
   ; numbers of images are passed in.
   
   n_images = n_elements(imgarr[*,0,0])

   number_replaced_per_image = lonarr(n_images)

   ; create darkarr replacing exceedences with clipped values
   darkarr = dipf_clip_images(imgarr, maxdifference)   
   
   ; get the array of indices and replace with nan
   all_pixels_to_replace = dipf_replace_with_nan( darkarr, imgarr, total_replaced )

    ;calculate the mean image using data that was not replaced
   mean_of_images = mean(darkarr, dim=1, /nan) ; mean of the image excluding filled values

   ; it might be possible to use fancy indexing and rebinning instead
   ; looping over each image uses less memory
   
   ; replace flagged pixels with the mean
   ; only needed for calculating the median
   for i=0, n_images-1 do begin
      pixels_to_replace = where( ~ finite(reform(darkarr[i,*,*])), num_pixels_to_replace)
      number_replaced_per_image[i] = num_pixels_to_replace
      if num_pixels_to_replace gt 0 then begin
         d = darkarr[i,*,*]                                       ; work with just one image
         d[pixels_to_replace] = mean_of_images[pixels_to_replace] ; replace with mean
         darkarr[i,*,*] = temporary(d) ; reassign back into dark array
      endif
   endfor
  
   ; return a cleaned copy of imgarr
  return,darkarr
end
