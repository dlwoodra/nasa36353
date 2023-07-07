; include file for 36.353

; to use, specify as @config36353 in the code
; in IDL these lines are then executed in place in that code the same as if it was typed into there
; include files ("batch" files in IDL) cannot support loops or any multiline commands or code blocks

; define constants used all over the place

;version = '1_1' ; restrict to +/-41 sec around apogee 6/30/23
version = '1_2' ; dark adjust MEGS-B using flight L0D code
numberstr = '353' ; 36.353
theyd = 2021252 ; sept 9, 2021 17:25:00 UT
humandatestr='Sept 9, 2021'
au = 0.98589 ; from lisird lasp_vsop87_1au_correction_PT1M 9/9/21 17:30 UTC
; call this interactively to get 1-AU factor from LISIRD
; s=get_lisird_data(dataset='lasp_vsop87_1au_correction_PT1M',mintime='2021-09-09T17:30:00',maxtime='2021-09-09T17:35:00',/jd)
; earth is far from sun, so need to increase irradiance

apogeekmstr = '293.98'
apogeesecstr = '280.05'

; used in msis00e
ft7 = 101.1
ft7a = 88.
fap = 7.0 ; Frederiksburg

sza=31.05 ; deg

;MEGS-A indices
dark1idx_a=[81,82,84,89,90] ; remove_megsa_spikes
ffidx_a=[133,134]
solaridx_a=[94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129]
dark2idx_a=[130,131,136,137]

;MEGS-B indices
dark1idx_b=[81,83,84,86,89,90] ; remove_megsb_spikes
ffidx_b=[133,134]
solaridx_b=[94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129]
dark2idx_b=[130,138,139]
