Comments on GPT2(/GMF2)

This subroutine determines pressure, temperature, temperature lapse rate, water 
vapour pressure, hydrostatic and wet mapping function coefficients ah and aw, 
and geoid undulation for specific sites near the Earth surface. It is based on 
a 5 x 5 degree external grid file ('gpt2_5.grd') with mean values as well as 
sine and cosine amplitudes for the annual and semiannual variation of the 
coefficients.

Apart from the better resolution (in space and time) and the larger number of 
parameters, there are some technical changes which do not make the transition 
from the old GPT/GMF straightforward:

There is an external grid file which – in the present form of the subroutine – 
has to be in the same directory as the subroutine. This can easily be changed in 
the subroutine.

With the new subroutine, the mapping function coefficients are already provided 
by GPT2. (This needs to be confirmed by testing.) These mapping function 
coefficients are then input to the subroutine for the gridded Vienna Mapping 
Function, i.e. to vmf1_ht.for. Thus, there is no ‘GMF2’ but just a GPT2.

Basically, GPT2 has to be called only once per 24 hour session, since it only 
contains annual and semiannual variations of the parameters. That’s why the 
subroutine allows as input just one epoch (modified Julian date) but many 
stations (the length of the vectors is presently set to 64800 but could be 
changed as necessary). The idea behind that is that the gridfile has to be 
loaded only once per epoch.

