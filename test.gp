#set term gif animate optimize delay 10 size 480,480
#set output 'movie.gif'

set pm3d at b
set xr[-5:10]
set yr[-5:10]
set zr[0:1]
set cbr[0:1]
set isosamples 50

do for [i = 0:50 ] {
   t=i*0.05
   splot sqrt(1/(1+t*t))*exp(-(x-t)**2/(1+t*t))*sqrt(1/(1+t*t))*exp(-(y-2*t)**2/(1+t*t))
   }

#set out
#set terminal wxt enhanced