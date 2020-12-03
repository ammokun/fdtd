reset
set pm3d 
set pm3d map
set ticslevel 0
set cbrange[-1E5:1E5]
set palette defined (-1E5 "blue", 0 "white", 1E5 "red")
splot file u 1:2:3 with pm3d
set key bottom
