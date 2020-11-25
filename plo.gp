reset
set pm3d 
set pm3d map
set ticslevel 0
set cbrange[-2E6:2E6]
set palette defined (-2E6 "blue", 0 "white", 2E6 "red")
splot file with pm3d
set key bottom
