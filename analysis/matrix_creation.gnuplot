set pm3d map
set autoscale xfix
set autoscale yfix
unset key
#set terminal pngcairo  transparent enhanced font "arial,16" fontscale 1.0 size 600, 400
set terminal pngcairo enhanced font "arial,12" fontscale 1.0 size 600, 400
set out "Matrix-corr.png"
set ytics out nomirror
set xtics out nomirror
set mxtics; set mytics
set border lw 3
set cbrange [-1:1]
set cbtics 0.2 in nomirror 
set palette defined (-1 "#FF00FF", -0.8 "#000080", -0.6 "#00FFFF", -0.4 "#808080", -0.2 "#FFFFFF", 0 "#FFFFFF", 0.2 "#FFFFFF", 0.4 "#808080", 0.6 "#008000", 0.8 "#FFFF00", 1 "#FF0000" )
plot "matrix.dat" matrix with image

#stats "matrix_correl_C1.dat" matrix using 3 nooutput
#cbmax = (abs(STATS_min) > abs(STATS_max) ? abs(STATS_min) : abs(STATS_max))
#set cbrange [-cbmax:cbmax]
quit
