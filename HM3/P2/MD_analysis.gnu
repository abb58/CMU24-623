# Plot to analyze the .ener file from the MD

plot "LDmj_sim.ener" u 2:7 lw 2 title "TE" w l, \
"LDmj_sim.ener" u 2:5 lw 2 title "PE" w l, \
"LDmj_sim.ener" u 2:6 lw 2 title "KE" w l
set xtics font "Times-Roman, 18"
set ytics font "Times-Roman, 18"
set xlabel "time" font "Times-Roman, 20" 
set ylabel "energy" font "Times-Roman, 20" 
set title "MD simulation (Energy) of Lennard-Jones liquid" font "Times-Roman, 30"
set term png size 1024, 768
set output "LDmj_sim_ener.png"
replot

plot "LDmj_sim.ener" u 2:8 lw 2 title "Px" w l, \
"LDmj_sim.ener" u 2:9 lw 2 title "Py", \
"LDmj_sim.ener" u 2:10 lw 2 title "Pz" w l
set xtics font "Times-Roman, 18"
set ytics font "Times-Roman, 18"
set xlabel "time" font "Times-Roman, 20" 
set ylabel "momentum" font "Times-Roman, 20" 
set title "MD simulation (Momentum) of Lennard-Jones liquid" font "Times-Roman, 30"
set term png size 1024, 768
set output "LDmj_sim_mom.png"
replot

set autoscale
plot "LDmj_sim.ener" u 2:3 lw 2 title "T*" w l
set xtics font "Times-Roman, 18"
set ytics font "Times-Roman, 18"
set xlabel "time" font "Times-Roman, 20" 
set ylabel "Temperature" font "Times-Roman, 20"
set title "MD simulation (T) of Lennard-Jones liquid" font "Times-Roman, 30"
set term png size 1024, 768
set output "LDmj_sim_T.png"
replot

set autoscale
plot "LDmj_sim.ener" u 2:4 lw 2 title "P*" w l
set xtics font "Times-Roman, 18"
set ytics font "Times-Roman, 18"
set xlabel "time" font "Times-Roman, 20" 
set ylabel "Pressure" font "Times-Roman, 20"
set title "MD simulation (P) of Lennard-Jones liquid" font "Times-Roman, 30"
set term png size 1024, 768
set output "LDmj_sim_P.png"
replot
