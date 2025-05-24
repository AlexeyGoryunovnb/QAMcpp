set terminal pngcairo enhanced font 'arial,10' fontscale 1.0 size 700, 500
set output 'DispPError.png'
set title 'SNR vs Bit Error Probability for QAM (p=2)'
set xlabel 'SNR (dB)'
set ylabel 'Bit Error Probability'
set logscale y
set grid
set key top right
plot 'qam_results(QAM64).dat' with linespoints title 'Simulated Bit Error', \
     7*erfc(sqrt(10**(x/10)))*(2**(2+1)-7*erfc(sqrt(10**(x/10)))/2)/(4**3) title 'Theoretical Upper Bound', \
     7*erfc(sqrt(10**(x/10)))*(2**(2+1)-7*erfc(sqrt(10**(x/10)))/2)/(4**3)/3 title 'Theoretical Lower Bound', \

