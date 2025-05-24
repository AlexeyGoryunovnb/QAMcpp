set terminal pngcairo enhanced font 'arial,10' fontscale 1.0 size 700, 500
set output 'DispPError.png'
set title 'SNR vs Bit Error Probability for QAM (p=3)'
set xlabel 'SNR (dB)'
set ylabel 'Bit Error Probability'
set logscale y
set grid
set key top right
plot 'qam_results(QAM256).dat' with linespoints title 'Simulated Bit Error', \
     15*erfc(sqrt(10**(x/10)))*(2**(3+1)-15*erfc(sqrt(10**(x/10)))/2)/(4**4) title 'Theoretical Upper Bound', \
     15*erfc(sqrt(10**(x/10)))*(2**(3+1)-15*erfc(sqrt(10**(x/10)))/2)/(4**4)/4 title 'Theoretical Lower Bound', \
