# xpshape
NMR eXchange Peak SHAPE simulation of the Bloch-McConnell equation

# Usage
<pre>
   Usage: xpshape [-model AB|ACB]
                  [-dw #(Hz)] [-kab #(1/s)] [-kba #(1/s)]
                  [-R20a #(1/s)] [-R20b #(1/s)]
                  [-a0 conc(M)] [-c0 conc(M)]
                  [-swh Hz(600)] [-fid pts(8192)]

   simulate NMR peakshape for two-site chemical exchange

   -model,AB  : A = B     (dw,kab,kba,R20a,R20b required)
          ACB : A + C = B (dw,kab,kba,R20a,R20b,a0,c0 required)

   -dw,  (A.frequency= -dw/2, B.frequency= +dw/2)
   -R20a,(exchange-free R2 of A state)
   -R20b,(exchange-free R2 of B state)
   -kab, (rate constant of A->B)
   -kba, (rate constant of B->A)
   -a0,  (initial concentration of A)
   -c0,  (initial concentration of C)
   -swh, (spectrum = -swh/2 ... +swh/2)
   -fid, (must be ...,1024,2048,4096,8192,16384,32768,65536,..)
</pre>

# Example
<pre>
### Kd=10mM, [P]=2uM, [L]=5-200uM
./xpshape -model ACB -dw 100 -kab 1e+6 -kba 1e+4 -R20a 0.1 -R20b 150 -a0 2e-6 -c0 5e-6 > tins.Kd.10mM.dat
./xpshape -model ACB -dw 100 -kab 1e+6 -kba 1e+4 -R20a 0.1 -R20b 150 -a0 2e-6 -c0 10e-6 >> tins.Kd.10mM.dat
./xpshape -model ACB -dw 100 -kab 1e+6 -kba 1e+4 -R20a 0.1 -R20b 150 -a0 2e-6 -c0 20e-6 >> tins.Kd.10mM.dat
./xpshape -model ACB -dw 100 -kab 1e+6 -kba 1e+4 -R20a 0.1 -R20b 150 -a0 2e-6 -c0 50e-6 >> tins.Kd.10mM.dat
./xpshape -model ACB -dw 100 -kab 1e+6 -kba 1e+4 -R20a 0.1 -R20b 150 -a0 2e-6 -c0 100e-6 >> tins.Kd.10mM.dat
./xpshape -model ACB -dw 100 -kab 1e+6 -kba 1e+4 -R20a 0.1 -R20b 150 -a0 2e-6 -c0 200e-6 >> tins.Kd.10mM.dat
### Kd=10uM, [P]=2uM, [L]=5-200uM
./xpshape -model ACB -dw 100 -kab 1e+6 -kba 1e+1 -R20a 0.1 -R20b 150 -a0 2e-6 -c0 5e-6   >  tins.Kd.10uM.dat
./xpshape -model ACB -dw 100 -kab 1e+6 -kba 1e+1 -R20a 0.1 -R20b 150 -a0 2e-6 -c0 10e-6  >> tins.Kd.10uM.dat
./xpshape -model ACB -dw 100 -kab 1e+6 -kba 1e+1 -R20a 0.1 -R20b 150 -a0 2e-6 -c0 20e-6  >> tins.Kd.10uM.dat
./xpshape -model ACB -dw 100 -kab 1e+6 -kba 1e+1 -R20a 0.1 -R20b 150 -a0 2e-6 -c0 50e-6  >> tins.Kd.10uM.dat
./xpshape -model ACB -dw 100 -kab 1e+6 -kba 1e+1 -R20a 0.1 -R20b 150 -a0 2e-6 -c0 100e-6 >> tins.Kd.10uM.dat
./xpshape -model ACB -dw 100 -kab 1e+6 -kba 1e+1 -R20a 0.1 -R20b 150 -a0 2e-6 -c0 200e-6 >> tins.Kd.10uM.dat
</pre>
