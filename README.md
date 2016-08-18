# xpshape
NMR eXchange Peak SHAPE simulation of the Bloch-McConnell equation

![Bloch-McConnell equation](https://github.com/sunghunbae/xpshape/raw/master/reference/Bloch-McConnell.png)

### References
1. Ulrich L. Gunther and Brian Schaffhausen (2002), J. Biomol. NMR 22, 201-209. NMRKIN: Simulating line shapes from two-dimensional spectra of proteins upon ligand binding. 
2. Hiroshi Matsuo et al. (1999) J. Am. Chem. Soc 121, 9903-9904. Identification by NMR Spectroscopy of Residues at Contact Surfaces in Large, Slowly Exchanging Macromolecular Complexes
3. Harden M. McConnell (1958) J. Chem. Physics. 28, 430-431. Reaction Rates by Nuclear Magnetic Resonance

### Usage
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

### Example: protein ligand binding
<pre>
# Kd=10mM, [P]=2uM, [L]=5-200uM
./xpshape -model ACB -dw 100 -kab 1e+6 -kba 1e+4 -R20a 0.1 -R20b 150 -a0 2e-6 -c0 5e-6 > tins.Kd.10mM.dat
./xpshape -model ACB -dw 100 -kab 1e+6 -kba 1e+4 -R20a 0.1 -R20b 150 -a0 2e-6 -c0 10e-6 >> tins.Kd.10mM.dat
./xpshape -model ACB -dw 100 -kab 1e+6 -kba 1e+4 -R20a 0.1 -R20b 150 -a0 2e-6 -c0 20e-6 >> tins.Kd.10mM.dat
./xpshape -model ACB -dw 100 -kab 1e+6 -kba 1e+4 -R20a 0.1 -R20b 150 -a0 2e-6 -c0 50e-6 >> tins.Kd.10mM.dat
./xpshape -model ACB -dw 100 -kab 1e+6 -kba 1e+4 -R20a 0.1 -R20b 150 -a0 2e-6 -c0 100e-6 >> tins.Kd.10mM.dat
./xpshape -model ACB -dw 100 -kab 1e+6 -kba 1e+4 -R20a 0.1 -R20b 150 -a0 2e-6 -c0 200e-6 >> tins.Kd.10mM.dat
# Kd=10uM, [P]=2uM, [L]=5-200uM
./xpshape -model ACB -dw 100 -kab 1e+6 -kba 1e+1 -R20a 0.1 -R20b 150 -a0 2e-6 -c0 5e-6   >  tins.Kd.10uM.dat
./xpshape -model ACB -dw 100 -kab 1e+6 -kba 1e+1 -R20a 0.1 -R20b 150 -a0 2e-6 -c0 10e-6  >> tins.Kd.10uM.dat
./xpshape -model ACB -dw 100 -kab 1e+6 -kba 1e+1 -R20a 0.1 -R20b 150 -a0 2e-6 -c0 20e-6  >> tins.Kd.10uM.dat
./xpshape -model ACB -dw 100 -kab 1e+6 -kba 1e+1 -R20a 0.1 -R20b 150 -a0 2e-6 -c0 50e-6  >> tins.Kd.10uM.dat
./xpshape -model ACB -dw 100 -kab 1e+6 -kba 1e+1 -R20a 0.1 -R20b 150 -a0 2e-6 -c0 100e-6 >> tins.Kd.10uM.dat
./xpshape -model ACB -dw 100 -kab 1e+6 -kba 1e+1 -R20a 0.1 -R20b 150 -a0 2e-6 -c0 200e-6 >> tins.Kd.10uM.dat
</pre>


### Example output
<pre>
# peak shape simulation [A + C = B]
# dw     = 628.319 (rad/s)
# kab    = 1e+06 (1/s)
# kba    = 10000 (1/s)
# R20a   = 0.1 (1/s)
# R20b   = 150 (1/s)
# kon    = 1e+06 (1/s)
# koff   = 10000 (1/s)
# Kd     = 0.01 (M)
# pb     = 0.00049965
# [A]ini = 2e-06 (M)
# [C]ini = 5e-06 (M)
# [A]    = 1.999e-06 (M)
# [C]    = 4.999e-06 (M)
# [B]    = 9.99301e-10 (M)
# FID    = 8192 (complex points)
# SWH    = 600 (Hz)
# X-axis unit: (Hz)
-300 0.00512234
-299.927 0.00512234
-299.854 0.00512233
-299.78 0.00512232
...... omitted ......
</pre>
![output plot using XMGRACE](https://github.com/sunghunbae/xpshape/raw/master/example/example.png)
