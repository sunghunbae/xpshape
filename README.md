# xpshape
NMR eXchange Peak SHAPE (simulation of Bloch-McConnell equation)

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
