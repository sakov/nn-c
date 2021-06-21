## nn ##
(Natural Neighbours interpolation)

**nn** is a C code for Natural Neighbours interpolation of 2D scattered data. It provides a C library and a command line 
utility **nnbathy**. 

Algorithmically, it was initially loosely based on the Dave Watson's description of nngridr; code-wise it is an 
independent development. You may see a comparison of performance of a (rather old) version of **nn** with nngridr in

  Quanfu Fan, Alon Efrat, Vladlen Koltun, Shankar Krishnan, and Suresh 
  Venkatasubramanian. Hardware-assisted Natural Neighbor Interpolation. 
  In Proc. 7th Workshop on Algorithm Engineering and Experiments (ALENEX), 2005.
  [pdf](http://nn-c.googlecode.com/files/fan05a.pdf)

**nn** is coded for robustness (to handle degenerate data) and scalability (to handle millions of data points), subject 
to using double precision calculations. For the underlying Delaunay triangulation it calls exact arithmetic code from 
[triangle](http://www.cs.cmu.edu/~quake/triangle.html). From v2 it is possible to run **nnbathy** on multiple CPUs with
triangulation stored in shared memory.


Checkout **nn** by running `git clone https://github.com/sakov/nn-c` or `svn checkout https://github.com/sakov/nn-c`.
