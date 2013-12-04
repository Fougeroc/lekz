Lyapunov exponents of the KZ cocycle
====================================

Algorithms to compute numerical approximations of the Lyapunov exponents of
the Kontsevich-Zorich cocycle for various SL(2,R)-invariant measures.

Installation and usage
----------------------

To build the package, do:

    $ python setup.py build_ext --inplace

To use it, simply launch Python/IPython/Sage in that repository and then

    >>> import gpc
    >>>  p = gpc.GeneralizedPermutationCyclicCover([[0,0,1],[1,2,2]], [7,0,5], 10)
    >>> p.print_locus()
    Q_0(-1^4) --> Q_5(8^2, 0^10)
    >>> p.lyapunov_exponents_H_plus()
    Q_0(-1^4) --> Q_5(8^2, 0^10)
    sample of 100 experiments
    32768 iterations (~2^15)
    ellapsed time 00:01:08
    Lexp Rauzy-Zorich: 0.546868 (std. dev. = 0.002734, conf. rad. 0.01 = 0.000704)
    theta1           : 0.999892 (std. dev. = 0.000064, conf. rad. 0.01 = 0.000017)
    theta2           : 0.407110 (std. dev. = 0.003116, conf. rad. 0.01 = 0.000803)
    theta3           : 0.372974 (std. dev. = 0.003146, conf. rad. 0.01 = 0.000810)
    theta4           : 0.326335 (std. dev. = 0.003030, conf. rad. 0.01 = 0.000781)
    theta5           : 0.285807 (std. dev. = 0.003059, conf. rad. 0.01 = 0.000788)


Structure
--------
The file gpc.py contains a simply Python file that implements combinatorics of
covers. The rest of the code is contained in the repository src/. The file
src/lekz.pyx contains the binding to the C-code and all others (
lyapunov\_exponents.h, generalized\_permutation.c, quad\_cyclic\_cover.c and random.c)
contain the C-code.

Todo
----
Explore the connected components of covers. For square tiled surface, there
seems to be a finitude of primitive cover. There might be the same phenomenon
in general.
