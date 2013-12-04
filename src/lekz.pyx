r"""
Python bindings for various computation of Lyapunov exponents.
"""

from libc.stdlib cimport malloc,free

cdef extern from "lyapunov_exponents.h":
    ctypedef struct quad_cyclic_cover:
        pass
    ctypedef struct generalized_permutation:
        pass

    # initialisation/allocation/free
    generalized_permutation * new_generalized_permutation(int *perm, int *twin, int k, int n)
    quad_cyclic_cover * new_quad_cyclic_cover(generalized_permutation * gp, int * sigma, size_t degree, size_t nb_vectors)
    void set_random_lengths_quad_cyclic_cover(quad_cyclic_cover * qcc)
    void set_random_vectors(quad_cyclic_cover * qcc)

    #int check_generalized_permutation(generalized_permutation *p)
    #int check_quad_cyclic_cover(quad_cyclic_cover * qcc)

    void free_generalized_permutation(generalized_permutation ** gp)
    void free_quad_cyclic_cover(quad_cyclic_cover ** qcc)

    # print
    #void print_generalized_permutation(generalized_permutation * p)
    void print_quad_cyclic_cover(quad_cyclic_cover * qcc)
    void print_vectors(quad_cyclic_cover * qcc)

    # algorithms
    #void renormalize_length_quad_cyclic_cover(quad_cyclic_cover * qcc)
    #void rauzy_induction_H_plus_quad_cyclic_cover(quad_cyclic_cover *qcc)

    int init_GS(size_t dim)
    void free_GS()
    #void orthogonalize_GS(quad_cyclic_cover * qcc, double * theta)

    void lyapunov_exponents_H_plus(quad_cyclic_cover *qcc, double *theta, size_t nb_induction)
    void top_lyapunov_exponents_H_plus(quad_cyclic_cover *qcc, double *theta, size_t nb_iterations)

def lyapunov_exponents_H_plus_cyclic_cover(
        gp, twin, k, n, sigma, degree,
        nb_vectors, nb_experiments, nb_iterations):
    r"""
    Compute the Lyapunov exponents of the H^+ part of the KZ-cocycle for
    covering locii.

    We assume that all the inputs are clean. If not, it may cause some SEGFAULT
    which would interrupt python!

    INPUT:

    - ``gp`` -- a generalized permutation given as a list of integers

    - ``twin`` -- the twin data of the gp

    - ``k`` -- the length of the top interval

    - ``n`` -- the length of gp

    - ``sigma`` -- covering data

    - ``nb_vectors`` -- the number of vectors to use

    - ``nb_experiments`` -- number of experimets

    - ``nb_iterations`` -- the number of iterations of the Rauzy-Zorich
      induction to perform

    - ``verbose`` -- if ``True`` print additional information concerning the
      mean and standard deviation
    """
    cdef int *p, *t, *s  # permutation, twin, sigma
    cdef generalized_permutation *gp_c
    cdef quad_cyclic_cover *qcc
    cdef double * theta


    # convert the data of into C values
    p = <int *> malloc(2 * n * sizeof(int))
    t = <int *> malloc(2 * n * sizeof(int))
    s = <int *> malloc(n * sizeof(int))

    for i from 0 <= i < n:
        p[i] = gp[i]
        t[i] = twin[i]
        s[i] = sigma[i]
    for i from n <= i < 2*n:
        p[i] = gp[i]
        t[i] = twin[i]
    theta = <double *> malloc((nb_vectors+1) * sizeof(double))

    gp_c = new_generalized_permutation(p, t, k, n)
    qcc = new_quad_cyclic_cover(gp_c,s,degree,nb_vectors)

    free_generalized_permutation(&(gp_c))
    free(p)
    free(t)
    free(s)

    res = [[] for _ in xrange(nb_vectors+1)]
    if nb_vectors == 1:
        for i in xrange(nb_experiments):
            top_lyapunov_exponents_H_plus(qcc, theta, nb_iterations)
            for j in xrange(2):
                res[j].append(theta[j])

    else:
        init_GS(nb_vectors)
        for i in xrange(nb_experiments):
            lyapunov_exponents_H_plus(qcc, theta, nb_iterations)
            for j in xrange(nb_vectors+1):
                res[j].append(theta[j])

    free_quad_cyclic_cover(&qcc)
    free(theta)

    return res
