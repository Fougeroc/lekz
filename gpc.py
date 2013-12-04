# we first check if we are in sage. If we are, we avoid reimplementing highly
# optimized stuff.
IN_SAGE = False
try:
    import sage
    IN_SAGE = True
except ImportError:
    pass

if not IN_SAGE:
    def gcd(p,q):
        r"""
        Computation of the gcd using Euclid algorithm.

        EXAMPLES:

            sage: gcd(2,4)
            2
            sage: gcd(12,6)
            6
            sage: gcd(3,5)
            1
        """
        if p < q:
            p,q = q,p
        while q != 0:
            p,q = q,p%q
        return p
else:
    from sage.rings.arith import gcd

def mean_and_std_dev(l):
    r"""
    Return the mean and standard deviation of the floatting point numbers in
    the list l.

    The implementation is very naive and should not be used for large list
    (>1000) of numbers.

    .. NOTE::

    mean and std are implemented in Sage but are quite buggy!
    """
    from math import sqrt
    m = sum(l) / len(l)
    if len(l) == 1:
        d = 0
    else:
        d = sum((x-m)**2 for x in l) / (len(l)-1)
    return m,sqrt(d)

def perm_to_twin(p):
    r"""
    Return the twin associated to the permutation ``p``.

    EXAMPLES::

        sage: import gpc
        sage: gpc.perm_to_twin([[0,0,1],[1,2,2]])
        [[(0, 1), (0, 0), (1, 0)], [(0, 2), (1, 2), (1, 1)]]
    """
    d = {}
    for k in xrange(2):
        for i,j in enumerate(p[k]):
            if j not in d:
                d[j] = []
            d[j].append((k,i))

    t = [[None] * len(p[0]), [None] * len(p[1])]
    for l in d.itervalues():
        assert(len(l) == 2)
        (k0,i0),(k1,i1) = l
        t[k0][i0] = (k1,i1)
        t[k1][i1] = (k0,i0)

    return t

def twin_canonical_orientation(twin):
    r"""
    Return a choice of orientation. Each letters in a pair is paired with one +1
    and one -1.

    TESTS::

        sage: gp = iet.GeneralizedPermutation('a a b','b c c')
        sage: canonical_orientation_from_twin(gp._twin)
        ([1, -1, 1], [-1, -1, 1])

        sage: gp = iet.GeneralizedPermutation('a a','b b c c')
        sage: canonical_orientation_from_twin(gp._twin)
        ([1, -1], [-1, 1, -1, 1])
    """
    or_top = []
    for i,(j,k) in enumerate(twin[0]):
        if j == 1 or i < k:
            or_top.append(1)
        else:
            or_top.append(-1)
    or_bot = []
    for i,(j,k) in enumerate(twin[1]):
        if j == 0 or i < k:
            or_bot.append(-1)
        else:
            or_bot.append(1)
    return or_top, or_bot

# formatting function
or_and_label_to_str = lambda l,o: ' %s'%l if o == 1 else '-%s'%l
nice_deg_str = lambda d,m: str(d) if m == 1 else str(d) + '^' + str(m)

def stratum_str(angles):
    r"""
    Return a nice string for the stratum from a list of angles.

    EXAMPLES::

        sage: print stratum_str([1,1,1,1])
        Q_0(-1^4)
        sage: print stratum_str([6,1,1,1,1])
        Q_1(4, -1^4)
    """
    degrees = {}
    for a in angles:
        d = a-2
        if d not in degrees:
            degrees[d] = 0
        degrees[d] += 1
    g = sum(m*d for (m,d) in degrees.iteritems()) // 4 + 1

    return "Q_%d("%g + ', '.join(nice_deg_str(d,degrees[d]) for d in sorted(degrees,reverse=True)) + ')'



def order_mod_n(k,n):
    r"""
    Return the (additive) order of ``k`` mod ``n`` which is ``n/gcd(k,n)``.

    EXAMPLES::

        sage: order_mod_n(0,2)
        1
        sage: order_mod_n(1,2)
        2
        sage: [order_mod_n(i,4) for i in xrange(4)]
        1 4 2 4
    """
    if k == 0:
        return 1
    return n // gcd(k,n)

class GeneralizedPermutationCyclicCover:
    r"""
    A cyclic cover of a generalized permutation.

    This class is used to model a an affine measure in the moduli space of
    translation surfaces that is obtained by cyclic cover of a Masur-Veech one.
    """
    def __init__(self, p, sigma=None, degree=None):
        r"""
        INPUT:

        - ``p`` -- a permutation. The labels used must be the first ``n``
          integers starting at ``0``.

        - ``sigma`` -- an optional covering data that describe the gluings
          for each interval. This data must be a list of positive integers.

        - ``degree`` -- an optional degree for the covering.
        """
        self._p = [map(int,p[0]), map(int,p[1])]
        self._n = (len(p[0]) + len(p[1])) // 2

        assert all(i >= 0 and i < self._n for i in p[0])
        assert all(i >= 0 and i < self._n for i in p[1])
        self._twin = perm_to_twin(p)
        self._orientation = twin_canonical_orientation(self._twin)

        if sigma is None: # not a cover
            assert degree is None
            sigma = [0] * self._n
            degree = 1

        self._sigma = map(int,sigma)
        self._degree = int(degree)

        assert len(sigma) == self._n
        for i in self._sigma:
            assert 0 <= i and i < degree

    def __repr__(self):
        o0 = self._orientation[0]
        p0 = self._p[0]
        ntop = len(p0)
        o1 = self._orientation[1]
        p1 = self._p[1]
        nbot = len(p1)
        return (' '.join(or_and_label_to_str(p0[i],o0[i]) for i in xrange(ntop))
                + '\n' +
                ' '.join(or_and_label_to_str(p1[i],o1[i]) for i in xrange(nbot)))

    def vertices_quotient(self):
        r"""
        Return the list of vertices in the quotient. Each vertex is modelized by
        the list of edges that are crossed by considering a small loop around
        that vertex. Each edge is identified by a pair ``(label,orientation)``.
        """
        p = self._p
        o = self._orientation

        ntop = len(p[0])
        nbot = len(p[1])
        go_right = {}
        for i in xrange(ntop-1):
            go_right[(p[0][i],o[0][i])] = (p[0][i+1],o[0][i+1])
        go_right[(p[0][-1],o[0][-1])] = (p[1][-1],o[1][-1])
        for i in xrange(nbot-1):
            go_right[(p[1][i+1],o[1][i+1])] = (p[1][i],o[1][i])
        go_right[(p[1][0],o[1][0])] = (p[0][0],o[0][0])

        cycles = []
        wait = set([(l,1) for l in range(self._n)] + [(l,-1) for l in range(self._n)])
        while wait:
            l = wait.pop()
            cycle = [l]
            ll = go_right[(l[0],-l[1])]
            while ll != l:
                cycle.append(ll)
                wait.remove(ll)
                ll = go_right[(ll[0],-ll[1])]
            cycles.append(cycle)

        return cycles

    def angles_and_ramified_angles(self):
        r"""
        Return the a pair of lists that are respectively the list of angles of
        the stratum in the quotient and the list of angles in the ambient
        stratum.

        EXAMPLES::

            sage: p = gpc.GeneralizedPermutationCyclicCover([[0,1,1],[2,2,3,3,4,4,0]],[1,0,1,1,0],2)
            sage: p.angles_and_ramified_angles()
            ([3, 1, 1, 1, 1, 1], [6, 2, 1, 1, 1, 1, 2, 2])
        """
        p = self._p
        o = self._orientation
        cycles = self.vertices_quotient()

        left = (p[0][0],o[0][0])
        right = (p[1][-1],o[1][-1])
        angles = []
        for cycle in cycles:
            angles.append(sum(1 for l in cycle if l != left and l != right))

        ramified_angles = []
        for a,cycle in zip(angles,cycles):
            monodromy = sum(o * self._sigma[l] for (l,o) in cycle)
            k = order_mod_n(monodromy,self._degree)
            ramified_angles.extend([k*a] * (self._degree // k))

        return angles, ramified_angles

    def genus(self):
        r"""
        Return the genus of the cover.
        """
        _,angles = self.angles_and_ramified_angles()
        return sum(a-2 for a in angles)/4 + 1

    def locus_string(self):
        r"""
        Return a string that describe the locus.
        """
        angles, ramified_angles = self.angles_and_ramified_angles()
        if self._degree == 1:
            return stratum_str(angles)
        else:
            return stratum_str(angles) + " --> " + stratum_str(ramified_angles)

    def print_locus(self):
        r"""
        Return the covering locus.

        Warning: there is no check of orientability. Hence Abelian strata
        appears as components of quadratic strata here!

        EXAMPLES:

        The orientation cover of the pillocase::

            sage: GPCC = GeneralizedPermutationCyclicCover
            sage: p = GPCC([[1,0,0],[2,2,1]],[1,0,1],2)
            sage: p.locus()
            Q_0(-1^4) --> Q_1(0^4)

        The locus in genus 2 associated to the windtree::

            sage: p = GPCC([[0,1,1],[2,2,3,3,4,4,0]],[0,1,1,1,0],2)
            sage: p.locus()
            Q_0(-1^5, 1) --> Q_1(-1, 0^4, 1)
        """
        print self.locus_string()

    def lyapunov_exponents_H_plus(self, nb_vectors=None, nb_experiments=100,
            nb_iterations=32768, verbose=True, output_file=None):
        r"""
        Compute the H^+ Lyapunov exponents for this covering locus.

        It calls the C-library lyap_exp interfaced with Cython. The computation
        might be significantly faster if ``nb_vectors=1`` (or if it is not
        provided but genus is 1).

        INPUT:

        - ``nb_vectors`` -- the number of exponents to compute. The number of
          vectors must not exceed the dimension of the space!

         - ``nb_experiments`` -- the number of experiments to perform. It might
           be around 100 (default value) in order that the estimation of
           confidence interval is accurate enough.

         - ``nb_iterations`` -- the number of iteration of the Rauzy-Zorich
           algorithm to perform for each experiments. The default is 2^15=32768
           which is rather small but provide a good compromise between speed and
           quality of approximation.

        - ``verbose`` -- if ``True`` provide additional informations rather than
          returning only the Lyapunov exponents (i.e. ellapsed time, confidence
          intervals, ...)

        - ``output_file`` -- if provided (as a file object or a string) output
          the additional information in the given file rather than on the
          standard output.
        """
        import time
        import lekz    # the cython bindings
        if nb_vectors is None:
            nb_vectors = self.genus()

        if output_file is None:
            from sys import stdout
            output_file = stdout
        elif isinstance(output_file, str):
            output_file = open(output_file, "w")

        nb_vectors = int(nb_vectors)
        nb_experiments = int(nb_experiments)
        nb_iterations = int(nb_iterations)

        if verbose:
            output_file.write(self.locus_string())
            output_file.write("\n")

        if nb_vectors <= 0:
            raise ValueError("the number of vectors must be positive")
        if nb_experiments <= 0:
            raise ValueError("the number of experiments must be positive")
        if nb_iterations <= 0:
            raise ValueError("the number of iterations must be positive")

        gp = self._p[0] + self._p[1]
        k = len(self._p[0])
        n = self._n
        twin = [interval * k + position for (interval,position) in self._twin[0] + self._twin[1]]

        t0 = time.time()
        res = lekz.lyapunov_exponents_H_plus_cyclic_cover(
                   gp, twin, k, n, self._sigma, self._degree,
                   nb_vectors, nb_experiments, nb_iterations)
        t1 = time.time()

        res_final = []

        m,d = mean_and_std_dev(res[0])
        if verbose:
            from math import log, floor, sqrt
            output_file.write("sample of %d experiments\n"%nb_experiments)
            output_file.write("%d iterations (~2^%d)\n"%(
                    nb_iterations,
                    floor(log(nb_iterations) / log(2))))
            output_file.write("ellapsed time %s\n"%time.strftime("%H:%M:%S",time.gmtime(t1-t0)))
            output_file.write("Lexp Rauzy-Zorich: %f (std. dev. = %f, conf. rad. 0.01 = %f)\n"%(
                    m,d, 2.576*d/sqrt(nb_experiments)))
        for i in xrange(1,nb_vectors+1):
            m,d = mean_and_std_dev(res[i])
            if verbose:
                output_file.write("theta%d           : %f (std. dev. = %f, conf. rad. 0.01 = %f)\n"%(
                    i,m,d, 2.576*d/sqrt(nb_experiments)))
            res_final.append(m)

        return res_final


def square_tiled_cyclic_cover(a0,a1,a2,a3,N):
    r"""
    Return the permutation associated to the square tiled cyclic cover
    M_N(a0,a1,a2,a3).
    """
    raise NotImplementedError
    if (a0+a1+a2+a3) %N:
        raise ValueError("a0 + a1 + a2 + a3 is not congruent to 0 mod N")

    return GeneralizedPermutationCyclicCover([[0,1,1],[2,2,0]], [], N)
