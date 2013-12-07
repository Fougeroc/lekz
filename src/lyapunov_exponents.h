#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#define EPSILON_LENGTH_PROJECTION .0000000000000001L
#define LENGTH_ERROR_TOLERANCE    .000000000001L
# define LOG2l     0.6931471805599453094172321214581766L

//#define USE_KAHAN_SUMMATION

/* vectors convention                                            */
/* 0 <= i < nb_vectors                                           */
/* 0 <= j < nb_intervals                                         */
/* 0 <= k < degree                                               */
/* elt at pos (i,j,k) = v[k + degree * (j + nb_intervals * i)]   */

/* actually, it would be good to store the teichmueller time */
/* inside qcc and transparently renormalize the lengths      */
/* depending on lbot and ltop                                */

/****************************/
/* custom random functions */
/***************************/

inline double drand(void);
inline long double ldrand(void);

/***************************************************/
/* data and functions for generalized permutations */
/***************************************************/

typedef struct Xgeneralized_permutation{
int * perm;
int * twin;
int k,n;
} generalized_permutation;

generalized_permutation * new_generalized_permutation(int *perm, int *twin, int k, int n);
int check_generalized_permutation(generalized_permutation *p);
void free_generalized_permutation(generalized_permutation ** gp);
void print_generalized_permutation(generalized_permutation * p);

/**************************************************************/
/* data and functions cyclic cover of quadratic differentials */
/**************************************************************/

typedef struct Xlabel{
long double length;      /* length of the subinterval                                                 */
int same_interval;       /* a boolean that tells us if the two intervals belong to the same interval  */
int sigma;               /* value of the group element on that label (seen as a pi_1 representation)  */
double * v;              /* a vector of size degree x nb_vectors */
} label;

/* Note: the interval datatype does not depend on the extension we take!! */
typedef struct Xinterval{
int orientation;                /* orientation of the interval (0 or 1) */
label * lab;                    /* all the data for that interval       */
struct Xinterval *twin;         /* the twin interval                    */
struct Xinterval *prev, *next;  /* the guy on the left and on the right */
} interval;

typedef struct{
size_t nb_labels;        /* number of labels                              */
size_t degree;           /* degree of the cover                           */
size_t nb_vectors;       /* number of vectors in use                      */
interval *top,*bot;      /* the leftmost intervals                        */
label * labels;          /* array of labels                               */
interval * intervals;    /* array of intervals                            */
long double length;      /* length of the top and bot intervals           */
double * buffer;         /* a buffer of size degree x nb_vectors          */
} quad_cyclic_cover;

quad_cyclic_cover * new_quad_cyclic_cover(generalized_permutation * gp, int * sigma, size_t degree, size_t nb_vectors);
void free_quad_cyclic_cover(quad_cyclic_cover ** qcc);
int check_quad_cyclic_cover(quad_cyclic_cover * qcc);
void print_quad_cyclic_cover(quad_cyclic_cover * qcc);
void set_random_lengths_quad_cyclic_cover(quad_cyclic_cover * qcc);
void renormalize_length_quad_cyclic_cover(quad_cyclic_cover * qcc);
void randomize_length_quad_cyclic_cover(quad_cyclic_cover * qcc);

void lyapunov_exponents_H_plus(quad_cyclic_cover *qcc, double *theta, size_t nb_induction);
void top_lyapunov_exponents_H_plus(quad_cyclic_cover *qcc, double *theta, size_t nb_iterations);

/* one step of Rauzy induction with Zorich acceleration */
/* depending on the method chosen, the vectors are updated differently */
void rauzy_induction_quad_cyclic_cover(quad_cyclic_cover *qcc);
void rauzy_induction_H_plus_quad_cyclic_cover(quad_cyclic_cover *qcc);

/***************************/
/* a bit of linear algebra */
/***************************/

void set_random_vectors(quad_cyclic_cover * qcc);
void print_vectors(quad_cyclic_cover * qcc);

int init_GS(size_t dim);
void free_GS(void);
void orthogonalize_GS(quad_cyclic_cover * qcc, double * theta);
void check_orthogonality(quad_cyclic_cover * qcc);

