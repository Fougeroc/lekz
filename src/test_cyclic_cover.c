/********************************************************************/
/* instruction to add test example:                                 */
/*  1 - increase the macro NB_TEST by 1                             */
/*  2 - follow the examples to define pj,tj,sj,dj,kj,nj (where j is */
/*      the appropriate integer)                                    */             
/*  3 - in the function init_examples, update the variables P,T,S,D */
/*      K and N to contain the new value                            */
/********************************************************************/

#include "lyapunov_exponents.h"

#define NB_TEST 5
#define NB_VECTORS 3

/* example Q(-1,-1,-1,-1) with no cover                               */
/* no theta+                                                          */
static int p0[6] = {0,0,1,1,2,2};   /* permutation data               */
static int t0[6] = {1,0,3,2,5,4};   /* twin data, i.e. p[i] = p[t[i]] */
static int s0[3] = {0,0,0};         /* covering data                  */
static int d0 = 1;                  /* degree of the cover            */
static int k0 = 3;                  /* length of the top interval     */
static int n0 = 3;                  /* number of labels               */

/* example Q(-1,-1,-1,-1) whose double cover is H(0)                  */
/* theta+ = 1                                                         */
static int * p1 = p0;
static int * t1 = t0;
static int s1[3] = {1,0,1};
static int d1 = 2;
static int k1 = 3;
static int n1 = 3;

/* example in Q(2,-1,-1) whose double cover is in H(1,1)                 */
/* for each label, the value of sigma should coincide with same_interval */
/* theta+ = 1, 1/2                                                       */
static int p2[8] = {0,1,2,2,1,3,3,0};
static int t2[8] = {7,4,3,2,1,6,5,0};
static int s2[4] = {0,0,1,1};
static int d2 = 2;
static int k2 = 4;
static int n2 = 4;

/* example in Q(1,-1^5) */
/* - orientation cover (end up in H(2))            */
/*   theta+ = 1, 1/3                               */
/* - ramified over 4 poles (end up in Q(1^2,-1^2)) */
/*   theta+ = 2/3?                                 */
static int p3[10] = {1,0,0,2,2,3,3,4,4,1};  static int *p4 = p3;
static int t3[10] = {9,2,1,4,3,6,5,8,7,0};  static int *t4 = t3;
static int s3[5]  = {1,0,1,1,1};            static int s4[5] = {1,0,1,1,0};
static int d3 = 2, d4 = 2;
static int k3 = 3, k4 = 3;
static int n3 = 5, n4 = 5;

static int * P[NB_TEST];
static int * T[NB_TEST];
static int * S[NB_TEST];
int D[NB_TEST];
int K[NB_TEST];
int N[NB_TEST];

generalized_permutation * gp;
quad_cyclic_cover * qcc;

void init_examples(void)
{
	P[0] = p0; P[1] = p1; P[2] = p2; P[3] = p3; P[4] = p4;
	T[0] = t0; T[1] = t1; T[2] = t2; T[3] = t3; T[4] = t4;
	S[0] = s0; S[1] = s1; S[2] = s2; S[3] = s3; S[4] = s4;
	D[0] = d0; D[1] = d1; D[2] = d2; D[3] = d3; D[4] = d4;
	K[0] = k0; K[1] = k1; K[2] = k2; K[3] = k3; K[4] = k4;
	N[0] = n0; N[1] = n1; N[2] = n2; N[3] = n3; N[4] = n4;
}

int test_allocation(int *p, int *t, int *s, int d, int k, int n)
{
	size_t i;

	fprintf(stdout, "TEST: allocation\n");
   	fprintf(stdout, "================\n");
	fprintf(stdout,"Test malloc/free for generalized permutation...");
	fflush(stdout);
	for(i=0; i<100; ++i)
	{
		gp = new_generalized_permutation(p,t,k,n);
		if(gp == NULL)
		{
			fprintf(stderr, "failed allocation\n");
			return 1;
		}
		if(check_generalized_permutation(gp))
		{
			fprintf(stderr, "error\n");
			return 1;
		}
		free_generalized_permutation(&gp);
	}
	fprintf(stdout, "done\n");

	fprintf(stdout, "Test malloc/free for quadratic cover...");
	fflush(stdout);
	for(i=0; i<100; ++i)
	{
		gp = new_generalized_permutation(p,t,k,n);
		if(gp == NULL)
		{
			fprintf(stderr, "failed allocation\n");
			return 1;
		}
		qcc = new_quad_cyclic_cover(gp,s,d, NB_VECTORS);
		if(qcc == NULL)
		{
			fprintf(stderr, "failed allocation\n");
			return 1;
		}

		init_GS(NB_VECTORS);
		set_random_vectors(qcc);
		check_orthogonality(qcc);

		free_generalized_permutation(&gp);
		free_quad_cyclic_cover(&qcc);
		free_GS();
	}
	fprintf(stdout,"done\n\n");
	fflush(stdout);

	return 0;
}

int test_qcc_H_plus(int *p, int *t, int *s, int d, int k, int n)
{
	size_t i,j,nb_ren=0;
	double * theta;

	gp = new_generalized_permutation(p,t,k,n);
	if(check_generalized_permutation(gp))
	{
		fprintf(stderr,"Error with generalized permutation\n");
		return 1;
	}
	qcc = new_quad_cyclic_cover(gp, s, d, NB_VECTORS);
	if(qcc == NULL)
	{
		fprintf(stderr, "Allocation error with quad. cylc. cov.\n");
		return 1;
	}
	set_random_lengths_quad_cyclic_cover(qcc);
	if(check_quad_cyclic_cover(qcc))
	{
		fprintf(stderr, "Error checking qcc.\n");
		return 1;
	}
	print_quad_cyclic_cover(qcc);

	init_GS(qcc->degree * qcc->nb_labels);

	theta = (double *) malloc((qcc->nb_vectors + 1) * sizeof(double));
	if(theta == NULL)
	{
		fprintf(stderr, "Allocation error for vectors\n");
		return 1;
	}

	for(i=0; i<NB_VECTORS+1; ++i) theta[i] = 0.;
	set_random_vectors(qcc);

	for(i=0; i<5000; ++i)
	{
		rauzy_induction_H_plus_quad_cyclic_cover(qcc);

		if(qcc->length < 0.00001)
		{
			theta[0] -= logl(qcc->length);
			fprintf(stdout,"renormalize at step %4u, length=%.12f\n",(unsigned int) i, (double) qcc->length);

			orthogonalize_GS(qcc, theta);

			printf("RZ Lexp: %f\n",theta[0] / i);
			printf("theta:");
			for(j=1; j<NB_VECTORS+1; ++j) printf(" %f", theta[j] / (2*theta[0]));
			printf("\n");

			
			if((nb_ren+1) % 100 == 0)  // a bit of salt
			{
				qcc->top->lab->length += ((double) (.5-drand())) / 0x40000000000;
				qcc->bot->lab->length += ((double) (.5-drand())) / 0x40000000000;
			}
			renormalize_length_quad_cyclic_cover(qcc);

			nb_ren += 1;
		}

		if(check_quad_cyclic_cover(qcc))
		{
			fprintf(stderr, "Error checking qcc at iteration step %d\n",(int)i);
			return 1;
		}

	}


	print_quad_cyclic_cover(qcc);

	fprintf(stdout, "free gp...");
	fflush(stdout);
	free_generalized_permutation(&gp);
	fprintf(stdout, "done\n");

	fprintf(stdout, "free qcc...");
	fflush(stdout);
	free_quad_cyclic_cover(&qcc);
	fprintf(stdout, "done\n");

	fprintf(stdout, "theta at %p\n",theta);
	fflush(stdout);
	free(theta);
	free_GS();

	return 0;
}

int main(void)
{
	size_t i;

//	for(i=0; i < 10; ++i) test_vectors();

	init_examples();

//	for(i=0; i<NB_TEST; i++)
//		if(test_allocation(P[i],T[i],S[i],D[i],K[i],N[i]))
//		{
//			fprintf(stderr,"error in allocation test %lu\n",i);
//			exit(EXIT_FAILURE);
//		}

	for(i=0; i < NB_TEST; ++i)
	{
		if(test_qcc_H_plus(P[i],T[i],S[i],D[i],K[i],N[i]))
		{
			fprintf(stderr,"error in test 0\n");
			exit(EXIT_FAILURE);
		}
	}

	return 0;
}
