/* Some highly non optimized linear algebra routines */
/* 1) random orthogonal vectors                      */
/* 2) orthogonalization                              */

#include "lyapunov_exponents.h"

/* 0 <= i < nb_vectors   */
/* 0 <= j < degree       */
/* 0 <= k < nb_intervals */
/*  index at pos (i,j,k) = v[k + nb_intervals * (j + degree * i)] */

/* global variable used in GS orthogonalisation                                 */
/* they must be dynamic so we keep them here in order to avoid many malloc/free */
double * scal     = NULL;
double * scal_new = NULL;
size_t scal_size  = 0;

void print_vectors(quad_cyclic_cover *qcc)
{
	size_t i,j,k;

	for(i=0; i<qcc->nb_vectors; ++i)
	{
		printf("v[%lu]=",i);
		for(j=0; j < qcc->nb_labels; ++j)
		{
			for(k=0; k < qcc->degree; ++k)
				printf(" %f", (qcc->labels)[j].v[i + qcc->nb_vectors * k]);
			printf(" |");
		}
		printf("\n");
	}
}

void set_random_vectors(quad_cyclic_cover *qcc)
/* set random orthogonal frame */
/* warning: before calling this function call init_GS(dim) */
{
	size_t i,j;
	double *vv,norm;

	for(j=0; j<qcc->nb_labels; ++j)
	{
		vv = (qcc->labels)[j].v;
		for(i=0; i<qcc->nb_vectors*qcc->degree; i++)
			vv[i] = drand() - .5;
	}
	if(qcc->nb_vectors == 1)
	{
		for(i=0; i<qcc->nb_labels; i++)
			for(j=0; j<qcc->degree; j++)
			{
				if((qcc->labels)[i].v[j] > 0) norm += (qcc->labels)[i].v[j];
				else norm -= (qcc->labels)[i].v[j];
			}
	}
	else orthogonalize_GS(qcc, NULL);
}

int init_GS(size_t nb_vectors)
/* allocate scal and scal_new in order that it contains at least 2*nb_vectors elements */
{
	if(scal == NULL)
	{
		scal_size = 256;
		/* warning: scal and scal_new are switched in the GS orthogonalization */
		/* so it is safer to allocate them independently                       */
		scal = (double *) malloc(scal_size * sizeof(double));
		scal_new = (double *) malloc(scal_size * sizeof(double));
		//printf("info: GS allocation with size %lu at %p and %p\n",scal_size, scal, scal_new);
	}

	else if(scal_size < nb_vectors)
	{
		while(scal_size < nb_vectors) scal_size = (scal_size * 3) / 2;
		scal = (double *) realloc(scal, scal_size * sizeof(double));
		scal_new = (double *) realloc(scal_new, scal_size * sizeof(double));
		//printf("info: GS reallocation with size %lu at %p and %p\n",scal_size,scal,scal_new);
	}

	if(scal == NULL)
	{
		scal_size = 0;
		return 1;
	}

	return 0;
}

void free_GS(void)
{
	scal_size = 0;
	if(scal != NULL)
	{
		free(scal);
		free(scal_new);
		scal = NULL;
		scal_new = NULL;
	}
}

void orthogonalize_GS(quad_cyclic_cover *qcc, double *theta)
/* orthogonalization using Gram-Schmidt process                                                */
/* warning: this method should be called AFTER being sure that init_GS has been called         */
/* it theta is not NULL it is updated by 2*log of the diagonal from position 1 to nb_vectors   */
/* warning: it is 2*log(entry) and not log(entry)                                              */
/* element at position (i,j,k): (qcc->labels)[j].v[i + nb_vectors*k];                          */
{
	size_t i,ii,j,k;

	double norm,sqnorm,c;
	double *tmp=NULL;
	double *vv;

	/* some check that we will remove */
	if(qcc->nb_vectors < 2)
	{
		fprintf(stderr, "calling GS with nb_vectors < 2 is not possible.\n");
		exit(EXIT_FAILURE);
	}
	if(qcc->nb_vectors > qcc->degree * qcc->nb_labels)
	{
		fprintf(stderr, "too much vectors for the dimension!\n");
		exit(EXIT_FAILURE);
	}
	if(scal == NULL || scal_new == NULL)
	{
		fprintf(stderr, "you must call init_GS before calling orthogonalize_GS.\n");
		exit(EXIT_FAILURE);
	}

	/* put v in the holonomy free subspace */
	/* compute <v_0,v_1> in scal_new[0] */
	/* compute <v_0,v_0> in sqnorm */
	scal_new[0] = 0.;
	sqnorm = 0.;
	for(j=0; j < qcc->nb_labels; ++j)
	{
		vv = (qcc->labels)[j].v;
		for(k=0; k<qcc->degree; ++k)
		{
			scal_new[0] += vv[0 + qcc->nb_vectors*k] * vv[1 + qcc->nb_vectors*k];
			sqnorm  += vv[0 + qcc->nb_vectors*k] * vv[0 + qcc->nb_vectors*k];
		}
	}


	/* vector by vector orhtonormalisation */
	for(i=1; i < qcc->nb_vectors-1; ++i)
	{
		/* sqnorm contains <v_(i-1),v_(i-1)> */
		/* v_0, v_1, ..., v_(i-2) are normalized */
		/* v_0, v_1, ..., v_(i-1) are orthogonal */
		/* scal_new contains the i scalar products <v_0,v_i>, <v_1,v_i>, ..., <v_(i-1),v_i> */
		tmp = scal; scal = scal_new; scal_new = tmp;
		for(ii=0; ii<=i; ++ii)
			scal_new[ii] = 0.;

		c = scal[i-1]/sqnorm;
		if(theta != NULL) theta[i] += log(sqnorm);
		norm = sqrt(sqnorm);
		sqnorm = 0.;

		/* c = <v_i,v_(i-1)> / <v_(i-1,v_(i-1)> */
		for(j=0; j < qcc->nb_labels; ++j)
		{
			vv = (qcc->labels)[j].v;
			for(k=0; k < qcc->degree; ++k)
			{
				/* subtract the part of the span of v_0, v_1, ..., v_(i-2) */
				for(ii=0; ii<i-1; ++ii)
					vv[i + qcc->nb_vectors * k] -= scal[ii] * vv[ii + qcc->nb_vectors * k];

				/* subtract the part of the span of v_(i-1) */
				vv[i + qcc->nb_vectors * k] -= c * vv[ii + qcc->nb_vectors*k];

				/* normalize v_(i-1) */
				vv[(i-1) + qcc->nb_vectors*k] /= norm;

				/* compute scalar products and norms for next loop */
				/* sqnorm = <v_i, v_i> */
				/* scal_new[ii] = <v_(i+1), v_ii>*/
				sqnorm += vv[i + qcc->nb_vectors * k] * vv[i + qcc->nb_vectors *k];
				for(ii=0; ii<=i; ++ii)
					scal_new[ii] += vv[ii + qcc->nb_vectors*k] * vv[i+1 + qcc->nb_vectors*k];
			}
		}
	}

	/* here i = NB_VECTORS-1 */
	/* renormalize v_(i-1)   */
	/* put v_i in the orthogonal of the span of v_0, v_1, ..., v_(i-1) */
	c = scal_new[i-1] / sqnorm;
	if(theta != NULL) theta[i] += log(sqnorm);
	norm = sqrt(sqnorm);
	sqnorm = .0;

	for(j=0; j< qcc->nb_labels; ++j)
	{
		vv = (qcc->labels)[j].v;
		for(k=0; k<qcc->degree; ++k)
		{
			for(ii=0; ii<i-1; ++ii)
				vv[i + qcc->nb_vectors*k] -= scal_new[ii] * vv[ii + qcc->nb_vectors*k];

			vv[i + qcc->nb_vectors*k] -= c * vv[i-1 + qcc->nb_vectors*k];

			vv[i-1 + qcc->nb_vectors * k] /= norm;

			sqnorm += vv[i + qcc->nb_vectors * k] * vv[i + qcc->nb_vectors * k];
		}
	}

	/* we renormalize v_i*/
	if(theta != NULL) theta[i+1] += log(sqnorm);
	norm = sqrt(sqnorm);
	for(j=0; j < qcc->nb_labels; ++j)
	{
		vv = (qcc->labels)[j].v;
		for(k=0; k < qcc->degree; ++k)
			vv[i + qcc->nb_vectors * k] /= norm;
	}
}

void check_orthogonality(quad_cyclic_cover *qcc)
{
	double s;
	double *vv;
	size_t i1,i2,j,k;

	for(i1=0; i1<qcc->nb_vectors; ++i1)
	{
		for(i2=i1; i2<qcc->nb_vectors; ++i2)
		{
			s = 0;
			for(j = 0; j < qcc->nb_labels; ++j)
			{
				vv = (qcc->labels)[j].v;
				for(k=0; k < qcc->degree; ++k)
					s += vv[i1 + qcc->nb_vectors*k] * vv[i2 + qcc->nb_vectors*k];
			}
			printf("<v%d,v%d> = %f\t",(int) i1, (int) i2, s);
		}
		printf("\n");
	}
}
