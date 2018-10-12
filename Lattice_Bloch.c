/****************************************************************************
*                                                                          *
* File    : Lattice_Bloch.c                                                *
*                                                                          *
* Purpose : Calculates the 2D Bloch Green function  from the 0th cell      *
*                                                                          *
* History : Date      Reason                                               *
*           01/02/17  Created by Rozalia Lukacs                            *
*           30/06/2017 Single precission is set, wavenumber cycle removed  *
*           02/10/2017 Single precisiion removed, double precission set;   *
*                      wavenumber cycle addded                             *
*                                                                          *
****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <mathimf.h>
#include <mkl_types.h>
#include <mkl_lapack.h>
#include <mkl.h>

// ZGESV prototype - solves the linear system of equation with double complex numbers

lapack_int LAPACKE_zgesv(int matrix_layout, lapack_int n, lapack_int nrhs, lapack_complex_double * A, lapack_int lda, lapack_int * ipiv, lapack_complex_double * B, int ldb);

//Define the constants

#define CC 3e8
#define EM_const 0.57721566490153286060651209008240243104215933593992 
#define M_PI 3.14159265358979323846
#define nano 0.000000001
#define Nx 100  //resolution in the x direction
#define Ny 100  //resolution in the y dimension
#define N Nx*Ny // total resolution of the Green function
#define NRHS 1  // parameter of the linear equation solution function
#define LDA N   // parameter of the ---||---
#define LDB N   // parameter of the ---||---
#define LAPACK_ROW_MAJOR 101
#define LAPACK_COL_MAJOR 102

//Sign function for int numbers

int signnum_typical(int x) {
  if (x > 0.0) return 1;
  if (x < 0.0) return -1;
  return 0;
}

/****************************************************************************
*                                                                          *
* Function: main                                                           *
*                                                                          *
* Purpose : Calculates the wave function of a 2D system with Bloch 	    *
*            boundary condition in x direction.                            *
*                                                                          *
* History : Date      Reason                                               *
*           01/02/17  Created                                              *
*                                                                          *
****************************************************************************/

int main(int argc, char *argv[])
{
	int funcGreen, MaxNrCell, source, ArrayOfCircles, Zebra, NCx, NCy;
	int i, j, l, m, MM, nx, ny, en, en1, en2, enmax;
	float* p = NULL;
	float* q = NULL;
	float lambda, lambdamin, lambdamax;
	float phi, k, kx, ky, K;
	 R, F, T, a, b,  Z;
	float dx, dy, x_null, y_null, x[N], y[N];
	float xcircle, ycircle, rn;
	float dummy3, dummy4, dummy1, dummy2;
	double _Complex DiagE;
	float _Complex n_index, n_index_film, Gp, dummy5, dummy6; 
	double _Complex v[N], Gtilde[N][N], G[N][N], eikr[N], M[N][N], psi_complex[N*NRHS];
	lapack_complex_double BT[N*NRHS], MT[N*N];
	lapack_int n = N, nrhs = NRHS, lda = LDA, ldb = LDB, info, IPIV[N];

	FILE *f;
	
	// read in the input parameters from the stdio -- I do this by adding the input file on the stdio

	scanf("%d\n",&funcGreen);  // the type of the Green function
	scanf("%d\n",&MaxNrCell);  // the maximum number of simulation cells, this have to be odd number
	scanf("%d\n",&source);     // the type of the source - ex. plane wave
	scanf("%f\n",&lambdamin);  // the left value of the lambda interval
	scanf("%f\n",&lambdamax);  // the right value of the lambda interval
	scanf("%d\n",&ArrayOfCircles); // do we have cyrcles yes=1, no=0
	scanf("%d\n",&Zebra);          // do we have a film yes=1; no=0
	scanf("%d\n",&NCx);            // number of disks in x direction
	scanf("%d\n",&NCy);            // number of disks in y direction
	scanf("%lf\n",&phi);           // the angle of the k vector of the plane wave -- incoming direction
	scanf("%f\n",&R);              // the radius of the disks
	scanf("%f\n",&T);			   // the thickness of the film
	scanf("%f\n",&dummy1); // real part of refractive index of disk
	scanf("%f\n",&dummy2); // imaginary part of refractive index of disk
	scanf("%f\n",&dummy3); // real part of refractive index of film
	scanf("%f\n",&dummy4); // imaginary part of refractive index of film
	
	printf("Started the program %d \n",N);
	
	// Define the geometry of the system
	
	a=0;   //initialize the variables
	b=0;
	
	if (ArrayOfCircles == 1)
		{
			F=0.0*R;
			a=2*NCx*R;
			b=2*NCy*R+2*F;
			//printf("The a	%f	b	%f\n", a,b);
		}

	if (Zebra == 1)
	{
		//a=2*T;
		b=b+(1.0*T);	
		//printf("The a	%12.12lf	b	%12.12lf	T	%12.12lf\n", a,b,T);
	}

	// construct the refractive indexes of the film and the disks
	
	n_index=dummy1+I*dummy2;
	n_index_film=dummy3+I*dummy4;	
	//printf("The n_index	%f	%f	n_index_film	%f	%f\n", creal(n_index),cimag(n_index),creal(n_index_film),cimag(n_index_film));

	// define the step in the space
	
	dx=a/Nx;
	dy=b/Ny;
	//printf("The dx	%12.12lf	dy	%12.12f\n", dx,dy);
	
	// define the origo of the simulation cell
	x_null=0;
	//y_null=0;
	//x_null=0;
	y_null=1.0*T;
	//printf("The x_null	%f	y_null	%12.12lf\n", x_null, y_null);	

	//define the grid in cell # 0

	for (i = 0; i < Nx; i++)
	{
		for (l = 0; l < Ny; l++)
		{
			j=l*Nx+i;
			x[j]=(i+0.5)*dx;
			y[j]=(l+0.5)*dy;
			//printf("The j	%d	x	%12.12f	y	%12.12f\n", j, x[j], y[j]);
		}
	}

	//Establish the potential

	for (m=0; m < N; m++)
	{
		if (ArrayOfCircles == 1)
		{		
			for (nx=1; nx <= NCx; nx++)
				for (ny=1; ny <= NCy; ny++)
				{
					xcircle=(2*nx-1)*R+x_null;
					ycircle=(2*ny-1)*R+y_null;

					rn=sqrt((x[m]-xcircle)*(x[m]-xcircle)+(y[m]-ycircle)*(y[m]-ycircle));

					if (rn < R)
						v[m]=1.0-n_index*n_index;
					else
						v[m]=0.0;

				}
		}
		if (Zebra == 1)
		{
			if ((y[m] > 0.0*T)&&(y[m] < 1.0*T))
				{
				v[m]=1.0-n_index_film*n_index_film;
				}
			else
				if (ArrayOfCircles == 0)
					v[m]=0.0;

		}
	//printf("The real part of V	%12.12f	imaginary part	%12.12f\n",creal(v[m]), cimag(v[m]));
	}

    // start the wavenumber cycle

	for (lambda = lambdamin; lambda <= lambdamax; lambda = lambda + nano)
	{
		
		// Define the energy
		
		k=((2*M_PI)/lambda);
		kx=k*cos(phi);
		ky=k*sin(phi);
		K=sqrt((kx*kx)+(ky*ky));
		//printf("The kx	%12.12lf ky	%12.12lf  k	%12.12lf	and K	%12.12lf\n",kx,ky,k,K);

		// Define the Bloch momentum p_n and q_n This part is not done!!!! Now I use p0=kx

		enmax=1;

		p=(float*) realloc (p, enmax * sizeof(float));
		q=(float*) realloc (q, enmax * sizeof(float));

		p[0]=kx;
		q[0]=sqrt(K*K-p[0]*p[0]);
	/*
		en1=floor((-a/(2*M_PI))*(kx+K));
		en2=floor((-a/(2*M_PI))*(kx-K)); // number of open channels
		
		printf("The number of open chanels en1	%d	en2	%d\n", en1, en2);

		if (signnum_typical(en1)==signnum_typical(en2))
			
			if (abs(en1) <= abs(en2))
					enmax=abs(en2)-abs(en1);
			else
					enmax=abs(en1)-abs(en2);
		else
			enmax=abs(en1)+abs(en2);

		printf("The maximum number of open chanels enmax	%d\n", enmax);

		p=(float*) realloc (p, enmax * sizeof(float));
		q=(float*) realloc (q, enmax * sizeof(float));

		for (en=0; en < enmax; en++)
			{
				p[en]=kx+(en*2*M_PI)/a;
				printf("The p[en]	%f",p[en]);
				q[en]=sqrt(K*K-p[en]*p[en]);
				//printf("The en	%d p[en]	%f	q[en]	%f\n",en,p[en],q[en]);

			}
			
		// add one closed chanels

		//q[en2+1]=kx+en*(2*M_PI)/a;
		//p[en2+1]=sqrt(K*K-p[en]*p[en]);
	*/
		//Calculate the diagonal elements of the Green function
		
		DiagE=((I/4)-(1/(2*M_PI))*(log(k*dx/2)+EM_const)-(0.25*log(0.5)-0.75+(M_PI/8))/M_PI)*(dx*dy);
		//printf("The diagonal elements of the Green function %12.18lf	i*%12.18lf\n",creal(DiagE),cimag(DiagE));
		
		//Testing output for testing purposes
		
		f = fopen("Output_for_testing_100x100_NRcell_51.txt","a");		
		
		fprintf(f,"%12.12lf\n",lambda);
		fprintf(f,"%f\n",creal(n_index));
		fprintf(f,"%12.12lf\n",dx);
		fprintf(f,"%12.12lf\n",dy);
		fprintf(f,"%12.12lf\n",a);
		fprintf(f,"%12.12lf\n",b);


		for(j=0; j < N; j++)
			fprintf(f,"%d	%16.16lf	%16.16lf\n", j, x[j], y[j]);
		fclose(f);

		// calculate the Bloch Green function

		for (MM=-((MaxNrCell-1)/2); MM <= ((MaxNrCell-1)/2); MM++)
		{
			for (j=0; j < N; j++ )
			{
			 //printf("The jth row  %d	", j);
				for (m=0; m < N; m++)
				{

					if ((j == m)&&(MM==0))
						G[j][m]=DiagE;
					else
					{
						Z=k*sqrt((x[j]-((MM*a)+x[m]))*(x[j]-((MM*a)+x[m]))+(y[j]-y[m])*(y[j]-y[m]));
						G[j][m]=j0(Z)+I*y0(Z);
					}
					Gtilde[j][m]=0.0+I*0.0;
					//printf("%2.16f  +  %2.16f i	",creal(G[j][m]), cimag(G[j][m]));
					
					/*if there are several cells the Bloch momentum can not take any value.
					The energy is redistributed among open chanels. The contribution from all
					must be included.*/

					if ( MaxNrCell == 1 )
						Gtilde[j][m]=Gtilde[j][m]+cexp(I*MM*(p[0]*a))*G[j][m];
					else
						{
						Gp=cexp(I*MM*(p[0]*a))*G[j][m];
						Gtilde[j][m]=Gtilde[j][m]+Gp;
						}
						// This branch is not finished yet
						/*for (en=0; en < enmax; en++)
						{
							Gp=cexp(I*MM*(p[en]*a))*G[j][m];
							printf("The en %d	p[en]	%f",en,p[en]);
							//dummy1=creal(Gp)*creal(G[j][m])-cimag(Gp)*cimag(G[j][m])+I*(creal(Gp)*cimag(G[j][m])+cimag(Gp)*creal(G[j][m]));
							//dummy2=creal(dummy1)+creal(Gtilde[j][m])+I*(cimag(dummy1)+cimag(Gtilde[j][m]));
							Gtilde[j][m]=Gtilde[j][m]+Gp;		
						}*/
				//printf("%10.16f  +  %10.16f i	",creal(Gtilde[j][m]), cimag(Gtilde[j][m]));
				//printf("\n");
				}
			}
		}

		for (j=0; j< N; j++)
		{
			//printf("The j row %d	", j);
			for (m=0; m < N; m++)
			{
				//printf("%10.10f  +  %10.10f i	",creal(Gtilde[j][m]), cimag(Gtilde[j][m]));
				Gtilde[j][m]=G[j][m]*v[m];
				//printf("%10.16f  +  %10.16f i	\n",creal(Gtilde[j][m]), cimag(Gtilde[j][m]));	
			}
		//printf("\n");
		}
		
		// create the source -- now a plane wave

		for (j=0; j < N; j++)
		{
			eikr[j]=cexp(I*((kx*x[j])+(ky*y[j])));
			//printf("The jth	%d eikr[j].real	%12.12f	eikr[j].imag	%12.12f\n",j,creal(eikr[j]), cimag(eikr[j]));
		}

		// initialize matrices for the Lippman-Swinger equation

		for (j=0; j< N; j++)
		{
		//printf("The j sor %d	", j);
			for (m=0; m < N; m++)
			{
				//dummy1=-(k*k*dx*dy*0.25)*cimag(Gtilde[j][m])+I*(k*k*dx*dy*0.25)*creal(Gtilde[j][m]);
				
				if (j==m)
					M[j][m]=1.0+I*(k*k*dx*dy*0.25)*(Gtilde[j][m]);
					//M[j][m]=1.0+creal(dummy1)+I*cimag(dummy1);
				else
					M[j][m]=0.0+I*(k*k*dx*dy*0.25)*(Gtilde[j][m]);
					//M[j][m]=0.0+creal(dummy1)+I*cimag(dummy1);
				
				//printf("%10.16f  +  %10.16f i	",creal(M[j][m]), cimag(M[j][m]));
			}
		//printf("\n");
		}

		//prepare your matrix for a FORTRAN routine

		for (j = 0; j < N; j++) /* to call a Fortran routine from C we */
		{       		       /* have to transform the matrix */
			for (m = 0; m < N; m++)
			{
			MT[m + N * j].real = creal(M[j][m]);
			MT[m + N * j].imag = cimag(M[j][m]);
			//printf("The j, m	%d	%d	MT.real	%12.16f	MT.imag	%12.16f\n", m , j, MT[j + N * m].real, MT[j + N * m].imag);
			}

			BT[j].real = creal(eikr[j]);
			BT[j].imag = cimag(eikr[j]);
			//printf("Az egyenlet elott j	%d	BT.real	%12.12f	BT.imag	%12.12f\n", j, BT[j].real, BT[j].imag);
		}
		
		info=LAPACKE_zgesv( LAPACK_COL_MAJOR, n, nrhs, MT, lda, IPIV, BT, ldb ); // this is the LAPACK routine, which solves the system of linear equations with complex double values
		
		printf("Output of ZGESV = %d\n",info);
		
		if (info ==0)
		   {
		   for (j=0; j < N*NRHS; j++)
			{
				psi_complex[j]=BT[j].real+I*BT[j].imag;
				//printf("The j	%d	Psi_complex	%12.12f	+ %12.12f i\n", j, BT[j].real, BT[j].imag);
			}
		   }
		else
		   {
		   for (j=0; j < N*NRHS; j++)
				psi_complex[j]=0.0;
		   }

		// write the calculated wavefunction in a file

		f = fopen("psi_complex_100x100_NRcell_51_interval.txt","a");

		fprintf(f,"%lf\n",lambda);

		for(j=0; j < N*NRHS; j++)
			fprintf(f,"%16.16lf	%16.16lf\n",creal(psi_complex[j]),cimag(psi_complex[j]));

		fclose(f);

		//write out values of x, y,v, eikr, etc. for testing purposes

		f = fopen("Output_for_testing_single_100x100.txt","a");		
		
		fprintf(f,"%12.12lf\n",lambda);
		fprintf(f,"%f\n",creal(n_index));
		fprintf(f,"%12.12lf\n",dx);
		fprintf(f,"%12.12lf\n",dy);
		fprintf(f,"%12.12lf\n",a);
		fprintf(f,"%12.12lf\n",b);

		for(j=0; j < N; j++)
			fprintf(f,"%16.16lf	%16.16lf	%16.16lf	%16.16lf	%16.16lf	%16.16lf\n", creal(v[j]), cimag(v[j]), x[j], y[j], creal(eikr[j]), cimag(eikr[j]));
		fclose(f);

	}
	printf("Hello, world!\n");
	return 0;
}

