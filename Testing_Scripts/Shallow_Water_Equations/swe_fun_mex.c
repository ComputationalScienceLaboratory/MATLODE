#include "mex.h"
#include <stdio.h>
#define g 9.81
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
		int M,N,i,j,len,l;	
		double *A, dx,dy;
		double *outArray;
		len=mxGetM(prhs[0]);
		A=mxGetPr(prhs[0]);
		M=mxGetScalar(prhs[1]);
		N=mxGetScalar(prhs[2]);
		dx=mxGetScalar(prhs[3]);
		dy=mxGetScalar(prhs[4]);
		double x2d[M][N], u2d[M][N], v2d[M][N], Hx[M-1][N-1],Ux[M-1][N-1],Vx[M-1][N-1];
		double Hy[M-1][N-1], Uy[M-1][N-1], Vy[M-1][N-1];
		l=M*N;
		for(i=0;i<M-1;i++)
		{
				for(j=0;j<N-1;j++)
				{
						Hx[i][j]=0;
						Hy[i][j]=0;
						Ux[i][j]=0;
						Uy[i][j]=0;
						Vx[i][j]=0;
						Vy[i][j]=0;
				}
		}
		for (i=0;i<M;i++)
		{
				for(j=0;j<N;j++)
				{
						int p =i*M+j;
						x2d[i][j]=A[p];
						u2d[i][j]=A[l+p];
						v2d[i][j]=A[2*l+p];
				}
		}
		for (i=0;i<M-1;i++)
		{
				for(j=0;j<N-2;j++)
				{
						Hx[i][j] = (x2d[i+1][j+1]+x2d[i][j+1])/2;
						Ux[i][j] = (u2d[i+1][j+1]+u2d[i][j+1])/2;
						Vx[i][j] = (v2d[i+1][j+1]+v2d[i][j+1])/2;
				}
		}
		for (i=0;i<M-2;i++)
		{
				for(j=0;j<N-1;j++)
				{
						Hy[i][j] = (x2d[i+1][j+1]+x2d[i+1][j])/2;
						Uy[i][j] = (u2d[i+1][j+1]+u2d[i+1][j])/2;
						Vy[i][j] = (v2d[i+1][j+1]+v2d[i+1][j])/2;
				}
		}
		
		
		for (i=1;i<M-1;i++)
		{
				for(j=1;j<N-1;j++)
				{
						x2d[i][j] = -(1/dx)*(Ux[i][j-1]-Ux[i-1][j-1]) - 
								    (1/dy)*(Vy[i-1][j]-Vy[i-1][j-1]);
						u2d[i][j] = - (1/dx)*((Ux[i][j-1]*Ux[i][j-1]/Hx[i][j-1] + g/2*Hx[i][j-1]*Hx[i][j-1]) -
										    (Ux[i-1][j-1]*Ux[i-1][j-1]/Hx[i-1][j-1] + g/2*Hx[i-1][j-1]*Hx[i-1][j-1]))-
								    (1/dy)*((Vy[i-1][j]*Uy[i-1][j]/Hy[i-1][j]) - 
													    (Vy[i-1][j-1]*Uy[i-1][j-1]/Hy[i-1][j-1]));
						v2d[i][j] = - (1/dx)*((Ux[i][j-1]*Vx[i][j-1]/Hx[i][j-1]) - 
										    (Ux[i-1][j-1]*Vx[i-1][j-1]/Hx[i-1][j-1])) 
								    - (1/dy)*((Vy[i-1][j]*Vy[i-1][j]/Hy[i-1][j] + g/2*Hy[i-1][j]*Hy[i-1][j]) - 
										    (Vy[i-1][j-1]*Vy[i-1][j-1]/Hy[i-1][j-1] + g/2*Hy[i-1][j-1]*Hy[i-1][j-1]));
				}
		}
		for(i=0;i<M;i++)
		{
				x2d[i][0]=x2d[i][1];
				u2d[i][0]=u2d[i][1];
				v2d[i][0]=-v2d[i][1];
				x2d[i][N-1]=x2d[i][N-2];
				u2d[i][N-1]=u2d[i][N-2];
				v2d[i][N-1]=-v2d[i][N-2];
		}
		for(i=0;i<N;i++)
		{
				x2d[0][i]=x2d[1][i];
				u2d[0][i]=-u2d[1][i];
				v2d[0][i]=v2d[1][i];
				x2d[M-1][i]=x2d[M-2][i];
				u2d[M-1][i]=-u2d[M-2][i];
				v2d[M-1][i]=v2d[M-2][i];
		}
		plhs[0]=mxCreateDoubleMatrix(len,1,mxREAL);
		outArray=mxGetPr(plhs[0]);
		for(i=0;i<M;i++)
		{
				for(j=0;j<N;j++)
				{
						int p = i*M+j;
						outArray[p]=x2d[i][j];
						outArray[l+p]=u2d[i][j];
						outArray[2*l+p]=v2d[i][j];
				}
		}

		return;
		
		
}




					


