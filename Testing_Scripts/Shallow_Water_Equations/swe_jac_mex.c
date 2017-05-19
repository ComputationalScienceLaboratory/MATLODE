#include <stdio.h>
#include "mex.h"
#include "matrix.h"
#define LEN 40
#define NUMNONZEROS 50
#define g 9.81
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
		int M,N, i,j,len,l,sparseLength, tmp,m,boundary_len, tmp1,n;
		double *x,dx,dy,*VAL,a1,a2,a3,a4,a5;
		double *ROW, *COL;
		double a6,a7,a8,a9,a10, a11,a12;
		x=mxGetPr(prhs[0]);
		M=mxGetScalar(prhs[1]);
		N=mxGetScalar(prhs[2]);
		dx=mxGetScalar(prhs[3]);
		dy=mxGetScalar(prhs[4]);
		l=M*N;
		len = mxGetM(prhs[0]);
		sparseLength=l*LEN;
		boundary_len=(2*(M+N)-4)*NUMNONZEROS;
		unsigned int col[sparseLength], row[sparseLength],row1[boundary_len],col1[boundary_len];
		double val[sparseLength],val1[boundary_len];
		tmp=-1;
		tmp1=-1;
		m=M;
		n=N;
		/*printf ("%d %d\n",M,N);*/
		for(i=0;i<M-2;i++)
		{
				for(j=0;j<N-2;j++)
				{
						int k = (i+1)*M+j+1;
						tmp=tmp+1;
						row[tmp]=k;
						col[tmp]=k+l+M;
						val[tmp]=-1/(2*dx);
						tmp=tmp+1;
						row[tmp]=k;
						col[tmp]=k+l-m;
						val[tmp]=1/(2*dx);
						tmp = tmp+1;
						row[tmp]=k;
						col[tmp]=k+2*l+1;
						val[tmp]=-1/(2*dy);
						tmp=tmp+1;
						row[tmp]=k;
						col[tmp]=k+2*l-1;
						val[tmp]=1/(2*dy);
						tmp=tmp+1;
						if(i==0 || j==0)
						{
								if(i ==0 && j==0)
								{
										tmp1=tmp1+1;
										row1[tmp1]=k-n-1;
										col1[tmp1]=k+l+m;
										val1[tmp1]=-1/(2*dx);
										tmp1=tmp1+1;
										row1[tmp1]=k-n-1;
										col1[tmp1]=k+l-m;
										val1[tmp1]=1/(2*dx);
										tmp1=tmp1+1;
										row1[tmp1]=k-n-1;
										col1[tmp1]=k+2*l+1;
										val1[tmp1]=-1/(2*dy);
										tmp1=tmp1+1;
										row1[tmp1]=k-n-1;
										col1[tmp1]=k+2*l-1;
										val1[tmp1]=1/(2*dy);
								}
						
								if(i==0)
								{
										tmp1=tmp1+1;
										row1[tmp1]=k-n;
										col1[tmp1]=k+l+m;
										val1[tmp1]=-1/(2*dx);
										tmp1=tmp1+1;
										row1[tmp1]=k-n;
										col1[tmp1]=k+l-m;
										val1[tmp1]=1/(2*dx);
										tmp1=tmp1+1;
										row1[tmp1]=k-n;
										col1[tmp1]=k+2*l+1;
										val1[tmp1]=-1/(2*dy);
										tmp1=tmp1+1;
										row1[tmp1]=k-n;
										col1[tmp1]=k+2*l-1;
										val1[tmp1]=1/(2*dy);
								}
								if(j==0)
								{
										tmp1=tmp1+1;
										row1[tmp1]=k-1;
										col1[tmp1]=k+l+m;
										val1[tmp1]=-1/(2*dx);
										tmp1=tmp1+1;
										row1[tmp1]=k-1;
										col1[tmp1]=k+l-m;
										val1[tmp1]=1/(2*dx);
										tmp1=tmp1+1;
										row1[tmp1]=k-1;
										col1[tmp1]=k+2*l+1;
										val1[tmp1]=-1/(2*dy);
										tmp1=tmp1+1;
										row1[tmp1]=k-1;
										col1[tmp1]=k+2*l-1;
										val1[tmp1]=1/(2*dy);
								}
								if(i==0 && j==N-3)
								{
										tmp1=tmp1+1;
										row1[tmp1]=k-n+1;
										col1[tmp1]=k+l+m;
										val1[tmp1]=-1/(2*dx);
										tmp1=tmp1+1;
										row1[tmp1]=k-n+1;
										col1[tmp1]=k+l-m;
										val1[tmp1]=1/(2*dx);
										tmp1=tmp1+1;
										row1[tmp1]=k-n+1;
										col1[tmp1]=k+2*l+1;
										val1[tmp1]=-1/(2*dy);
										tmp1=tmp1+1;
										row1[tmp1]=k-n+1;
										col1[tmp1]=k+2*l-1;
										val1[tmp1]=1/(2*dy);
								}
								if(i==M-3 && j==0)
								{
										tmp1=tmp1+1;
										row1[tmp1]=k+n-1;
										col1[tmp1]=k+l+m;
										val1[tmp1]=-1/(2*dx);
										tmp1=tmp1+1;
										row1[tmp1]=k+n-1;
										col1[tmp1]=k+l-m;
										val1[tmp1]=1/(2*dx);
										tmp1=tmp1+1;
										row1[tmp1]=k+n-1;
										col1[tmp1]=k+2*l+1;
										val1[tmp1]=-1/(2*dy);
										tmp1=tmp1+1;
										row1[tmp1]=k+n-1;
										col1[tmp1]=k+2*l-1;
										val1[tmp1]=1/(2*dy);

								}
						}
						if(i==M-3 ||j==N-3)
						{
								if(i==M-3 && j== N-3)
								{
										tmp1=tmp1+1;
										row1[tmp1]=k+n+1;
										col1[tmp1]=k+l+m;
										val1[tmp1]=-1/(2*dx);
										tmp1=tmp1+1;
										row1[tmp1]=k+n+1;
										col1[tmp1]=k+l-m;
										val1[tmp1]=1/(2*dx);
										tmp1=tmp1+1;
										row1[tmp1]=k+n+1;
										col1[tmp1]=k+2*l+1;
										val1[tmp1]=-1/(2*dy);
										tmp1=tmp1+1;
										row1[tmp1]=k+n+1;
										col1[tmp1]=k+2*l-1;
										val1[tmp1]=1/(2*dy);
								}
								if(i==M-3)
								{
										tmp1=tmp1+1;
										row1[tmp1]=k+n;
										col1[tmp1]=k+l+m;
										val1[tmp1]=-1/(2*dx);
										tmp1=tmp1+1;
										row1[tmp1]=k+n;
										col1[tmp1]=k+l-m;
										val1[tmp1]=1/(2*dx);
										tmp1=tmp1+1;
										row1[tmp1]=k+n;
										col1[tmp1]=k+2*l+1;
										val1[tmp1]=-1/(2*dy);
										tmp1=tmp1+1;
										row1[tmp1]=k+n;
										col1[tmp1]=k+2*l-1;
										val1[tmp1]=1/(2*dy);
								}
								if(j==N-3)
								{
										tmp1=tmp1+1;
										row1[tmp1]=k+1;
										col1[tmp1]=k+l+m;
										val1[tmp1]=-1/(2*dx);
										tmp1=tmp1+1;
										row1[tmp1]=k+1;
										col1[tmp1]=k+l-m;
										val1[tmp1]=1/(2*dx);
										tmp1=tmp1+1;
										row1[tmp1]=k+1;
										col1[tmp1]=k+2*l+1;
										val1[tmp1]=-1/(2*dy);
										tmp1=tmp1+1;
										row1[tmp1]=k+1;
										col1[tmp1]=k+2*l-1;
										val1[tmp1]=1/(2*dy);
								}
						}
						k=k+l;
						a1 = (x[k+m]+x[k])/(x[k-l+m]+x[k-l]); /*x(50)+x(44)..*/
						a2 = (x[k]+x[k-m])/(x[k-l]+x[k-l-m]); /*(44)+x(38)..*/
						a3 = (x[k+l+1]+x[k+l])/(x[k-l+1]+x[k-l]); /*%x(81)+x(80)..*/
						a4 = (x[k+l-1]+x[k+l])/(x[k-l-1]+x[k-l]); /*%x(79)+x(80)..*/
						a5 = (x[k-l+m]+x[k-l]); /*%x(14)+x(8)*/
						a6 = (x[k-l]+x[k-l-m]); /*%x(8)+x(2)*/
						a7 = (x[k+1]+x[k])/(x[k-l+1]+x[k-l]); /*%x(45)+x(44)..*/
						a8 = (x[k]+x[k-1])/(x[k-l]+x[k-l-1]); /*%x(44)+x(43)..*/
						row[tmp]=k;
						col[tmp]=k+m;
						val[tmp]=-a1/dx;
						tmp=tmp+1;
						row[tmp]=k;
						col[tmp]=k;
						val[tmp]= (a1-a2)*(-1/dx)-(1/(2*dy))*(a3-a4);
						tmp=tmp+1;
						row[tmp]=k;
						col[tmp]=k-l+m;
						val[tmp] = (-1/dx)*(-a1*a1/2+g/4*a5);
						tmp=tmp+1;
						row[tmp]=k;
						col[tmp]=k-l;
						val[tmp] = (-1/dx)*((-a1*a1/2+g/4*a5)-(-a2*a2/2+g/4*a6))- (1/dy)*(-a3*a7/2 +a4*a8/2);
						tmp=tmp+1;
						row[tmp]=k;
						col[tmp]=k-m;
						val[tmp] = a2/(dx);
						tmp=tmp+1;
						row[tmp]=k;
						col[tmp]=k-l-m;
						val[tmp] = (-1/dx)*(a2*a2/2-g/4*a6);
						tmp=tmp+1;
						row[tmp]=k;
						col[tmp]=k+l+1;
						val[tmp] = -a7/(2*dy);
						tmp=tmp+1;
						row[tmp]=k;
						col[tmp]=k+l;
						val[tmp] = (-1/dy)*(a7/2-a8/2);
						tmp=tmp+1;
						row[tmp]=k;
						col[tmp]=k+1;
						val[tmp] = (-1/dy)*(a3/2);
						tmp=tmp+1;
						row[tmp]=k;
						col[tmp]=k+l-1;
						val[tmp]=(1/dy)*a8/2;
						tmp=tmp+1;
						row[tmp]=k;
						col[tmp]=k-1;
						val[tmp] = a4/(2*dy);
						tmp=tmp+1;
						row[tmp]=k;
						col[tmp]=k-l-1;
						val[tmp] = (-1/dy)*(a4*a8/2);
						tmp=tmp+1;
						row[tmp]=k;
						col[tmp]=k-l+1;
						val[tmp] = (1/dy)*(a3*a7/2);
						tmp=tmp+1;
						if(i==0||j==0)
						{
								if(i==0 && j==0)
								{
										tmp1=tmp1+1;
										row1[tmp1]=k-n-1;
										col1[tmp1]=k+m;
										val1[tmp1]=-a1/dx*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k-n-1;
										col1[tmp1]=k;
										val1[tmp1]= ((a1-a2)*(-1/dx)-(1/(2*dy))*(a3-a4))*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k-n-1;
										col1[tmp1]=k-l+m;
										val1[tmp1] = ((-1/dx)*(-a1*a1/2+g/4*a5))*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k-n-1;
										col1[tmp1]=k-l;
										val1[tmp1] = ((-1/dx)*((-a1*a1/2+g/4*a5)-(-a2*a2/2+g/4*a6))- (1/dy)*(-a3*a7/2 +a4*a8/2))*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k-n-1;
										col1[tmp1]=k-m;
										val1[tmp1] = (a2/(dx))*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k-n-1;
										col1[tmp1]=k-l-m;
										val1[tmp1] = ((-1/dx)*(a2*a2/2-g/4*a6))*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k-n-1;
										col1[tmp1]=k+l+1;
										val1[tmp1] = (-a7/(2*dy))*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k-n-1;
										col1[tmp1]=k+l;
										val1[tmp1] = ((-1/dy)*(a7/2-a8/2))*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k-n-1;
										col1[tmp1]=k+1;
										val1[tmp1] = ((-1/dy)*(a3/2))*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k-n-1;
										col1[tmp1]=k+l-1;
										val1[tmp1]=(1/dy)*a8/2*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k-n-1;
										col1[tmp1]=k-1;
										val1[tmp1] = a4/(2*dy)*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k-n-1;
										col1[tmp1]=k-l-1;
										val1[tmp1] = ((-1/dy)*(a4*a8/2))*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k-n-1;
										col1[tmp1]=k-l+1;
										val1[tmp1] = ((1/dy)*(a3*a7/2))*-1;
								}
						

								if(i==0)
								{
										tmp1=tmp1+1;
										row1[tmp1]=k-n;
										col1[tmp1]=k+m;
										val1[tmp1]=-a1/dx*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k-n;
										col1[tmp1]=k;
										val1[tmp1]= ((a1-a2)*(-1/dx)-(1/(2*dy))*(a3-a4))*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k-n;
										col1[tmp1]=k-l+m;
										val1[tmp1] = ((-1/dx)*(-a1*a1/2+g/4*a5))*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k-n;
										col1[tmp1]=k-l;
										val1[tmp1] = ((-1/dx)*((-a1*a1/2+g/4*a5)-(-a2*a2/2+g/4*a6))- (1/dy)*(-a3*a7/2 +a4*a8/2))*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k-n;
										col1[tmp1]=k-m;
										val1[tmp1] = (a2/(dx))*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k-n;
										col1[tmp1]=k-l-m;
										val1[tmp1] = ((-1/dx)*(a2*a2/2-g/4*a6))*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k-n;
										col1[tmp1]=k+l+1;
										val1[tmp1] = (-a7/(2*dy))*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k-n;
										col1[tmp1]=k+l;
										val1[tmp1] = ((-1/dy)*(a7/2-a8/2))*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k-n;
										col1[tmp1]=k+1;
										val1[tmp1] = ((-1/dy)*(a3/2))*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k-n;
										col1[tmp1]=k+l-1;
										val1[tmp1]=(1/dy)*a8/2*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k-n;
										col1[tmp1]=k-1;
										val1[tmp1] = a4/(2*dy)*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k-n;
										col1[tmp1]=k-l-1;
										val1[tmp1] = ((-1/dy)*(a4*a8/2))*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k-n;
										col1[tmp1]=k-l+1;
										val1[tmp1] = ((1/dy)*(a3*a7/2))*-1;

								}
								if(j==0)
								{
										tmp1=tmp1+1;
										row1[tmp1]=k-1;
										col1[tmp1]=k+m;
										val1[tmp1]=-a1/dx;
										tmp1=tmp1+1;
										row1[tmp1]=k-1;
										col1[tmp1]=k;
										val1[tmp1]= (a1-a2)*(-1/dx)-(1/(2*dy))*(a3-a4);
										tmp1=tmp1+1;
										row1[tmp1]=k-1;
										col1[tmp1]=k-l+m;
										val1[tmp1] = (-1/dx)*(-a1*a1/2+g/4*a5);
										tmp1=tmp1+1;
										row1[tmp1]=k-1;
										col1[tmp1]=k-l;
										val1[tmp1] = (-1/dx)*((-a1*a1/2+g/4*a5)-(-a2*a2/2+g/4*a6))- (1/dy)*(-a3*a7/2 +a4*a8/2);
										tmp1=tmp1+1;
										row1[tmp1]=k-1;
										col1[tmp1]=k-m;
										val1[tmp1] = a2/(dx);
										tmp1=tmp1+1;
										row1[tmp1]=k-1;
										col1[tmp1]=k-l-m;
										val1[tmp1] = (-1/dx)*(a2*a2/2-g/4*a6);
										tmp1=tmp1+1;
										row1[tmp1]=k-1;
										col1[tmp1]=k+l+1;
										val1[tmp1] = -a7/(2*dy);
										tmp1=tmp1+1;
										row1[tmp1]=k-1;
										col1[tmp1]=k+l;
										val1[tmp1] = (-1/dy)*(a7/2-a8/2);
										tmp1=tmp1+1;
										row1[tmp1]=k-1;
										col1[tmp1]=k+1;
										val1[tmp1] = (-1/dy)*(a3/2);
										tmp1=tmp1+1;
										row1[tmp1]=k-1;
										col1[tmp1]=k+l-1;
										val1[tmp1]=(1/dy)*a8/2;
										tmp1=tmp1+1;
										row1[tmp1]=k-1;
										col1[tmp1]=k-1;
										val1[tmp1] = a4/(2*dy);
										tmp1=tmp1+1;
										row1[tmp1]=k-1;
										col1[tmp1]=k-l-1;
										val1[tmp1] = (-1/dy)*(a4*a8/2);
										tmp1=tmp1+1;
										row1[tmp1]=k-1;
										col1[tmp1]=k-l+1;
										val1[tmp1] = (1/dy)*(a3*a7/2);

								}
								if(i==0 && j==N-3)
								{
										tmp1=tmp1+1;
										row1[tmp1]=k-n+1;
										col1[tmp1]=k+m;
										val1[tmp1]=-a1/dx*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k-n+1;
										col1[tmp1]=k;
										val1[tmp1]= ((a1-a2)*(-1/dx)-(1/(2*dy))*(a3-a4))*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k-n+1;
										col1[tmp1]=k-l+m;
										val1[tmp1] = ((-1/dx)*(-a1*a1/2+g/4*a5))*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k-n+1;
										col1[tmp1]=k-l;
										val1[tmp1] = ((-1/dx)*((-a1*a1/2+g/4*a5)-(-a2*a2/2+g/4*a6))- (1/dy)*(-a3*a7/2 +a4*a8/2))*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k-n+1;
										col1[tmp1]=k-m;
										val1[tmp1] = (a2/(dx))*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k-n+1;
										col1[tmp1]=k-l-m;
										val1[tmp1] = ((-1/dx)*(a2*a2/2-g/4*a6))*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k-n+1;
										col1[tmp1]=k+l+1;
										val1[tmp1] = (-a7/(2*dy))*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k-n+1;
										col1[tmp1]=k+l;
										val1[tmp1] = ((-1/dy)*(a7/2-a8/2))*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k-n+1;
										col1[tmp1]=k+1;
										val1[tmp1] = ((-1/dy)*(a3/2))*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k-n+1;
										col1[tmp1]=k+l-1;
										val1[tmp1]=(1/dy)*a8/2*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k-n+1;
										col1[tmp1]=k-1;
										val1[tmp1] = a4/(2*dy)*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k-n+1;
										col1[tmp1]=k-l-1;
										val1[tmp1] = ((-1/dy)*(a4*a8/2))*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k-n+1;
										col1[tmp1]=k-l+1;
										val1[tmp1] = ((1/dy)*(a3*a7/2))*-1;
								}
								if(i==M-3 && j==0)
								{
										tmp1=tmp1+1;
										row1[tmp1]=k+n-1;
										col1[tmp1]=k+m;
										val1[tmp1]=-a1/dx*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k+n-1;
										col1[tmp1]=k;
										val1[tmp1]= ((a1-a2)*(-1/dx)-(1/(2*dy))*(a3-a4))*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k+n-1;
										col1[tmp1]=k-l+m;
										val1[tmp1] = ((-1/dx)*(-a1*a1/2+g/4*a5))*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k+n-1;
										col1[tmp1]=k-l;
										val1[tmp1] = ((-1/dx)*((-a1*a1/2+g/4*a5)-(-a2*a2/2+g/4*a6))- (1/dy)*(-a3*a7/2 +a4*a8/2))*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k+n-1;
										col1[tmp1]=k-m;
										val1[tmp1] = (a2/(dx))*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k+n-1;
										col1[tmp1]=k-l-m;
										val1[tmp1] = ((-1/dx)*(a2*a2/2-g/4*a6))*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k+n-1;
										col1[tmp1]=k+l+1;
										val1[tmp1] = (-a7/(2*dy))*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k+n-1;
										col1[tmp1]=k+l;
										val1[tmp1] = ((-1/dy)*(a7/2-a8/2))*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k+n-1;
										col1[tmp1]=k+1;
										val1[tmp1] = ((-1/dy)*(a3/2))*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k+n-1;
										col1[tmp1]=k+l-1;
										val1[tmp1]=(1/dy)*a8/2*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k+n-1;
										col1[tmp1]=k-1;
										val1[tmp1] = a4/(2*dy)*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k+n-1;
										col1[tmp1]=k-l-1;
										val1[tmp1] = ((-1/dy)*(a4*a8/2))*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k+n-1;
										col1[tmp1]=k-l+1;
										val1[tmp1] = ((1/dy)*(a3*a7/2))*-1;
								}


						}
						if(i==M-3 || j==N-3)
						{
								if(i==M-3&&j==N-3)
								{
										tmp1=tmp1+1;
										row1[tmp1]=k+n+1;
										col1[tmp1]=k+m;
										val1[tmp1]=-a1/dx*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k+n+1;
										col1[tmp1]=k;
										val1[tmp1]= ((a1-a2)*(-1/dx)-(1/(2*dy))*(a3-a4))*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k+n+1;
										col1[tmp1]=k-l+m;
										val1[tmp1] = ((-1/dx)*(-a1*a1/2+g/4*a5))*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k+n+1;
										col1[tmp1]=k-l;
										val1[tmp1] = ((-1/dx)*((-a1*a1/2+g/4*a5)-(-a2*a2/2+g/4*a6))- (1/dy)*(-a3*a7/2 +a4*a8/2))*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k+n+1;
										col1[tmp1]=k-m;
										val1[tmp1] = (a2/(dx))*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k+n+1;
										col1[tmp1]=k-l-m;
										val1[tmp1] = ((-1/dx)*(a2*a2/2-g/4*a6))*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k+n+1;
										col1[tmp1]=k+l+1;
										val1[tmp1] = (-a7/(2*dy))*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k+n+1;
										col1[tmp1]=k+l;
										val1[tmp1] = ((-1/dy)*(a7/2-a8/2))*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k+n+1;
										col1[tmp1]=k+1;
										val1[tmp1] = ((-1/dy)*(a3/2))*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k+n+1;
										col1[tmp1]=k+l-1;
										val1[tmp1]=(1/dy)*a8/2*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k+n+1;
										col1[tmp1]=k-1;
										val1[tmp1] = a4/(2*dy)*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k+n+1;
										col1[tmp1]=k-l-1;
										val1[tmp1] = ((-1/dy)*(a4*a8/2))*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k+n+1;
										col1[tmp1]=k-l+1;
										val1[tmp1] = ((1/dy)*(a3*a7/2))*-1;
								}
								if(i==M-3)
								{
										tmp1=tmp1+1;
										row1[tmp1]=k+n;
										col1[tmp1]=k+m;
										val1[tmp1]=-a1/dx*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k+n;
										col1[tmp1]=k;
										val1[tmp1]= ((a1-a2)*(-1/dx)-(1/(2*dy))*(a3-a4))*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k+n;
										col1[tmp1]=k-l+m;
										val1[tmp1] = ((-1/dx)*(-a1*a1/2+g/4*a5))*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k+n;
										col1[tmp1]=k-l;
										val1[tmp1] = ((-1/dx)*((-a1*a1/2+g/4*a5)-(-a2*a2/2+g/4*a6))- (1/dy)*(-a3*a7/2 +a4*a8/2))*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k+n;
										col1[tmp1]=k-m;
										val1[tmp1] = (a2/(dx))*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k+n;
										col1[tmp1]=k-l-m;
										val1[tmp1] = ((-1/dx)*(a2*a2/2-g/4*a6))*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k+n;
										col1[tmp1]=k+l+1;
										val1[tmp1] = (-a7/(2*dy))*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k+n;
										col1[tmp1]=k+l;
										val1[tmp1] = ((-1/dy)*(a7/2-a8/2))*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k+n;
										col1[tmp1]=k+1;
										val1[tmp1] = ((-1/dy)*(a3/2))*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k+n;
										col1[tmp1]=k+l-1;
										val1[tmp1]=(1/dy)*a8/2*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k+n;
										col1[tmp1]=k-1;
										val1[tmp1] = a4/(2*dy)*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k+n;
										col1[tmp1]=k-l-1;
										val1[tmp1] = ((-1/dy)*(a4*a8/2))*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k+n;
										col1[tmp1]=k-l+1;
										val1[tmp1] = ((1/dy)*(a3*a7/2))*-1;
								}
								if(j==N-3)
								{
										tmp1=tmp1+1;
										row1[tmp1]=k+1;
										col1[tmp1]=k+m;
										val1[tmp1]=-a1/dx;
										tmp1=tmp1+1;
										row1[tmp1]=k+1;
										col1[tmp1]=k;
										val1[tmp1]= (a1-a2)*(-1/dx)-(1/(2*dy))*(a3-a4);
										tmp1=tmp1+1;
										row1[tmp1]=k+1;
										col1[tmp1]=k-l+m;
										val1[tmp1] = (-1/dx)*(-a1*a1/2+g/4*a5);
										tmp1=tmp1+1;
										row1[tmp1]=k+1;
										col1[tmp1]=k-l;
										val1[tmp1] = (-1/dx)*((-a1*a1/2+g/4*a5)-(-a2*a2/2+g/4*a6))- (1/dy)*(-a3*a7/2 +a4*a8/2);
										tmp1=tmp1+1;
										row1[tmp1]=k+1;
										col1[tmp1]=k-m;
										val1[tmp1] = a2/(dx);
										tmp1=tmp1+1;
										row1[tmp1]=k+1;
										col1[tmp1]=k-l-m;
										val1[tmp1] = (-1/dx)*(a2*a2/2-g/4*a6);
										tmp1=tmp1+1;
										row1[tmp1]=k+1;
										col1[tmp1]=k+l+1;
										val1[tmp1] = -a7/(2*dy);
										tmp1=tmp1+1;
										row1[tmp1]=k+1;
										col1[tmp1]=k+l;
										val1[tmp1] = (-1/dy)*(a7/2-a8/2);
										tmp1=tmp1+1;
										row1[tmp1]=k+1;
										col1[tmp1]=k+1;
										val1[tmp1] = (-1/dy)*(a3/2);
										tmp1=tmp1+1;
										row1[tmp1]=k+1;
										col1[tmp1]=k+l-1;
										val1[tmp1]=(1/dy)*a8/2;
										tmp1=tmp1+1;
										row1[tmp1]=k+1;
										col1[tmp1]=k-1;
										val1[tmp1] = a4/(2*dy);
										tmp1=tmp1+1;
										row1[tmp1]=k+1;
										col1[tmp1]=k-l-1;
										val1[tmp1] = (-1/dy)*(a4*a8/2);
										tmp1=tmp1+1;
										row1[tmp1]=k+1;
										col1[tmp1]=k-l+1;
										val1[tmp1] = (1/dy)*(a3*a7/2);
								}


						}
						k = k+l;
						a9 = (x[k+m]+x[k])/(x[k-2*l+m]+x[k-2*l]); /*%x(80)+x(86)..*/
						a10 =(x[k]+x[k-m])/(x[k-2*l]+x[k-2*l-m]); /*%x(80)+x(74)..*/
						a11 =(x[k-2*l]+x[k-2*l+1]); /*%x(9)+x(8)*/
						a12 =(x[k-2*l]+x[k-2*l-1]); /*%x(8)+x(7)*/
						row[tmp]=k;
						col[tmp]=k-l+m;
						val[tmp] = -a9/(2*dx);
						tmp=tmp+1;
						row[tmp]=k;
						col[tmp]=k-l;
						val[tmp] = (-1/(2*dx))*(a9-a10);
						tmp=tmp+1;
						row[tmp]=k;
						col[tmp]=k;
						val[tmp] = (-1/(2*dx))*(a1-a2)-(1/dy)*(a3-a4);
						tmp=tmp+1;
						row[tmp]=k;
						col[tmp]=k+m;
						val[tmp] = -a1/(2*dx);
						tmp=tmp+1;
						row[tmp]=k;
						col[tmp]=k-2*l+m;
						val[tmp] = (a1*a9)/(2*dx);
						tmp=tmp+1;
						row[tmp]=k;
						col[tmp]=k-2*l;
						val[tmp] = (-1/dx)*(-a1*a9/2+a2*a10/2)-(1/dy)*(-a3*a3/2+g/4*(a11)+a4*a4/2-g/4*a12);
						tmp=tmp+1;
						row[tmp]=k;
						col[tmp]=k-l-m;
						val[tmp] = (a10/(2*dx));
						tmp=tmp+1;
						row[tmp]=k;
						col[tmp]=k-m;
						val[tmp] = (a2/(2*dx));
						tmp=tmp+1;
						row[tmp]=k;
						col[tmp]=k-2*l-m;
						val[tmp] = (-1/(2*dx))*(a2*a10);
						tmp=tmp+1;
						row[tmp]=k;
						col[tmp]=k+1;
						val[tmp] = -a3/dy;
						tmp=tmp+1;
						row[tmp]=k;
						col[tmp]=k-1;
						val[tmp] = a4/dy;
						tmp=tmp+1;
						row[tmp]=k;
						col[tmp]=k-2*l+1;
						val[tmp] = (-1/dy)*(-a3*a3/2+g/4*a11);
						tmp=tmp+1;
						row[tmp]=k;
						col[tmp]=k-2*l-1;
						val[tmp] = (-1/dy)*(a4*a4/2-g/4*a12);
						if(i==0 || j==0)
						{
								if(i==0 && j==0)
								{
										tmp1=tmp1+1;
										row1[tmp1]=k-n-1;
										col1[tmp1]=k-l+m;
										val1[tmp1] = -a9/(2*dx)*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k-n-1;
										col1[tmp1]=k-l;
										val1[tmp1] = ((-1/(2*dx))*(a9-a10))*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k-n-1;
										col1[tmp1]=k;
										val1[tmp1] = ((-1/(2*dx))*(a1-a2)-(1/dy)*(a3-a4))*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k-n-1;
										col1[tmp1]=k+m;
										val1[tmp1] = (-a1/(2*dx))*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k-n-1;
										col1[tmp1]=k-2*l+m;
										val1[tmp1] = ((a1*a9)/(2*dx))*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k-n-1;
										col1[tmp1]=k-2*l;
										val1[tmp1] = ((-1/dx)*(-a1*a9/2+a2*a10/2)-(1/dy)*(-a3*a3/2+g/4*(a11)+a4*a4/2-g/4*a12))*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k-n-1;
										col1[tmp1]=k-l-m;
										val1[tmp1] = ((a10/(2*dx)))*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k-n-1;
										col1[tmp1]=k-m;
										val1[tmp1] = ((a2/(2*dx)))*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k-n-1;
										col1[tmp1]=k-2*l-m;
										val1[tmp1] = ((-1/(2*dx))*(a2*a10))*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k-n-1;
										col1[tmp1]=k+1;
										val1[tmp1] = (-a3/dy)*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k-n-1;
										col1[tmp1]=k-1;
										val1[tmp1] = (a4/dy)*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k-n-1;
										col1[tmp1]=k-2*l+1;
										val1[tmp1] = ((-1/dy)*(-a3*a3/2+g/4*a11))*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k-n-1;
										col1[tmp1]=k-2*l-1;
										val1[tmp1] = ((-1/dy)*(a4*a4/2-g/4*a12))*-1;

								}
								if(i==0)
								{
										tmp1=tmp1+1;
										row1[tmp1]=k-n;
										col1[tmp1]=k-l+m;
										val1[tmp1] = -a9/(2*dx);
										tmp1=tmp1+1;
										row1[tmp1]=k-n;
										col1[tmp1]=k-l;
										val1[tmp1] = (-1/(2*dx))*(a9-a10);
										tmp1=tmp1+1;
										row1[tmp1]=k-n;
										col1[tmp1]=k;
										val1[tmp1] = (-1/(2*dx))*(a1-a2)-(1/dy)*(a3-a4);
										tmp1=tmp1+1;
										row1[tmp1]=k-n;
										col1[tmp1]=k+m;
										val1[tmp1] = -a1/(2*dx);
										tmp1=tmp1+1;
										row1[tmp1]=k-n;
										col1[tmp1]=k-2*l+m;
										val1[tmp1] = (a1*a9)/(2*dx);
										tmp1=tmp1+1;
										row1[tmp1]=k-n;
										col1[tmp1]=k-2*l;
										val1[tmp1] = (-1/dx)*(-a1*a9/2+a2*a10/2)-(1/dy)*(-a3*a3/2+g/4*(a11)+a4*a4/2-g/4*a12);
										tmp1=tmp1+1;
										row1[tmp1]=k-n;
										col1[tmp1]=k-l-m;
										val1[tmp1] = (a10/(2*dx));
										tmp1=tmp1+1;
										row1[tmp1]=k-n;
										col1[tmp1]=k-m;
										val1[tmp1] = (a2/(2*dx));
										tmp1=tmp1+1;
										row1[tmp1]=k-n;
										col1[tmp1]=k-2*l-m;
										val1[tmp1] = (-1/(2*dx))*(a2*a10);
										tmp1=tmp1+1;
										row1[tmp1]=k-n;
										col1[tmp1]=k+1;
										val1[tmp1] = -a3/dy;
										tmp1=tmp1+1;
										row1[tmp1]=k-n;
										col1[tmp1]=k-1;
										val1[tmp1] = a4/dy;
										tmp1=tmp1+1;
										row1[tmp1]=k-n;
										col1[tmp1]=k-2*l+1;
										val1[tmp1] = (-1/dy)*(-a3*a3/2+g/4*a11);
										tmp1=tmp1+1;
										row1[tmp1]=k-n;
										col1[tmp1]=k-2*l-1;
										val1[tmp1] = (-1/dy)*(a4*a4/2-g/4*a12);
								}
								if(j==0)
								{
										tmp1=tmp1+1;
										row1[tmp1]=k-1;
										col1[tmp1]=k-l+m;
										val1[tmp1] = (-a9/(2*dx))*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k-1;
										col1[tmp1]=k-l;
										val1[tmp1] = ((-1/(2*dx))*(a9-a10))*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k-1;
										col1[tmp1]=k;
										val1[tmp1] = ((-1/(2*dx))*(a1-a2)-(1/dy)*(a3-a4))*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k-1;
										col1[tmp1]=k+m;
										val1[tmp1] = (-a1/(2*dx))*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k-1;
										col1[tmp1]=k-2*l+m;
										val1[tmp1] = ((a1*a9)/(2*dx))*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k-1;
										col1[tmp1]=k-2*l;
										val1[tmp1] = ((-1/dx)*(-a1*a9/2+a2*a10/2)-(1/dy)*(-a3*a3/2+g/4*(a11)+a4*a4/2-g/4*a12))*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k-1;
										col1[tmp1]=k-l-m;
										val1[tmp1] = ((a10/(2*dx)))*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k-1;
										col1[tmp1]=k-m;
										val1[tmp1] = ((a2/(2*dx)))*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k-1;
										col1[tmp1]=k-2*l-m;
										val1[tmp1] = ((-1/(2*dx))*(a2*a10))*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k-1;
										col1[tmp1]=k+1;
										val1[tmp1] = (-a3/dy)*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k-1;
										col1[tmp1]=k-1;
										val1[tmp1] = (a4/dy)*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k-1;
										col1[tmp1]=k-2*l+1;
										val1[tmp1] = ((-1/dy)*(-a3*a3/2+g/4*a11))*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k-1;
										col1[tmp1]=k-2*l-1;
										val1[tmp1] = ((-1/dy)*(a4*a4/2-g/4*a12))*-1;
											
								}
								if(i==0 && j==N-3)
								{
										tmp1=tmp1+1;
										row1[tmp1]=k-n+1;
										col1[tmp1]=k-l+m;
										val1[tmp1] = -a9/(2*dx)*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k-n+1;
										col1[tmp1]=k-l;
										val1[tmp1] = ((-1/(2*dx))*(a9-a10))*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k-n+1;
										col1[tmp1]=k;
										val1[tmp1] = ((-1/(2*dx))*(a1-a2)-(1/dy)*(a3-a4))*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k-n+1;
										col1[tmp1]=k+m;
										val1[tmp1] = (-a1/(2*dx))*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k-n+1;
										col1[tmp1]=k-2*l+m;
										val1[tmp1] = ((a1*a9)/(2*dx))*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k-n+1;
										col1[tmp1]=k-2*l;
										val1[tmp1] = ((-1/dx)*(-a1*a9/2+a2*a10/2)-(1/dy)*(-a3*a3/2+g/4*(a11)+a4*a4/2-g/4*a12))*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k-n+1;
										col1[tmp1]=k-l-m;
										val1[tmp1] = ((a10/(2*dx)))*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k-n+1;
										col1[tmp1]=k-m;
										val1[tmp1] = ((a2/(2*dx)))*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k-n+1;
										col1[tmp1]=k-2*l-m;
										val1[tmp1] = ((-1/(2*dx))*(a2*a10))*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k-n+1;
										col1[tmp1]=k+1;
										val1[tmp1] = (-a3/dy)*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k-n+1;
										col1[tmp1]=k-1;
										val1[tmp1] = (a4/dy)*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k-n+1;
										col1[tmp1]=k-2*l+1;
										val1[tmp1] = ((-1/dy)*(-a3*a3/2+g/4*a11))*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k-n+1;
										col1[tmp1]=k-2*l-1;
										val1[tmp1] = ((-1/dy)*(a4*a4/2-g/4*a12))*-1;
								}
								if(i==N-3 && j==0)
								{

										tmp1=tmp1+1;
										row1[tmp1]=k+n-1;
										col1[tmp1]=k-l+m;
										val1[tmp1] = -a9/(2*dx)*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k+n-1;
										col1[tmp1]=k-l;
										val1[tmp1] = ((-1/(2*dx))*(a9-a10))*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k+n-1;
										col1[tmp1]=k;
										val1[tmp1] = ((-1/(2*dx))*(a1-a2)-(1/dy)*(a3-a4))*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k+n-1;
										col1[tmp1]=k+m;
										val1[tmp1] = (-a1/(2*dx))*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k+n-1;
										col1[tmp1]=k-2*l+m;
										val1[tmp1] = ((a1*a9)/(2*dx))*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k+n-1;
										col1[tmp1]=k-2*l;
										val1[tmp1] = ((-1/dx)*(-a1*a9/2+a2*a10/2)-(1/dy)*(-a3*a3/2+g/4*(a11)+a4*a4/2-g/4*a12))*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k+n-1;
										col1[tmp1]=k-l-m;
										val1[tmp1] = ((a10/(2*dx)))*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k+n-1;
										col1[tmp1]=k-m;
										val1[tmp1] = ((a2/(2*dx)))*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k+n-1;
										col1[tmp1]=k-2*l-m;
										val1[tmp1] = ((-1/(2*dx))*(a2*a10))*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k+n-1;
										col1[tmp1]=k+1;
										val1[tmp1] = (-a3/dy)*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k+n-1;
										col1[tmp1]=k-1;
										val1[tmp1] = (a4/dy)*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k+n-1;
										col1[tmp1]=k-2*l+1;
										val1[tmp1] = ((-1/dy)*(-a3*a3/2+g/4*a11))*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k+n-1;
										col1[tmp1]=k-2*l-1;
										val1[tmp1] = ((-1/dy)*(a4*a4/2-g/4*a12))*-1;

								}

						}
						if(i==M-3 || j==N-3)
						{
								if(i==M-3 && j==N-3)
								{
										tmp1=tmp1+1;
										row1[tmp1]=k+n+1;
										col1[tmp1]=k-l+m;
										val1[tmp1] = (-a9/(2*dx))*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k+n+1;
										col1[tmp1]=k-l;
										val1[tmp1] = ((-1/(2*dx))*(a9-a10))*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k+n+1;
										col1[tmp1]=k;
										val1[tmp1] = ((-1/(2*dx))*(a1-a2)-(1/dy)*(a3-a4))*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k+n+1;
										col1[tmp1]=k+m;
										val1[tmp1] = (-a1/(2*dx))*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k+n+1;
										col1[tmp1]=k-2*l+m;
										val1[tmp1] = ((a1*a9)/(2*dx))*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k+n+1;
										col1[tmp1]=k-2*l;
										val1[tmp1] = ((-1/dx)*(-a1*a9/2+a2*a10/2)-(1/dy)*(-a3*a3/2+g/4*(a11)+a4*a4/2-g/4*a12))*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k+n+1;
										col1[tmp1]=k-l-m;
										val1[tmp1] = ((a10/(2*dx)))*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k+n+1;
										col1[tmp1]=k-m;
										val1[tmp1] = ((a2/(2*dx)))*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k+n+1;
										col1[tmp1]=k-2*l-m;
										val1[tmp1] = ((-1/(2*dx))*(a2*a10))*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k+n+1;
										col1[tmp1]=k+1;
										val1[tmp1] = (-a3/dy)*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k+n+1;
										col1[tmp1]=k-1;
										val1[tmp1] = (a4/dy)*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k+n+1;
										col1[tmp1]=k-2*l+1;
										val1[tmp1] = ((-1/dy)*(-a3*a3/2+g/4*a11))*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k+n+1;
										col1[tmp1]=k-2*l-1;
										val1[tmp1] = ((-1/dy)*(a4*a4/2-g/4*a12))*-1;
								}
								if(i==M-3)
								{
										tmp1=tmp1+1;
										row1[tmp1]=k+n;
										col1[tmp1]=k-l+m;
										val1[tmp1] = -a9/(2*dx);
										tmp1=tmp1+1;
										row1[tmp1]=k+n;
										col1[tmp1]=k-l;
										val1[tmp1] = (-1/(2*dx))*(a9-a10);
										tmp1=tmp1+1;
										row1[tmp1]=k+n;
										col1[tmp1]=k;
										val1[tmp1] = (-1/(2*dx))*(a1-a2)-(1/dy)*(a3-a4);
										tmp1=tmp1+1;
										row1[tmp1]=k+n;
										col1[tmp1]=k+m;
										val1[tmp1] = -a1/(2*dx);
										tmp1=tmp1+1;
										row1[tmp1]=k+n;
										col1[tmp1]=k-2*l+m;
										val1[tmp1] = (a1*a9)/(2*dx);
										tmp1=tmp1+1;
										row1[tmp1]=k+n;
										col1[tmp1]=k-2*l;
										val1[tmp1] = (-1/dx)*(-a1*a9/2+a2*a10/2)-(1/dy)*(-a3*a3/2+g/4*(a11)+a4*a4/2-g/4*a12);
										tmp1=tmp1+1;
										row1[tmp1]=k+n;
										col1[tmp1]=k-l-m;
										val1[tmp1] = (a10/(2*dx));
										tmp1=tmp1+1;
										row1[tmp1]=k+n;
										col1[tmp1]=k-m;
										val1[tmp1] = (a2/(2*dx));
										tmp1=tmp1+1;
										row1[tmp1]=k+n;
										col1[tmp1]=k-2*l-m;
										val1[tmp1] = (-1/(2*dx))*(a2*a10);
										tmp1=tmp1+1;
										row1[tmp1]=k+n;
										col1[tmp1]=k+1;
										val1[tmp1] = -a3/dy;
										tmp1=tmp1+1;
										row1[tmp1]=k+n;
										col1[tmp1]=k-1;
										val1[tmp1] = a4/dy;
										tmp1=tmp1+1;
										row1[tmp1]=k+n;
										col1[tmp1]=k-2*l+1;
										val1[tmp1] = (-1/dy)*(-a3*a3/2+g/4*a11);
										tmp1=tmp1+1;
										row1[tmp1]=k+n;
										col1[tmp1]=k-2*l-1;
										val1[tmp1] = (-1/dy)*(a4*a4/2-g/4*a12);

								}
								if(j==N-3)
								{
										tmp1=tmp1+1;
										row1[tmp1]=k+1;
										col1[tmp1]=k-l+m;
										val1[tmp1] = (-a9/(2*dx))*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k+1;
										col1[tmp1]=k-l;
										val1[tmp1] = ((-1/(2*dx))*(a9-a10))*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k+1;
										col1[tmp1]=k;
										val1[tmp1] = ((-1/(2*dx))*(a1-a2)-(1/dy)*(a3-a4))*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k+1;
										col1[tmp1]=k+m;
										val1[tmp1] = (-a1/(2*dx))*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k+1;
										col1[tmp1]=k-2*l+m;
										val1[tmp1] = ((a1*a9)/(2*dx))*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k+1;
										col1[tmp1]=k-2*l;
										val1[tmp1] = ((-1/dx)*(-a1*a9/2+a2*a10/2)-(1/dy)*(-a3*a3/2+g/4*(a11)+a4*a4/2-g/4*a12))*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k+1;
										col1[tmp1]=k-l-m;
										val1[tmp1] = ((a10/(2*dx)))*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k+1;
										col1[tmp1]=k-m;
										val1[tmp1] = ((a2/(2*dx)))*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k+1;
										col1[tmp1]=k-2*l-m;
										val1[tmp1] = ((-1/(2*dx))*(a2*a10))*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k+1;
										col1[tmp1]=k+1;
										val1[tmp1] = (-a3/dy)*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k+1;
										col1[tmp1]=k-1;
										val1[tmp1] = (a4/dy)*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k+1;
										col1[tmp1]=k-2*l+1;
										val1[tmp1] = ((-1/dy)*(-a3*a3/2+g/4*a11))*-1;
										tmp1=tmp1+1;
										row1[tmp1]=k+1;
										col1[tmp1]=k-2*l-1;
										val1[tmp1] = ((-1/dy)*(a4*a4/2-g/4*a12))*-1;
								}
						}

				}
		}
		/*printf("%d\n", tmp);*/
		tmp=tmp+1;
		tmp1=tmp1+1;
		/*printf("tmp1=%d\n",tmp1);*/
		int tt=tmp+tmp1;
		int AA[1]={tt};
		plhs[0]=mxCreateNumericArray(1,AA,mxDOUBLE_CLASS,mxREAL);
		ROW=mxGetPr(plhs[0]);
		plhs[1]=mxCreateNumericArray(1,AA,mxDOUBLE_CLASS,mxREAL);
		COL=mxGetPr(plhs[1]);
		plhs[2]=mxCreateDoubleMatrix(tt,1,mxREAL);
		VAL=mxGetPr(plhs[2]);
		for(i=0;i<tmp;i++)
		{
				ROW[i]=row[i]+1;
				COL[i]=col[i]+1;
				VAL[i]=val[i];
		}
		for(i=0;i<tmp1;i++)
		{
				ROW[tmp+i]=row1[i]+1;
				COL[tmp+i]=col1[i]+1;
				VAL[tmp+i]=val1[i];
		}

		return;
}










