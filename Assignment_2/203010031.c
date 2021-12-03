#include<stdlib.h>
#include<stdio.h>
#include<math.h>

void PJ (){
	FILE * f;
	int IM=41;
	int JM=81;
	double dxPJ=1.0/IM-1;
	double dyPJ=2.0/JM-1;
	double TPJ[JM][IM];
	double T_GPJ[JM][IM];
	double betaPJ=1.0;
	double btsqrPJ=betaPJ*betaPJ;
	int i,j,iterationPJ;
	iterationPJ=1;
	int count =0; int n=0;
	double error; 
	// Boundary conditions
	for(i=0;i<JM;i++){
		for(j=0;j<IM;j++){
			
	  	  	T_GPJ[0][j] = 250.0;
			T_GPJ[JM-1][j] = 0.0;
			T_GPJ[j][0]=0.0;
			T_GPJ[j][IM-1]=0.0;
			
			TPJ[0][j] = 250.0;
			TPJ[JM-1][j] = 0.0;
			TPJ[j][0]=0.0;
			TPJ[j][IM-1]=0.0;
			
		}
	}
//// initial conditions
	for(i=1;i < (JM-2);i++){
		for(j=1;j < (IM-2);j++){
			TPJ[i][j]=0.0;
		}
	}
	// iterations
	f=fopen("error versus iterationPJ for point jacobi.txt","w");
	fprintf(f,"Iteration ----- Error \n");
	error=1;
	iterationPJ=0;
	while (error > 0.001)
	{
		iterationPJ=iterationPJ +1;	
		for(i=1;i < (JM-1);i++)
		{
			for(j=1; j < (IM-1);j++)
			{	
				TPJ[i][j]=(1/(2*(1+(btsqrPJ))))*(T_GPJ[i+1][j] + T_GPJ[i-1][j] + (btsqrPJ*(T_GPJ[i][j+1] + T_GPJ[i][j-1])));	
				
			}
		}
		error=0;	
		for(i=1;i<JM-1;i++)
		{
			for(j=1;j<IM-1;j++)
			{
				error=error+fabs(TPJ[i][j]-T_GPJ[i][j]);
				T_GPJ[i][j]  = TPJ[i][j];
				
			}
		}
		fprintf(f,"%d    %lf\n",iterationPJ,error);				
	}
	printf("1.Point jacobi conver ges at an iterationPJ = %d with an error of %lf \n",iterationPJ,error);
	 fclose(f);
	
	f=fopen("Matrix point jacobi.dat","w");
	for(i=0;i<=(JM-1);i++)
	{
		for(j=0;j<=(IM-1);j++)
		{
			
			
			fprintf(f,"%d %d %lf\n ",j,i, TPJ[i][j]);
		}
	}
	fclose(f);
	
// center line values	
	f=fopen("Mid-line Y tenperature  jacobi.txt","w");
	for(i=0;i<=(JM-1);i++)
	{
		fprintf(f,"%d   %f\n ",i,TPJ[i][JM/2]);
	}
	fclose(f);
	
	f=fopen("Mid-line X tenperature  jacobi.txt","w");
	
	for(j=0;j<=(IM-1);j++)
	{	
		fprintf(f,"%d   %f\n",i,TPJ[IM/2][j]);
	}
	fclose(f);
}

void PGS (){
	FILE * f;
	int IM=41;
	int JM=81;
	double dxPGS=1.0/(IM-1);
	double dyPGS=2.0/(JM-1);
	double TPGS[JM][IM];
	double T_GPGS[JM][IM];
	double beta=1.0;
	double btsqrPGS=beta*beta;
	int i,j,iterationPGS;
	iterationPGS=1;
	double error; 
	// Boundary condition
	for(i=0;i<JM;i++){
		for(j=0;j<IM;j++){
	  	  	T_GPGS[0][j] = 250.0;
			T_GPGS[JM-1][j] = 0.0;
			T_GPGS[j][0]=0.0;
			T_GPGS[j][IM-1]=0.0;
			
			TPGS[0][j] = 250.0;
			TPGS[JM-1][j] = 0.0;
			TPGS[j][0]=0.0;
			TPGS[j][IM-1]=0.0;
			
		}
	}
// initial conditions
	for(i=1;i < (JM-2);i++)
	{
		for(j=1;j < (IM-2);j++){
			TPGS[i][j]=0.0;
		}
	}
	// iterations
	f=fopen("error versus iterationPGS for point jacobi.txt","w");
	fprintf(f,"Iteration ----- Error \n");
	error=1;
	iterationPGS=0;
	while (error > 0.001)
	{
		iterationPGS=iterationPGS +1;	
		for(i=1;i < (JM-1);i++)
		{
			for(j=1; j < (IM-1);j++)
			{	
				TPGS[i][j]=(1/(2*(1+(btsqrPGS))))*(T_GPGS[i+1][j] + TPGS[i-1][j] + (btsqrPGS*(T_GPGS[i][j+1] + TPGS[i][j-1])));		
			}
		}
		error=0.0;	
		for(i=1;i<JM-1;i++)
		{
			for(j=1;j<IM-1;j++)
			{
				error=error+fabs(TPGS[i][j]-T_GPGS[i][j]);
				T_GPGS[i][j]  = TPGS[i][j];
			}
		}
		fprintf(f,"%d    %lf\n",iterationPGS,error);
					
	}
	printf("Point Gauss seidel converges at an iterationPGS = %d with an error of %lf\n ",iterationPGS,error);
	 fclose(f);
	
	f=fopen("Matrix point jacobi.dat","w");
	
	for(i=0;i<=(JM-1);i++)
	{
		for(j=0;j<=(IM-1);j++)
		{
			
			
			fprintf(f,"%d %d %lf\n ",j,i, TPGS[i][j]);
		}
		
		//fprintf(f,"\n");
	}
	fclose(f);
	
// center line values	
	f=fopen("Mid-line Y tenperature  jacobi.txt","w");
	for(i=0;i<=(JM-1);i++)
	{
		fprintf(f,"%d   %f\n ",i,TPGS[i][JM/2]);
	}
	fclose(f);
	
	f=fopen("Mid-line X tenperature  jacobi.txt","w");
	
	for(j=0;j<=(IM-1);j++)
	{	
		fprintf(f,"%d   %f\n",i,TPGS[IM/2][j]);
	}
	fclose(f);
}


void PSOR (){
	FILE * f;
	int IM=41;
	int JM=81;
	double a,b;
	double dx=1.0/(IM-1);
	double dy=2.0/(JM-1);
	double T[JM][IM];
	double t[JM][IM];
	double T_G[JM][IM];
	double beta=1.0;
	double btsqr = beta*beta;
	double pi = 3.145;
	double omega;
	int i,j,iteration;
	iteration=1;
	double error; 
// Boundary conditions
	for(i=0;i<JM;i++){
		for(j=0;j<IM;j++){
	  	  	T_G[0][j] = 250.0;
			T_G[JM-1][j] = 0.0;
			T_G[j][0]=0.0;
			T_G[j][IM-1]=0.0;
			
			T[0][j] = 250.0;
			T[JM-1][j] = 0.0;
			T[j][0]=0.0;
			T[j][IM-1]=0.0;
		}
	}
// initial conditions
	for(i=1;i < (JM-2);i++){
		for(j=1;j < (IM-2);j++){
			T[i][j]=0.0;
		}
	}
	// iterations
	f=fopen("error versus iteration for Point successive over realxation.txt","w");
	fprintf(f,"Iteration ----- Error \n");
	error=2;
	while (error > 0.001)
	{
		iteration=iteration +1;	
		for(i=1;i < (JM-1);i++)
		{
			for(j=1; j < (IM-1);j++)
			{
				a =(((cos(pi/(IM-1)))+ (btsqr*cos(pi/(JM-1))))/(1+btsqr));
				b=a*a;
	  			omega = (2-(2*sqrt(1-b)))/b;
				T[i][j]=(1-omega)*T_G[i][j]+(omega/(2*(1+(btsqr))))*(T_G[i+1][j] + T[i-1][j] + (btsqr*(T_G[i][j+1] + T[i][j-1])));	
					
			}
		}
		error=0;	
		for(i=1;i<JM-1;i++)
		{
			for(j=1;j<IM-1;j++)
			{
				error=error+fabs(T[i][j]-T_G[i][j]);
				T_G[i][j]  = T[i][j];
				
			}
		}
		
		
		fprintf(f,"%d    %lf\n",iteration,error);
		
	}
	printf("Point successive over relaxation converges at an iteration = %d with an error of %lf\n ",iteration,error);
	fclose(f);
f=fopen("Matrix Point successive over realxation.dat","w");
	
	for(i=0;i<=(JM-1);i++)
	{
		for(j=0;j<=(IM-1);j++)
		{
			fprintf(f,"%d %d %lf\n ",j,i, T[i][j]);	
		}
		
	}
	fclose(f);
	
// center line values	
	f=fopen("Mid-line Y tenperature  Point successive over realxation.txt","w");
	for(i=0;i<=(JM-1);i++)
	{
		fprintf(f,"%d   %f\n ",i,T[i][JM/2]);
	}
	fclose(f);
	
	f=fopen("Mid-line X tenperature  Point successive over realxation.txt","w");
	
	for(j=0;j<=(IM-1);j++)
	{	
		fprintf(f,"%d   %f\n",i,T[IM/2][j]);
	}
	fclose(f);
}



int main()
{
	PJ();	
	PGS();
	PSOR();

	return 0;
}
