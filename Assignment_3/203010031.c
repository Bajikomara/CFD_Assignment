#include<stdio.h>
#include<stdlib.h>
#include<math.h>
/*In all the schemes
I is N.o of grid points 
dt  is the smallest time step
dy is distance between two grid points
dy=h/I-1
h is the distance between two plates
nu is kinematic viscosity
i is used in for loop to advaance
CFL is nu*(dt/(dy*dy))*/

void FTCS() {
	int I=61; 
	double h=0.04,dy=h/60.0;
	double nu=0.000217,CFL=0.5,y;
	double dt=(CFL*dy*dy)/nu;
	int N=2000;
	double U[I][N]; // velocity matrix to march in sapce and time
	int i,t;
	FILE * f;
	// Initial conditions
	for(i=0;i<I;i++)
	{
		
		U[i][0]=0.0;
			
	}
	U[0][0]=40.0;
	
// computing at various points
	for(t=0;t<(N);t++){
		//Boundary conditions
		U[0][t]=40.0;
		U[I-1][t]=0.0;
		// computational velocity matrix computing at various grid points and at various time steps
		for(i=1;i<(I-1);i++)
		{
			U[i][t+1]=U[i][t] + ((1*CFL)*(U[i+1][t]-2*U[i][t]+U[i-1][t]))	;
			
		}
		
	}
	// started writing to file
	f=fopen("Y VS U FTCS .csv","w");
	fprintf(f,"t=0.15,,,t=0.8,,,t=1.5\n\n\n");
	for(i=0;i<I;i++){
		y=i*dy;
		
		fprintf(f,"%lf,%lf,,%lf,%lf,,%lf,%lf\n",y,U[i][147],y,U[i][782],y,U[i][1465]);
	}
	fclose(f);	
}


void CNS(){
	int	I=61,t;
	double CFL;
	int i,j,n;
	double h=0.04,dy=0.04/60,y;
	double dt=0.01;
	double C[I],A[I][I],u[I],U[I];
	/* ======================================================
	Matrices in the form AX=D (crank nicholson matrix)
	Notations are followed as explained in Anderson text book 
	C[I] is a matrix representing D
	A[I][I] is a matrix reresenting A
	u[I] is a matrix representing X
	==========================================================*/
	double nu= 0.000217;
	double error;
	CFL=nu*(dt/(dy*dy));
	double a = CFL*0.5,b=1+CFL; // a and b are diagonal elemnts in a crank nicholson matrix "A"
	FILE * f;
	for(i=0;i<I;i++)
	{
		u[0]=40.0;
		u[i]=0.0;
	}	
	for(i=0;i<I-2;i++)
		{
			for(j=0;j<I-2;j++)
			{
				A[i][j]= 0.0;
			}
		}
	error = 10;	
	t=0;
	while (error>0.001){
		t=t+1;
		for(i=0;i<I-2;i++){
			if (i==0){
				C[i]= -a*u[i+1] -u[i+1] -a*(u[i+2] -(2*u[i+1]) + u[i]);
				
			}else if(i==I-3){
				C[i]= -a*u[I-1] -u[i+1] -a*(u[i+2]-(2*u[i+1]) + u[i] ) ;
			}else{
				C[i]= -u[i+1] -a*(u[i+2]-(2*u[i+1])  + u[i]) ;
			}
		}
		
		for (i=0;i<I-2;i++){
			if (i==0){
				A[i][i]=-b;
				A[i][i+1]=a;
			}else if (i==I-3){
				A[I-3][I-4]=a;
				A[I-3][I-3]=-b;
			}else {
				A[i][i-1]=a;
				A[i][i]=-b;
				A[i][i+1]=a;
			}
		}
		for(i=1;i<I-2;i++){
			A[i][i]=A[i][i] - (A[i][i-1]*A[i-1][i]/A[i-1][i-1]) ;
			C[i]= C[i] - (C[i-1]*A[i][i-1]/A[i-1][i-1]);
			A[i][i-1]=0;	
		}
		for(i=I-3;i>=0;i--){
			if(i==I-3){
				U[i+1]=C[i]/A[i][i];
			}else{
				U[i+1]=(C[i]-(A[i][i+1]*U[i+2]))/A[i][i];
			}
			
		}
		error=0.0;
		for(i=1;i<I-1;i++)
		{	
			
			error = error + (U[i]-u[i]);
			
		}
		
		for(i=1;i<I-1;i++)
		{
			u[i]= U[i];
		}
		if(t==15){
			f=fopen("Y VS U crank nicholson at t=0.15 .csv","w");
			for (i=0;i<I;i++){
				y=i*dy;
				fprintf(f,"%lf,%lf\n",y,u[i]);
			}
			fclose(f);
		}else if (t==80){
			f=fopen("Y VS U crank nicholson at t=0.8 .csv","w");
			
			for (i=0;i<I;i++){
				y=i*dy;
				fprintf(f,"%lf,%lf\n",y,u[i]);
			}
			fclose(f);
		}else if (t==150){
			f=fopen("Y VS U crank nicholson at t=1.5 .csv","w");
			
			for (i=0;i<I;i++){
				y=i*dy;
				fprintf(f,"%lf,%lf\n",y,u[i]);
			}
			fclose(f);
		}			
	}
}

void ES (){
	/*here 
	eta= y/(2*sqrt(nu*t)
	eta1=h/(2*sqrt(nu*t)
	erf is erf(eta) 
	erf1 is erf((2*eta1)+eta)
	erf2 is erf((2*eta1)+eta)
	UL,h,sum,z,erfc with 1,2 represents at respective erf calcultion
	m is used for summation of all terms*/
	int i,j,I=61,a,k,T;
	double eta,eta1,nu=0.000217,h,H=0.04,dy,z,sum=0,erf,erfc,y,m,UL,U[I];
	double UL1,h1,sum1=0,z1,erf1,erfc1;
	double UL2,h2,sum2=0,z2,erf2,erfc2;
	double  t;
	printf ("Enter the time step\n");
	scanf("%lf",&t);
	dy=H/(I-1);	
	FILE * f;
	
	for(i=0;i<(I);i++)
	{
		eta1=H/(2*(pow(nu*t,0.5)));
		y=i*dy;
		eta=y/(2*pow((nu*t),0.5));
		a=0;
		h=(eta-a)/100;
		sum=0.0;	
		for(j=1;j<100;j++)
		{
			z=a+(j*h);
			sum=sum+((2/1.772453)*exp(-z*z));
		}
		erf=(h/2.0)*(((2/1.772453)*exp(-a*a))+2*sum+((2/1.772453)*exp(-eta*eta)));
		
		erfc=1-erf;
		// sum of error functions
		m = 0.0;
		for(k=1;k<20;k++){
			// 1st error function in summation
			UL1=(y/(2*pow((nu*t),0.5)))+(2*k*eta1);

			a=0;
			h1=(UL1-a)/100;
			sum1=0.0;	
			for(j=1;j<100;j++)
			{
				z1=a+(j*h1);
				sum1=sum1+((2/1.772453)*exp(-z1*z1));
			}
			erf1=(h1/2.0)*(((2/1.772453)*exp(-a*a))+2*sum1+((2/1.772453)*exp(-UL1*UL1)));
			erfc1=1-erf1;
			//2nd erfc in summation 
			UL2= -(y/(2*pow((nu*t),0.5)))+(2*k*eta1);;
			a=0;
			h2=(UL2-a)/100;	
			sum2=0.0;
			for(j=1;j<100;j++)
			{
				z2=a+(j*h2);
				
				sum2=sum2+((2/1.772453)*exp(-z2*z2));
			}
			erf2=(h2/2.0)*(((2/1.772453)*exp(-a*a))+2*sum2+((2/1.772453)*exp(-UL2*UL2)));
			
			erfc2=1-erf2;
			
			m=m+(+erfc1-erfc2);
				
		}
		
		U[i]=40*(erfc+m);
	}
	f=fopen("Y VS U Exact solution.csv","w");
	fprintf(f,"At time step  %0.2lf\n\n\n",t);
	for(i=0;i<I;i++){
		y=i*dy;
		fprintf(f,"%lf,%lf\n",y,U[i]);
	}
	fclose(f);
	
}
	

int main(){
	int i;
	printf("Please enter the number for respective scheme to get output for respective scheme\n");
	printf("1.Forward time Centered space\n");
	printf("2.Crank nicholson\n");
	printf("3.Exact solution\n");
	scanf("%d",&i);
	switch (i){
		case 1 :{
			FTCS();
			printf("Files are written for FTCS");
			break;
		}
		case 2 :{
			CNS();
			printf("Files are written for crank nicholson");
			break;
		}
		case 3 :{
			ES();
			printf("Files are written for Exact solution");
			break;
		}
		
	}
	return 0;
}
