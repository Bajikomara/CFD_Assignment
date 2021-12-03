#include<stdio.h>
#include<stdlib.h>
#include<math.h>
void SW (){
	//variable declaration
	int N=103,x,i,j;
	double dx=0.01,dt,t;
	double u[N],U_1[N],U_2[N],U_3[N],Fps_1[N],Fns_1[N],Fps_2[N],Fns_2[N],Fps_3[N],Fns_3[N],a[N];
	double Temp[N],P[N],rho[N],M[N],E[N],sum[N];
	double gma=1.4,umax,tempmax,Mmax,Pmax,lamdmax;
	/* the below are applicable to both SW and VL schemes
	N is n.o of grid pints
	x is used ti calculate iterations
	i is used to calculate to march in time
	dx is the length of smallest element
	dt is the time
	t is used to calculate and reach the required time
	Fps is values on  F+ matrix
	Fns is values on  F- matrix
	U is values in Umatrix
	SUM[N] calculates fabs(u[i])+a[i]
	rest all are self explanatory
	*/
	
	FILE * f;// for writing to a file
	//initial conditions
	for(i=0;i<N;i++)
	{
		Temp[i]=300.0;
		u[i]=0.0;
		if(i>=0 && i<51 )
	  {
  		P[i]=5.0*101325.0; // left of diapragm pressure
	  }else if(i>=51 && i<=N)
	  {
	  	P[i]=1.0*101325.0; // right of diapragm pressure
	  }
	  rho[i]=P[i]/(287*Temp[i]);
	  U_1[i]=rho[i];
	  U_2[i]=rho[i]*u[i];
	  E[i]=(P[i]/0.4)+(0.5*rho[i]*u[i]*u[i]);
	  U_3[i]=E[i];
	  a[i]=sqrt(1.4*287.0*Temp[i]); 
	  M[i]=u[i]/a[i];
	}

	t=0.0;
	x=0;
	//calculating at each grid point
	while (t<0.00075){
		x+=1;
		//calculating lamda max
		for(i=1;i<N;i++){
			sum[i]=fabs(u[i])+a[i];
			if(lamdmax < sum[i]){
				lamdmax=sum[i];
			}
		}
		dt=0.5*dx/lamdmax;
		t+=dt;
		
		for(i=1;i<N-1;i++){
			
			if(M[i]<-1.0){
				for(j=0;j<3;j++){
					Fns_1[i+j-1] =rho[i+j-1]*u[i+j-1];
					Fns_2[i+j-1] =(rho[i+j-1]*u[i+j-1]*u[i+j-1])+P[i+j-1];
					Fns_3[i+j-1] =(E[i+j-1]+P[i+j-1])*u[i+j-1];
					
					Fps_1[i+j-1]=0;
					Fps_2[i+j-1]=0;
					Fps_3[i+j-1]=0;
				}
			}else if(M[i]>-1.0 && M[i]<=0.0 ){
				for(j=0;j<3;j++){
					Fps_1[i+j-1] = ((rho[i+j-1])/(2*1.4))*(u[i+j-1]+a[i+j-1]);
					Fps_2[i+j-1] = ((rho[i+j-1])/(2*1.4))*(u[i+j-1]+a[i+j-1])*(u[i+j-1]+a[i+j-1]);
					Fps_3[i+j-1] = ((rho[i+j-1])/(2*1.4))*(u[i+j-1]+a[i+j-1])*((0.5*u[i+j-1]*u[i+j-1])+((a[i+j-1]*a[i+j-1])/0.4)+(u[i+j-1]*a[i+j-1]));
					
					Fns_1[i+j-1] =(0.2857*rho[i+j-1]*u[i+j-1])+ ((rho[i+j-1])/(2*1.4))*(u[i+j-1]-a[i+j-1]);
					Fns_2[i+j-1] =(0.2857*rho[i+j-1]*u[i+j-1]*u[i+j-1])+((rho[i+j-1])/(2*1.4))*(u[i+j-1]-a[i+j-1])*(u[i+j-1]-a[i+j-1]);
					Fns_3[i+j-1] =(0.2857*rho[i+j-1]*u[i+j-1]*u[i+j-1]*u[i+j-1]*0.5)+((rho[i+j-1])/(2*1.4))*(u[i+j-1]-a[i+j-1])*((0.5*u[i+j-1]*u[i+j-1])+((a[i+j-1]*a[i+j-1])/0.4)-(u[i+j-1]*a[i+j-1]));
				}
				
			}else if(M[i]>0.0 && M[i]<=1.0){
				for(j=0;j<3;j++){
					Fps_1[i+j-1] =((0.4/1.4)*rho[i+j-1]*u[i+j-1])+(((0.5/1.4)*rho[i+j-1])*(u[i+j-1]+a[i+j-1]));
					Fps_2[i+j-1] =((0.4/1.4)*rho[i+j-1]*u[i+j-1]*u[i+j-1])+(((0.5/1.4)*rho[i+j-1])*((u[i+j-1]+a[i+j-1])*(u[i+j-1]+a[i+j-1])));
					Fps_3[i+j-1] =((0.4/1.4)*rho[i+j-1]*u[i+j-1]*u[i+j-1]*u[i+j-1]*0.5)+(((0.5/1.4)*rho[i+j-1])*(u[i+j-1]+a[i+j-1])*((0.5*u[i+j-1]*u[i+j-1])+((1/0.4)*a[i+j-1]*a[i+j-1])+(a[i+j-1]*u[i+j-1])));
					
					Fns_1[i+j-1]=(((0.5/1.4)*rho[i+j-1])*(u[i+j-1]-a[i+j-1]));
					Fns_2[i+j-1]=(((0.5/1.4)*rho[i+j-1])*((u[i+j-1]-a[i+j-1])*(u[i+j-1]-a[i+j-1])));
					Fns_3[i+j-1]=(((0.5/1.4)*rho[i+j-1])*(u[i+j-1]-a[i+j-1])*((0.5*u[i+j-1]*u[i+j-1])+((1/0.4)*a[i+j-1]*a[i+j-1])-(a[i+j-1]*u[i+j-1])));
				}
				
			}else if(M[i]>1.0){
				for(j=0;j<3;j++){
					Fps_1[i+j-1] =rho[i+j-1]*u[i+j-1];
					Fps_2[i+j-1] =(rho[i+j-1]*u[i+j-1]*u[i+j-1])+P[i+j-1];
					Fps_3[i+j-1] =(E[i+j-1]+P[i+j-1])*u[i+j-1];
					
					Fns_1[i+j-1]=0;	
					Fns_2[i+j-1]=0;
					Fns_3[i+j-1]=0;
				}
			}
			U_1[i]=rho[i]-((dt/dx)*(Fps_1[i]-Fps_1[i-1])) - ((dt/dx)*(Fns_1[i+1]-Fns_1[i]));
			U_2[i]=(rho[i]*u[i])-((dt/dx)*(Fps_2[i]-Fps_2[i-1])) - ((dt/dx)*(Fns_2[i+1]-Fns_2[i]));
			U_3[i]=E[i]-((dt/dx)*(Fps_3[i]-Fps_3[i-1])) - ((dt/dx)*(Fns_3[i+1]-Fns_3[i]));
			
		
		}
		for(i=1;i<N;i++){
			
			rho[i]=U_1[i];
			u[i]=U_2[i]/rho[i];
			
			P[i]=0.4*(E[i]-(0.5*rho[i]*u[i]*u[i]));
			E[i]=U_3[i];
			Temp[i]=P[i]/(287*rho[i]);
			a[i]=pow((gma*287*Temp[i]),0.5);
			M[i]=u[i]/a[i];
			
		}
		
	} 
	
	for(i=1;i<N;i++){
			
		if(umax < u[i]){
			umax=u[i];
		}
	}
	
	for(i=1;i<N;i++){
			
		if(Mmax < M[i]){
			Mmax=M[i];
		}
	}
	
	for(i=1;i<N;i++){
			
		if(tempmax < Temp[i]){
			tempmax=Temp[i];
		}
	}
	
	for(i=1;i<N;i++){
			
		if(Pmax < P[i]){
			Pmax=P[i];
		}
	}
	
	
	f=fopen("steger-Warming .csv","w");
	fprintf(f,"Grid point, Veocity,Speed of sound,Temperature,Pressure,Mach Number\n");
	for(i=1;i<N;i++){
		
		fprintf(f," %d,%lf,%lf,%lf,%lf,%lf\n",i,u[i],a[i],Temp[i],P[i],M[i]);
	} 
}

void VL (){
	int N=103,x,i,j;
	double dx=0.01,dt,t;
	double u[N],U_1[N],U_2[N],U_3[N],Fps_1[N],Fns_1[N],Fps_2[N],Fns_2[N],Fps_3[N],Fns_3[N],a[N];
	double Temp[N],P[N],rho[N],M[N],E[N],sum[N];
	double gma=1.4,umax,tempmax,Mmax,Pmax,lamdamax;
	FILE * f;
	//initial conditions
	for(i=0;i<N;i++)
	{
		Temp[i]=300.0;
		u[i]=0.0;
		if(i>=0 && i<51 )
	  {
  		P[i]=5.0*101325.0;
	  }else if(i>=51 && i<=N)
	  {
	  	P[i]=1.0*101325.0;
	  }
	  rho[i]=P[i]/(287*Temp[i]);
	  U_1[i]=rho[i];
	  U_2[i]=rho[i]*u[i];
	  E[i]=(P[i]/0.4)+(0.5*rho[i]*u[i]*u[i]);
	  U_3[i]=E[i];
	  a[i]=sqrt(1.4*287.0*Temp[i]); 
	  M[i]=u[i]/a[i];
	}

	t=0.0;
	x=0;
	while (t<0.00075){
		x+=1;
		for(i=1;i<N;i++){
			sum[i]=fabs(u[i])+a[i];
			if(lamdamax < sum[i]){
				lamdamax=sum[i];
			}
		}
		dt=0.5*dx/lamdamax;
		t+=dt;
		//printf("time = %lf\n",t);
		for(i=1;i<N-1;i++){
			
			if(M[i]<-1.0){
				for(j=0;j<3;j++){
					Fns_1[i+j-1] =rho[i+j-1]*u[i+j-1];
					Fns_2[i+j-1] =(rho[i+j-1]*u[i+j-1]*u[i+j-1])+P[i+j-1];
					Fns_3[i+j-1] =(E[i+j-1]+P[i+j-1])*u[i+j-1];
					
					Fps_1[i+j-1]=0;
					Fps_2[i+j-1]=0;
					Fps_3[i+j-1]=0;
				}
			}else if(M[i]>-1.0 && M[i]<=1.0 ){
				for(j=0;j<3;j++){
					Fps_1[i+j-1] = (0.25*rho[i+j-1]*a[i+j-1])*(1+M[i+j-1])*(1+M[i+j-1]);
					Fps_2[i+j-1] = ((0.25*rho[i+j-1]*a[i+j-1])*(1+M[i+j-1])*(1+M[i+j-1]))*((2*a[i+j-1]/1.4))*((0.2*M[i+j-1])+1);
					Fps_3[i+j-1] = ((0.25*rho[i+j-1]*a[i+j-1])*(1+M[i+j-1])*(1+M[i+j-1]))*((2*a[i+j-1]*a[i+j-1]/0.96)*((0.2*M[i+j-1])+1)*((0.2*M[i+j-1])+1));
					
					Fns_1[i+j-1] =-(0.25*rho[i+j-1]*a[i+j-1])*(-1+M[i+j-1])*(-1+M[i+j-1]);
					Fns_2[i+j-1] =-((0.25*rho[i+j-1]*a[i+j-1])*(-1+M[i+j-1])*(-1+M[i+j-1]))*((2*a[i+j-1]/1.4))*((0.2*M[i+j-1])-1);
					Fns_3[i+j-1] =-((0.25*rho[i+j-1]*a[i+j-1])*(-1+M[i+j-1])*(-1+M[i+j-1]))*((2*a[i+j-1]*a[i+j-1]/0.96)*((0.2*M[i+j-1])-1)*((0.2*M[i+j-1])-1));
				}
				
			
				
			}else if(M[i]>1.0){
				for(j=0;j<3;j++){
					Fps_1[i+j-1] =rho[i+j-1]*u[i+j-1];
					Fps_2[i+j-1] =(rho[i+j-1]*u[i+j-1]*u[i+j-1])+P[i+j-1];
					Fps_3[i+j-1] =(E[i+j-1]+P[i+j-1])*u[i+j-1];
					
					Fns_1[i+j-1]=0;	
					Fns_2[i+j-1]=0;
					Fns_3[i+j-1]=0;
				}
			}
			U_1[i]=rho[i]-((dt/dx)*(Fps_1[i]-Fps_1[i-1])) - ((dt/dx)*(Fns_1[i+1]-Fns_1[i]));
			U_2[i]=(rho[i]*u[i])-((dt/dx)*(Fps_2[i]-Fps_2[i-1])) - ((dt/dx)*(Fns_2[i+1]-Fns_2[i]));
			U_3[i]=E[i]-((dt/dx)*(Fps_3[i]-Fps_3[i-1])) - ((dt/dx)*(Fns_3[i+1]-Fns_3[i]));
			
		
		}
		for(i=1;i<N;i++){
			
			rho[i]=U_1[i];
			u[i]=U_2[i]/rho[i];
			//P[i]=0.4*(U_3[i]-(0.5*(U_2[i]*U_2[i])/U_1[i]));
			
			P[i]=0.4*(E[i]-(0.5*rho[i]*u[i]*u[i]));
			E[i]=U_3[i];
			Temp[i]=P[i]/(287*rho[i]);
			a[i]=pow((gma*287*Temp[i]),0.5);
			M[i]=u[i]/a[i];
			
		}
		
	} 
	
	for(i=1;i<N;i++){
			
		if(umax < u[i]){
			umax=u[i];
		}
	}
	//umax=u[1];
	for(i=1;i<N;i++){
			
		if(Mmax < M[i]){
			Mmax=M[i];
		}
	}
	
	for(i=1;i<N;i++){
			
		if(tempmax < Temp[i]){
			tempmax=Temp[i];
		}
	}
	
	for(i=1;i<N;i++){
			
		if(Pmax < P[i]){
			Pmax=P[i];
		}
	}
	
	

	f=fopen("Vanleer .csv","w");
	fprintf(f,"Grid point, Veocity,Speed of sound,Temperature,Pressure,Mach Number\n");
	for(i=1;i<N;i++){
		
		fprintf(f," %d,%lf,%lf,%lf,%lf,%lf\n",i,u[i],a[i],Temp[i],P[i],M[i]);
	} 
}

int main(){
	int i;
	printf("Enter the flux scheme to compute:\n");
	printf("1.steger-Warming\n");
	printf("2.Lan Veer\n");
	scanf("%d",&i);
	switch (i){
		case 1 :{
			SW();
			printf("Files are written for steger-Warming");
			break;
		}
		case 2 :{
			VL();
			printf("Files are written for Lan Veer");
			break;
		}
	}
    return 0;
}




