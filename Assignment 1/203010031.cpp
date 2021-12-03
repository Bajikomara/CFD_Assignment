#include<iostream>
#include<math.h>
#include<stdlib.h>
#include<fstream>
#include<new>
#include<string>
using namespace std;
int main(){
	ofstream myfile;
	myfile.open("20301031 .csv");
	int i,t,a=1,type,ini_sol,ts;
	double *U,*u,nu,Nu,*Lu,*Bu,sigm=0.4;
	string scheme;
	U=new double[500];
	Lu=new double[500];
	Bu=new double[500];
	u=new double[500];
	double dx=0.01,dt;
	cout<<"Enter the value of nu "<<endl;
	cin>>nu;
	cout<<"Enter the time steps"<<endl;
	cin>>ts;
	for(ini_sol=1;ini_sol<=5;ini_sol++){
		cout<<"Enter the initial solution type "<<endl;
		cin>>ini_sol;
		if(ini_sol==1){
			myfile<<ini_sol<<"st initial solution "<<endl;
			myfile<<"nu is "<<nu<<endl;
			for(type=1;type<=6;type++){
				cout<<"Enter the  Finite Difference Schemes"<<endl;
				cin>>type;
				if(type==1){
					scheme = "Forward Time Forward Space (FTFS) scheme";
					myfile<<"type is "<<scheme<<endl;	
					for(t=0;t<=ts;t++){
						dt=t*nu*dx;
						cout<<"time step is:"<<t<<endl;
						myfile<<"time step is "<<t<<endl;
						for(i=1;i<=100;i+=1){
							if(i<=20){
            					u[i]=1;
        			
        					}else{
            					u[i]=0;
        					}
        					Nu=a*(dt/dx);
        					U[i]=u[i]-((dt/dx)*(u[i+1]-u[i]));;
        					cout<<U[i]<<endl;
        					myfile<<i*0.01<<","<<U[i]<<endl;
        					//cout<<"At time step "<<t<<" and "<<i<<" iteration is : "<<U[i]<<endl;	
						}		
					}
				}else if(type==2){
					scheme = "Forward Time Central Space (FTCS) scheme";
					myfile<<"type is  "<<scheme<<endl;	
					for(t=0;t<=ts;t++){
						dt=t*nu*dx;
						cout<<"time step is:"<<t<<endl;
						myfile<<"time step is "<<t<<endl;
						for(i=1;i<=100;i+=1){
							if(i<=20){
            					u[i]=1;
        			
        					}else{
            					u[i]=0;
        					}
        					Nu=a*(dt/dx);
        					U[i]=u[i]-((0.5*Nu)*(u[i+1]-u[i-1]));
        					cout<<U[i]<<endl;
        					myfile<<i*0.01<<","<<U[i]<<endl;
        					//cout<<"At time step "<<t<<" and "<<i<<" iteration is : "<<U[i]<<endl;	
						}
					}	
				}else if(type==3){
					scheme = "Forward Time Backward Space (FTBS) scheme";
					myfile<<"type is  "<<scheme<<endl;	
					for(t=0;t<=ts;t++){
						dt=t*nu*dx;
						cout<<"time step is:"<<t<<endl;
						myfile<<"time step is "<<t<<endl;
						for(i=1;i<=100;i+=1){
							if(i<=20){
            					u[i]=1;
        			
        					}else{
            					u[i]=0;
        					}
        					Nu=a*(dt/dx);
        					U[i]=u[i]-((0.5*Nu)*(u[i]-u[i-1]));
        					cout<<U[i]<<endl;
        					myfile<<i*0.01<<","<<U[i]<<endl;
        					//cout<<"At time step "<<t<<" and "<<i<<" iteration is : "<<U[i]<<endl;	
						}
					}
				}else if(type==4){
					scheme = "Lax-Wendroff (LW) scheme";
					myfile<<"type is  "<<scheme<<endl;	
					for(t=0;t<=ts;t++){
						dt=t*nu*dx;
						cout<<"time step is:"<<t<<endl;
						myfile<<"time step is "<<t<<endl;
						for(i=1;i<=100;i+=1){
							if(i<=20){
            					u[i]=1;
        			
        					}else{
            					u[i]=0;
        					}
        					Nu=a*(dt/dx);
        					Lu[i]=u[i]-((0.5*Nu)*(u[i+1]-u[i-1]))-((0.5*Nu*Nu)*(u[i+1]+u[i-1]-(2*u[i])));
        					cout<<U[i]<<endl;
        					myfile<<i*0.01<<","<<Lu[i]<<endl;
        				//cout<<"At time step "<<t<<" and "<<i<<" iteration is : "<<U[i]<<endl;	
						}
					}
				}else if(type==5){
					scheme = "Beam-Warming (BW) scheme";
					myfile<<"type is  "<<scheme<<endl;	
					for(t=0;t<=ts;t++){
						dt=t*nu*dx;
						cout<<"time step is:"<<t<<endl;
						myfile<<"time step is "<<t<<endl;
						for(i=1;i<=100;i+=1){
							if(i<=20){
            					u[i]=1;
        			
        					}else{
            					u[i]=0;
        					}
        					Nu=a*(dt/dx);
        					Bu[i]=u[i]-((0.5*Nu)*(3*u[i]-4*u[i-1]+u[i-2]))-((0.5*Nu*Nu)*(u[i+1]+u[i-2]-(2*u[i-1])));
        					cout<<U[i]<<endl;
        					myfile<<i*0.01<<","<<Bu[i]<<endl;
        					//cout<<"At time step "<<t<<" and "<<i<<" iteration is : "<<U[i]<<endl;	
						}
					}
			
				}else if(type==6){
					scheme = "Fromm (FR) scheme";
					myfile<<"type is  "<<scheme<<endl;	
					for(t=0;t<=ts;t++){
						dt=t*nu*dx;
						cout<<"time step is:"<<t<<endl;
						myfile<<"time step is "<<t<<endl;
						for(i=1;i<=100;i+=1){
							if(i<=20){
            					u[i]=1;
        			
        					}else{
            					u[i]=0;
        					}
        					Nu=a*(dt/dx);
        					U[i]=0.5*(Lu[i]+Bu[i]);
        					cout<<U[i]<<endl;
        					myfile<<i*0.01<<","<<U[i]<<endl;
        				//cout<<"At time step "<<t<<" and "<<i<<" iteration is : "<<U[i]<<endl;	
						}
					}
				}
			}
		}else if(ini_sol==2){
			myfile<<ini_sol<<"st initial solution "<<endl;
			
			myfile<<"nu is "<<nu<<endl;
			for(type=1;type<=6;type++){
				cout<<"Enter the  Finite Difference Schemes"<<endl;
				cin>>type;
				if(type==1){
					scheme = "Forward Time Forward Space (FTFS) scheme";
					myfile<<"type is "<<scheme<<endl;	
					for(t=0;t<=ts;t++){
						dt=t*nu*dx;
						cout<<"time step is:"<<t<<endl;
						myfile<<"time step is "<<t<<endl;
						for(i=1;i<=100;i+=1){
							if(i>=0 && i<5 ){
            					u[i]=0;
        			
        					}else if (i>=5 && i<35){
            					u[i]=sin((4*M_PI*(((i*0.01)-0.05)/0.3)));
        					}else{
        						u[i]=0;
							}
        					Nu=a*(dt/dx);
        					U[i]=u[i]-((dt/dx)*(u[i+1]-u[i]));;
        					cout<<U[i]<<endl;
        					myfile<<i*0.01<<","<<U[i]<<endl;
        					//cout<<"At time step "<<t<<" and "<<i<<" iteration is : "<<U[i]<<endl;	
						}		
					}
				}else if(type==2){
					scheme = "Forward Time Central Space (FTCS) scheme";
					myfile<<"type is  "<<scheme<<endl;	
					for(t=0;t<=ts;t++){
						dt=t*nu*dx;
						cout<<"time step is:"<<t<<endl;
						myfile<<"time step is "<<t<<endl;
						for(i=1;i<=100;i+=1){
							if(i>=0 && i<5 ){
            					u[i]=0;
        			
        					}else if (i>=5 && i<35){
            					u[i]=sin((4*M_PI*(((i*0.01)-0.05)/0.3)));
        					}else{
        						u[i]=0;
							}
        					Nu=a*(dt/dx);
        					U[i]=u[i]-((0.5*Nu)*(u[i+1]-u[i-1]));
        					cout<<U[i]<<endl;
        					myfile<<i*0.01<<","<<U[i]<<endl;
        					//cout<<"At time step "<<t<<" and "<<i<<" iteration is : "<<U[i]<<endl;	
						}
					}	
				}else if(type==3){
					scheme = "Forward Time Backward Space (FTBS) scheme";
					myfile<<"type is  "<<scheme<<endl;	
					for(t=0;t<=ts;t++){
						dt=t*nu*dx;
						cout<<"time step is:"<<t<<endl;
						myfile<<"time step is "<<t<<endl;
						for(i=1;i<=100;i+=1){
							if(i>=0 && i<5 ){
            					u[i]=0;
        			
        					}else if (i>=5 && i<35){
            					u[i]=sin((4*M_PI*(((i*0.01)-0.05)/0.3)));
        					}else{
        						u[i]=0;
							}
        					Nu=a*(dt/dx);
        					U[i]=u[i]-((0.5*Nu)*(u[i]-u[i-1]));
        					cout<<U[i]<<endl;
        					myfile<<i*0.01<<","<<U[i]<<endl;
        					//cout<<"At time step "<<t<<" and "<<i<<" iteration is : "<<U[i]<<endl;	
						}
					}
				}else if(type==4){
					scheme = "Lax-Wendroff (LW) scheme";
					myfile<<"type is  "<<scheme<<endl;	
					for(t=0;t<=ts;t++){
						dt=t*nu*dx;
						cout<<"time step is:"<<t<<endl;
						myfile<<"time step is "<<t<<endl;
						for(i=1;i<=100;i+=1){
							if(i>=0 && i<5 ){
            					u[i]=0;
        			
        					}else if (i>=5 && i<35){
            					u[i]=sin((4*M_PI*(((i*0.01)-0.05)/0.3)));
        					}else{
        						u[i]=0;
							}
        					Nu=a*(dt/dx);
        					Lu[i]=u[i]-((0.5*Nu)*(u[i+1]-u[i-1]))-((0.5*Nu*Nu)*(u[i+1]+u[i-1]-(2*u[i])));
        					cout<<U[i]<<endl;
        					myfile<<i*0.01<<","<<Lu[i]<<endl;
        				//cout<<"At time step "<<t<<" and "<<i<<" iteration is : "<<U[i]<<endl;	
						}
					}
				}else if(type==5){
					scheme = "Beam-Warming (BW) scheme";
					myfile<<"type is  "<<scheme<<endl;	
					for(t=0;t<=ts;t++){
						dt=t*nu*dx;
						cout<<"time step is:"<<t<<endl;
						myfile<<"time step is "<<t<<endl;
						for(i=1;i<=100;i+=1){
						if(i>=0 && i<5 ){
            					u[i]=0;
        			
        					}else if (i>=5 && i<35){
            					u[i]=sin((4*M_PI*(((i*0.01)-0.05)/0.3)));
        					}else{
        						u[i]=0;
							}
        					Nu=a*(dt/dx);
        					Bu[i]=u[i]-((0.5*Nu)*(3*u[i]-4*u[i-1]+u[i-2]))-((0.5*Nu*Nu)*(u[i+1]+u[i-2]-(2*u[i-1])));
        					cout<<U[i]<<endl;
        					myfile<<i*0.01<<","<<Bu[i]<<endl;
        					//cout<<"At time step "<<t<<" and "<<i<<" iteration is : "<<U[i]<<endl;	
						}
					}
			
				}else if(type==6){
					scheme = "Fromm (FR) scheme";
					myfile<<"type is  "<<scheme<<endl;	
					for(t=0;t<=ts;t++){
						dt=t*nu*dx;
						cout<<"time step is:"<<t<<endl;
						myfile<<"time step is "<<t<<endl;
						for(i=1;i<=100;i+=1){
							if(i>=0 && i<5 ){
            					u[i]=0;
        			
        					}else if (i>=5 && i<35){
            					u[i]=sin((4*M_PI*(((i*0.01)-0.05)/0.3)));
        					}else{
        						u[i]=0;
							}
        					Nu=a*(dt/dx);
        					U[i]=0.5*(Lu[i]+Bu[i]);
        					cout<<U[i]<<endl;
        					myfile<<i*0.01<<","<<U[i]<<endl;
        				//cout<<"At time step "<<t<<" and "<<i<<" iteration is : "<<U[i]<<endl;	
						}
					}
				}
			}

			
		}else if (ini_sol==3){
			myfile<<ini_sol<<"st initial solution "<<endl;
			
			myfile<<"nu is "<<nu<<endl;
			for(type=1;type<=6;type++){
				cout<<"Enter the  Finite Difference Schemes"<<endl;
				cin>>type;
				if(type==1){
					scheme = "Forward Time Forward Space (FTFS) scheme";
					myfile<<"type is "<<scheme<<endl;	
					for(t=0;t<=ts;t++){
						dt=t*nu*dx;
						cout<<"time step is:"<<t<<endl;
						myfile<<"time step is "<<t<<endl;
						for(i=1;i<=100;i+=1){
							if(i>=0 && i<5 ){
            					u[i]=0;
        			
        					}else if (i>=5 && i<35){
            					u[i]=sin((8*M_PI*(((i*0.01)-0.05)/0.3)));
        					}else{
        						u[i]=0;
							}
        					Nu=a*(dt/dx);
        					U[i]=u[i]-((dt/dx)*(u[i+1]-u[i]));;
        					cout<<U[i]<<endl;
        					myfile<<i*0.01<<","<<U[i]<<endl;
        					//cout<<"At time step "<<t<<" and "<<i<<" iteration is : "<<U[i]<<endl;	
						}		
					}
				}else if(type==2){
					scheme = "Forward Time Central Space (FTCS) scheme";
					myfile<<"type is  "<<scheme<<endl;	
					for(t=0;t<=ts;t++){
						dt=t*nu*dx;
						cout<<"time step is:"<<t<<endl;
						myfile<<"time step is "<<t<<endl;
						for(i=1;i<=100;i+=1){
							if(i>=0 && i<5 ){
            					u[i]=0;
        			
        					}else if (i>=5 && i<35){
            					u[i]=sin((8*M_PI*(((i*0.01)-0.05)/0.3)));
        					}else{
        						u[i]=0;
							}
        					Nu=a*(dt/dx);
        					U[i]=u[i]-((0.5*Nu)*(u[i+1]-u[i-1]));
        					cout<<U[i]<<endl;
        					myfile<<i*0.01<<","<<U[i]<<endl;
        					//cout<<"At time step "<<t<<" and "<<i<<" iteration is : "<<U[i]<<endl;	
						}
					}	
				}else if(type==3){
					scheme = "Forward Time Backward Space (FTBS) scheme";
					myfile<<"type is  "<<scheme<<endl;	
					for(t=0;t<=ts;t++){
						dt=t*nu*dx;
						cout<<"time step is:"<<t<<endl;
						myfile<<"time step is "<<t<<endl;
						for(i=1;i<=100;i+=1){
							if(i>=0 && i<5 ){
            					u[i]=0;
        			
        					}else if (i>=5 && i<35){
            					u[i]=sin((8*M_PI*(((i*0.01)-0.05)/0.3)));
        					}else{
        						u[i]=0;
							}
        					Nu=a*(dt/dx);
        					U[i]=u[i]-((0.5*Nu)*(u[i]-u[i-1]));
        					cout<<U[i]<<endl;
        					myfile<<i*0.01<<","<<U[i]<<endl;
        					//cout<<"At time step "<<t<<" and "<<i<<" iteration is : "<<U[i]<<endl;	
						}
					}
				}else if(type==4){
					scheme = "Lax-Wendroff (LW) scheme";
					myfile<<"type is  "<<scheme<<endl;	
					for(t=0;t<=ts;t++){
						dt=t*nu*dx;
						cout<<"time step is:"<<t<<endl;
						myfile<<"time step is "<<t<<endl;
						for(i=1;i<=100;i+=1){
							if(i>=0 && i<5 ){
            					u[i]=0;
        			
        					}else if (i>=5 && i<35){
            					u[i]=sin((8*M_PI*(((i*0.01)-0.05)/0.3)));
        					}else{
        						u[i]=0;
							}
        					Nu=a*(dt/dx);
        					Lu[i]=u[i]-((0.5*Nu)*(u[i+1]-u[i-1]))-((0.5*Nu*Nu)*(u[i+1]+u[i-1]-(2*u[i])));
        					cout<<U[i]<<endl;
        					myfile<<i*0.01<<","<<Lu[i]<<endl;
        				//cout<<"At time step "<<t<<" and "<<i<<" iteration is : "<<U[i]<<endl;	
						}
					}
				}else if(type==5){
					scheme = "Beam-Warming (BW) scheme";
					myfile<<"type is  "<<scheme<<endl;	
					for(t=0;t<=ts;t++){
						dt=t*nu*dx;
						cout<<"time step is:"<<t<<endl;
						myfile<<"time step is "<<t<<endl;
						for(i=1;i<=100;i+=1){
						if(i>=0 && i<5 ){
            					u[i]=0;
        			
        					}else if (i>=5 && i<35){
            					u[i]=sin((8*M_PI*(((i*0.01)-0.05)/0.3)));
        					}else{
        						u[i]=0;
							}
        					Nu=a*(dt/dx);
        					Bu[i]=u[i]-((0.5*Nu)*(3*u[i]-4*u[i-1]+u[i-2]))-((0.5*Nu*Nu)*(u[i+1]+u[i-2]-(2*u[i-1])));
        					cout<<U[i]<<endl;
        					myfile<<i*0.01<<","<<Bu[i]<<endl;
        					//cout<<"At time step "<<t<<" and "<<i<<" iteration is : "<<U[i]<<endl;	
						}
					}
			
				}else if(type==6){
					scheme = "Fromm (FR) scheme";
					myfile<<"type is  "<<scheme<<endl;	
					for(t=0;t<=ts;t++){
						dt=t*nu*dx;
						cout<<"time step is:"<<t<<endl;
						myfile<<"time step is "<<t<<endl;
						for(i=1;i<=100;i+=1){
							if(i>=0 && i<5 ){
            					u[i]=0;
        			
        					}else if (i>=5 && i<35){
            					u[i]=sin((8*M_PI*(((i*0.01)-0.05)/0.3)));
        					}else{
        						u[i]=0;
							}
        					Nu=a*(dt/dx);
        					U[i]=0.5*(Lu[i]+Bu[i]);
        					cout<<U[i]<<endl;
        					myfile<<i*0.01<<","<<U[i]<<endl;
        				//cout<<"At time step "<<t<<" and "<<i<<" iteration is : "<<U[i]<<endl;	
						}
					}
				}
			}
		
			
		}else if (ini_sol==4){
			myfile<<ini_sol<<"st initial solution "<<endl;
			
			myfile<<"nu is "<<nu<<endl;
			for(type=1;type<=6;type++){
				cout<<"Enter the  Finite Difference Schemes"<<endl;
				cin>>type;
				if(type==1){
					scheme = "Forward Time Forward Space (FTFS) scheme";
					myfile<<"type is "<<scheme<<endl;	
					for(t=0;t<=ts;t++){
						dt=t*nu*dx;
						cout<<"time step is:"<<t<<endl;
						myfile<<"time step is "<<t<<endl;
						for(i=1;i<=100;i+=1){
							if(i>=0 && i<5 ){
            					u[i]=0;
        			
        					}else if (i>=5 && i<35){
            					u[i]=sin((12*M_PI*(((i*0.01)-0.05)/0.3)));
        					}else{
        						u[i]=0;
							}
        					Nu=a*(dt/dx);
        					U[i]=u[i]-((dt/dx)*(u[i+1]-u[i]));;
        					cout<<U[i]<<endl;
        					myfile<<i*0.01<<","<<U[i]<<endl;
        					//cout<<"At time step "<<t<<" and "<<i<<" iteration is : "<<U[i]<<endl;	
						}		
					}
				}else if(type==2){
					scheme = "Forward Time Central Space (FTCS) scheme";
					myfile<<"type is  "<<scheme<<endl;	
					for(t=0;t<=ts;t++){
						dt=t*nu*dx;
						cout<<"time step is:"<<t<<endl;
						myfile<<"time step is "<<t<<endl;
						for(i=1;i<=100;i+=1){
							if(i>=0 && i<5 ){
            					u[i]=0;
        			
        					}else if (i>=5 && i<35){
            					u[i]=sin((12*M_PI*(((i*0.01)-0.05)/0.3)));
        					}else{
        						u[i]=0;
							}
        					Nu=a*(dt/dx);
        					U[i]=u[i]-((0.5*Nu)*(u[i+1]-u[i-1]));
        					cout<<U[i]<<endl;
        					myfile<<i*0.01<<","<<U[i]<<endl;
        					//cout<<"At time step "<<t<<" and "<<i<<" iteration is : "<<U[i]<<endl;	
						}
					}	
				}else if(type==3){
					scheme = "Forward Time Backward Space (FTBS) scheme";
					myfile<<"type is  "<<scheme<<endl;	
					for(t=0;t<=ts;t++){
						dt=t*nu*dx;
						cout<<"time step is:"<<t<<endl;
						myfile<<"time step is "<<t<<endl;
						for(i=1;i<=100;i+=1){
							if(i>=0 && i<5 ){
            					u[i]=0;
        			
        					}else if (i>=5 && i<35){
            					u[i]=sin((12*M_PI*(((i*0.01)-0.05)/0.3)));
        					}else{
        						u[i]=0;
							}
        					Nu=a*(dt/dx);
        					U[i]=u[i]-((0.5*Nu)*(u[i]-u[i-1]));
        					cout<<U[i]<<endl;
        					myfile<<i*0.01<<","<<U[i]<<endl;
        					//cout<<"At time step "<<t<<" and "<<i<<" iteration is : "<<U[i]<<endl;	
						}
					}
				}else if(type==4){
					scheme = "Lax-Wendroff (LW) scheme";
					myfile<<"type is  "<<scheme<<endl;	
					for(t=0;t<=ts;t++){
						dt=t*nu*dx;
						cout<<"time step is:"<<t<<endl;
						myfile<<"time step is "<<t<<endl;
						for(i=1;i<=100;i+=1){
							if(i>=0 && i<5 ){
            					u[i]=0;
        			
        					}else if (i>=5 && i<35){
            					u[i]=sin((12*M_PI*(((i*0.01)-0.05)/0.3)));
        					}else{
        						u[i]=0;
							}
        					Nu=a*(dt/dx);
        					Lu[i]=u[i]-((0.5*Nu)*(u[i+1]-u[i-1]))-((0.5*Nu*Nu)*(u[i+1]+u[i-1]-(2*u[i])));
        					cout<<U[i]<<endl;
        					myfile<<i*0.01<<","<<Lu[i]<<endl;
        				//cout<<"At time step "<<t<<" and "<<i<<" iteration is : "<<U[i]<<endl;	
						}
					}
				}else if(type==5){
					scheme = "Beam-Warming (BW) scheme";
					myfile<<"type is  "<<scheme<<endl;	
					for(t=0;t<=ts;t++){
						dt=t*nu*dx;
						cout<<"time step is:"<<t<<endl;
						myfile<<"time step is "<<t<<endl;
						for(i=1;i<=100;i+=1){
						if(i>=0 && i<5 ){
            					u[i]=0;
        			
        					}else if (i>=5 && i<35){
            					u[i]=sin((12*M_PI*(((i*0.01)-0.05)/0.3)));
        					}else{
        						u[i]=0;
							}
        					Nu=a*(dt/dx);
        					Bu[i]=u[i]-((0.5*Nu)*(3*u[i]-4*u[i-1]+u[i-2]))-((0.5*Nu*Nu)*(u[i+1]+u[i-2]-(2*u[i-1])));
        					cout<<U[i]<<endl;
        					myfile<<i*0.01<<","<<Bu[i]<<endl;
        					//cout<<"At time step "<<t<<" and "<<i<<" iteration is : "<<U[i]<<endl;	
						}
					}
			
				}else if(type==6){
					scheme = "Fromm (FR) scheme";
					myfile<<"type is  "<<scheme<<endl;	
					for(t=0;t<=ts;t++){
						dt=t*nu*dx;
						cout<<"time step is:"<<t<<endl;
						myfile<<"time step is "<<t<<endl;
						for(i=1;i<=100;i+=1){
							if(i>=0 && i<5 ){
            					u[i]=0;
        			
        					}else if (i>=5 && i<35){
            					u[i]=sin((12*M_PI*(((i*0.01)-0.05)/0.3)));
        					}else{
        						u[i]=0;
							}
        					Nu=a*(dt/dx);
        					U[i]=0.5*(Lu[i]+Bu[i]);
        					cout<<U[i]<<endl;
        					myfile<<i*0.01<<","<<U[i]<<endl;
        				//cout<<"At time step "<<t<<" and "<<i<<" iteration is : "<<U[i]<<endl;	
						}
					}
				}
			}
		
		}else if (ini_sol==5){
			myfile<<ini_sol<<"st initial solution "<<endl;
			
			myfile<<"nu is "<<nu<<endl;
			for(type=1;type<=6;type++){
				cout<<"Enter the  Finite Difference Schemes"<<endl;
				cin>>type;
				if(type==1){
					scheme = "Forward Time Forward Space (FTFS) scheme";
					myfile<<"type is "<<scheme<<endl;	
					for(t=0;t<=ts;t++){
						dt=t*nu*dx;
						cout<<"time step is:"<<t<<endl;
						myfile<<"time step is "<<t<<endl;
						for(i=1;i<=100;i+=1){
							if(i>=0 && i<5 ){
            					u[i]=0;
        			
        					}else if (i>=5 && i<35){
            					u[i]=exp(-50*(((i*0.01-0.2)*(i*0.01-0.2))/(sigm*sigm)));
        					}else{
        						u[i]=0;
							}
        					Nu=a*(dt/dx);
        					U[i]=u[i]-((dt/dx)*(u[i+1]-u[i]));;
        					cout<<U[i]<<endl;
        					myfile<<i*0.01<<","<<U[i]<<endl;
        					//cout<<"At time step "<<t<<" and "<<i<<" iteration is : "<<U[i]<<endl;	
						}		
					}
				}else if(type==2){
					scheme = "Forward Time Central Space (FTCS) scheme";
					myfile<<"type is  "<<scheme<<endl;	
					for(t=0;t<=ts;t++){
						dt=t*nu*dx;
						cout<<"time step is:"<<t<<endl;
						myfile<<"time step is "<<t<<endl;
						for(i=1;i<=100;i+=1){
							if(i>=0 && i<5 ){
            					u[i]=0;
        			
        					}else if (i>=5 && i<35){
            					u[i]=exp(-50*(((i*0.01-0.2)*(i*0.01-0.2))/(sigm*sigm)));
        					}else{
        						u[i]=0;
							}
        					Nu=a*(dt/dx);
        					U[i]=u[i]-((0.5*Nu)*(u[i+1]-u[i-1]));
        					cout<<U[i]<<endl;
        					myfile<<i*0.01<<","<<U[i]<<endl;
        					//cout<<"At time step "<<t<<" and "<<i<<" iteration is : "<<U[i]<<endl;	
						}
					}	
				}else if(type==3){
					scheme = "Forward Time Backward Space (FTBS) scheme";
					myfile<<"type is  "<<scheme<<endl;	
					for(t=0;t<=ts;t++){
						dt=t*nu*dx;
						cout<<"time step is:"<<t<<endl;
						myfile<<"time step is "<<t<<endl;
						for(i=1;i<=100;i+=1){
							if(i>=0 && i<5 ){
            					u[i]=0;
        			
        					}else if (i>=5 && i<35){
            					u[i]=exp(-50*(((i*0.01-0.2)*(i*0.01-0.2))/(sigm*sigm)));
        					}else{
        						u[i]=0;
							}
        					Nu=a*(dt/dx);
        					U[i]=u[i]-((0.5*Nu)*(u[i]-u[i-1]));
        					cout<<U[i]<<endl;
        					myfile<<i*0.01<<","<<U[i]<<endl;
        					//cout<<"At time step "<<t<<" and "<<i<<" iteration is : "<<U[i]<<endl;	
						}
					}
				}else if(type==4){
					scheme = "Lax-Wendroff (LW) scheme";
					myfile<<"type is  "<<scheme<<endl;	
					for(t=0;t<=ts;t++){
						dt=t*nu*dx;
						cout<<"time step is:"<<t<<endl;
						myfile<<"time step is "<<t<<endl;
						for(i=1;i<=100;i+=1){
							if(i>=0 && i<5 ){
            					u[i]=0;
        			
        					}else if (i>=5 && i<35){
            					u[i]=exp(-50*(((i*0.01-0.2)*(i*0.01-0.2))/(sigm*sigm)));
        					}else{
        						u[i]=0;
							}
        					Nu=a*(dt/dx);
        					Lu[i]=u[i]-((0.5*Nu)*(u[i+1]-u[i-1]))-((0.5*Nu*Nu)*(u[i+1]+u[i-1]-(2*u[i])));
        					cout<<U[i]<<endl;
        					myfile<<i*0.01<<","<<Lu[i]<<endl;
        				//cout<<"At time step "<<t<<" and "<<i<<" iteration is : "<<U[i]<<endl;	
						}
					}
				}else if(type==5){
					scheme = "Beam-Warming (BW) scheme";
					myfile<<"type is  "<<scheme<<endl;	
					for(t=0;t<=ts;t++){
						dt=t*nu*dx;
						cout<<"time step is:"<<t<<endl;
						myfile<<"time step is "<<t<<endl;
						for(i=1;i<=100;i+=1){
						if(i>=0 && i<5 ){
            					u[i]=0;
        			
        					}else if (i>=5 && i<35){
            					u[i]=exp(-50*(((i*0.01-0.2)*(i*0.01-0.2))/(sigm*sigm)));
        					}else{
        						u[i]=0;
							}
        					Nu=a*(dt/dx);
        					Bu[i]=u[i]-((0.5*Nu)*(3*u[i]-4*u[i-1]+u[i-2]))-((0.5*Nu*Nu)*(u[i+1]+u[i-2]-(2*u[i-1])));
        					cout<<U[i]<<endl;
        					myfile<<i*0.01<<","<<Bu[i]<<endl;
        					//cout<<"At time step "<<t<<" and "<<i<<" iteration is : "<<U[i]<<endl;	
						}
					}
			
				}else if(type==6){
					scheme = "Fromm (FR) scheme";
					myfile<<"type is  "<<scheme<<endl;	
					for(t=0;t<=ts;t++){
						dt=t*nu*dx;
						cout<<"time step is:"<<t<<endl;
						myfile<<"time step is "<<t<<endl;
						for(i=1;i<=100;i+=1){
							if(i>=0 && i<5 ){
            					u[i]=0;
        			
        					}else if (i>=5 && i<35){
            					u[i]=exp(-50*(((i*0.01-0.2)*(i*0.01-0.2))/(sigm*sigm)));
        					}else{
        						u[i]=0;
							}
        					Nu=a*(dt/dx);
        					U[i]=0.5*(Lu[i]+Bu[i]);
        					cout<<U[i]<<endl;
        					myfile<<i*0.01<<","<<U[i]<<endl;
        				//cout<<"At time step "<<t<<" and "<<i<<" iteration is : "<<U[i]<<endl;	
						}
					}
				}
			}
		
		}
	}
	return 0;
}
	
