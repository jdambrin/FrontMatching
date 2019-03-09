#include "mex.h"
#include <math.h>
#include <functional>
#include <queue>
#include <vector>
#include <iostream>
#include <limits>
#define inf INFINITY

typedef std::pair<int, int> intPair;

template<typename T> class custom_priority_queue : public std::priority_queue<T, std::vector<T> >
{
public:
    
    bool remove(const T& value) {
        auto it = std::find(this->c.begin(), this->c.end(), value);
        if (it != this->c.end()) {
            this->c.erase(it);
            std::make_heap(this->c.begin(), this->c.end(), this->comp);
            return true;
        }
        else {
            return false;
        }
    }
} ;

struct node {
    int i,j;
    double val;
    bool operator<(node other) const
    {
        return val > other.val;
    }
    
    bool operator==(node other) const
    {
        return (val == other.val) && (i==other.i) && (j==other.j);
    }
    
} ;

typedef struct node node;

double ff(double sig, double dx, double dy, double M_H, double M_V, double Mc, double tau, double V_H, double V_V, double theta){
    return (sig*(pow(dx,2)*theta*M_V+pow(dy,2)*(theta-1)*M_H)+(V_H-V_V)*tau*pow(Mc,2))/(tau*pow(Mc,2));
}

double ffprime(double sig, double dx, double dy, double M_H, double M_V, double Mc, double tau, double V_H, double V_V, double theta){
    return sig*((dx*dx+dy*dy)/(tau*Mc)-pow(dx*dx*theta+dy*dy*(theta-1),2)/(tau*tau*tau*Mc)-2*(dx*dx*theta+dy*dy*(theta-1))*(M_H-M_V)/(tau*Mc*Mc)+2*pow(M_H-M_V,2)*tau/(Mc*Mc*Mc));
}

void newtonEikonal(double* V, double* M, double V_H, double V_V, double M_H, double M_V, double dx, double dy, double sig){
    double err=1;
    double theta=0.5;
    double theta_o;
    double tau;
    double Mc;
    double f;
    double fprime;
    int k=0;
    //   std::cout << "in" << std::endl;
    while(err>1e-10){
        k++;
        theta_o=theta;
        tau=sqrt(pow(dx*theta,2)+pow(dy*(1-theta),2));
        Mc=theta*M_H+(1-theta)*M_V;
        
        if(theta<0){
            fprime=ffprime(sig,dx,dy,M_H,M_V,M_V,dy,V_H,V_V,0);
            f=ff(sig,dx,dy,M_H,M_V,M_V,dy,V_H,V_V,0)+theta*fprime;
        }
        else if(theta>1){
            fprime=ffprime(sig,dx,dy,M_H,M_V,M_H,dx,V_H,V_V,1);
            f=ff(sig,dx,dy,M_H,M_V,M_H,dx,V_H,V_V,1)+(theta-1)*fprime;
        }
        else{
            f=ff(sig,dx,dy,M_H,M_V,Mc,tau,V_H,V_V,theta);
            fprime=ffprime(sig,dx,dy,M_H,M_V,Mc,tau,V_H,V_V,theta);
        }
        
        // f=sig*tau/Mc+theta*V_H+(1-theta)*V_V;
        
        theta=theta-f/fprime;
        
        int stat;
        if(isnan(theta)){
            std::cout << sig << "\t" << M_H << "\t" << M_V << "\t" << V_H << "\t" << V_V <<  std::endl;
            mexErrMsgIdAndTxt("MonProg:newtonEikonal", "problème de convergence du Newton");
        }
        
        err=fabs(theta_o-theta);
    }
    theta=fmin(theta,1);
    theta=fmax(theta,0);
        
    tau=sqrt(pow(dx*theta,2)+pow(dy*(1-theta),2));
    Mc=theta*M_H+(1-theta)*M_V;

    // std::cout << k << "\t" << theta << std::endl;
    *V = sig*tau/Mc + theta*V_H + (1-theta)*V_V;
    *M = Mc;

}

void updateHeap(custom_priority_queue<node> *narrow, std::vector<intPair > *update, double **d, double **mob, int **status, node **Tnarrow, double dx, double dy, double Nx, double Ny){
    
    for(int k=0;k<(*update).size();k++){
        
        int i=(*update).at(k).first;
        int j=(*update).at(k).second;
        int im=i-1;
        int ip=i+1;
        int jm=j-1;
        int jp=j+1;
        double V,V1,V2,V3,V4,V_H,V_V,Pmax,M,M1,M2,M3,M4,M_H,M_V;
        
        double Vc[4]={inf, inf, inf, inf};
        double Mc[4]={inf, inf, inf, inf};
        
        int Hneig[2]={0,0};
        int Vneig[2]={0,0};
        
        double sig = copysign(1.0, d[i][j]);
        
        if(status[im][j]==3 && im!=0){
            V1=d[im][j];
            M1=mob[im][j];
            Hneig[0]=1;
        }
        else
            V1=inf;
        
        if(status[ip][j]==3 && ip!=Nx+1){
            V2=d[ip][j];
            M2=mob[ip][j];
            Hneig[1]=1;
        }
        else
            V2=inf;
        
        if(status[i][jm]==3 && jm!=0){
            V3=d[i][jm];
            M3=mob[i][jm];
            Vneig[0]=1;
        }
        else
            V3=inf;
        
        if(status[i][jp]==3 && jp!=Ny+1){
            V4=d[i][jp];
            M4=mob[i][jp];
            Vneig[1]=1;
        }
        else
            V4=inf;
        
        if(Vneig[0]+Vneig[1]==0){
            if(Hneig[0]==1){
                Vc[0]=V1+sig*dx/M1;
                Mc[0]=M1;
            }
            
            if(Hneig[1]==1){
                Vc[1]=V2+sig*dx/M2;
                Mc[1]=M2;
            }
        }
        else if(Hneig[0]+Hneig[1]==0){
            if(Vneig[0]==1){
                Vc[0]=V3+sig*dy/M3;
                Mc[0]=M3;
            }
            if(Vneig[1]==1){
                Vc[1]=V4+sig*dy/M4;
                Mc[1]=M4;
            }
        }
        else{
            if(Hneig[0]*Vneig[0]==1){
                V_H=V1; M_H=M1;
                V_V=V3; M_V=M3;
                newtonEikonal(&Vc[0],&Mc[0],V_H,V_V,M_H,M_V,dx,dy,sig);
            }
            if(Hneig[1]*Vneig[0]==1){
                V_H=V2; M_H=M2;
                V_V=V3; M_V=M3;
                newtonEikonal(&Vc[1],&Mc[1],V_H,V_V,M_H,M_V,dx,dy,sig);
            }
            if(Hneig[1]*Vneig[1]==1){
                V_H=V2; M_H=M2;
                V_V=V4; M_V=M4;
                newtonEikonal(&Vc[2],&Mc[2],V_H,V_V,M_H,M_V,dx,dy,sig);
            }
            if(Hneig[0]*Vneig[1]==1){
                V_H=V1; M_H=M1;
                V_V=V4; M_V=M4;
                newtonEikonal(&Vc[3],&Mc[3],V_H,V_V,M_H,M_V,dx,dy,sig);
            }
        }
        
        V=Vc[0];
        M=Mc[0];
        
        for(int i=1 ; i<4 ; i++){
            if(abs(Vc[i])<abs(V)){
                V=Vc[i];
                M=Mc[i];
            }
            
        }
        
        mob[i][j]=M;
        
        node El;
        El.i=i;
        El.j=j;
        El.val=fabs(V);
        (*narrow).push(El);
        Tnarrow[i][j]=El;
    }
    
    (*update).erase((*update).begin(),(*update).end());
    
}


void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    
    // Input (d,dx,dy)
    
    double multiplier;              /* input scalar */
    double *inMatrix;               /* 1xN input matrix */
    double *inMatrix2;               /* 1xN input matrix */
    size_t ncols;                   /* size of matrix */
    double *outMatrix;              /* output matrix */
    double *outMatrix2;              /* output matrix */
    
    /* create a pointer to the real data in the input matrix  */
    inMatrix = mxGetPr(prhs[0]);
    
    /* get dimensions of the input matrix */
    
    int Nx = mxGetM(prhs[0]);
    int Ny = mxGetN(prhs[0]);
    
    double **d = (double **)malloc((Nx+2)*sizeof(double *));
    for (int i=0; i<Nx+2; i++)
        d[i] = (double *)malloc((Ny+2) * sizeof(double));
    
    for(int i=1;i<Nx+1;i++)
        for(int j=1;j<Ny+1;j++)
            d[i][j]=inMatrix[(j-1)+Ny*(i-1)];
    
    for(int i=0;i<Nx+2;i++){
        d[i][0]=inf;
        d[i][Ny+1]=inf;
    }
    
    for(int j=0;j<Ny+2;j++){
        d[0][j]=inf;
        d[Nx+1][j]=inf;
    }
    
    // Mobility matrix
    
    inMatrix2 = mxGetPr(prhs[1]);
    
    double **mob = (double **)malloc((Nx+2)*sizeof(double *));
    for (int i=0; i<Nx+2; i++)
        mob[i] = (double *)malloc((Ny+2) * sizeof(double));
    
    for(int i=1;i<Nx+1;i++)
        for(int j=1;j<Ny+1;j++)
            mob[i][j]=inMatrix2[(j-1)+Ny*(i-1)];
    
    for(int i=0;i<Nx+2;i++){
        mob[i][0]=1;
        mob[i][Ny+1]=1;
    }
    
    for(int j=0;j<Ny+2;j++){
        mob[0][j]=1;
        mob[Nx+1][j]=1;
    }
    
    
    
//     for(int i=1;i<Nx+1;i++){
//         for(int j=1;j<Ny+1;j++){
//             std::cout << d[i][j] << "\t" ;
//             }
//         std::cout << "\n";
//     }
    
    double dx,dy;
    
    dx=*mxGetPr(prhs[2]);
    dy=*mxGetPr(prhs[3]);
    
    
    int **status = (int **)malloc((Nx+2)*sizeof(int *));
    for (int i=0; i<Nx+2; i++)
        status[i] = (int *)malloc((Ny+2) * sizeof(int));
    
    node **Tnarrow = (node **)malloc((Nx+2)*sizeof(node *));
    for (int i=0; i<Nx+2; i++)
        Tnarrow[i] = (node *)malloc((Ny+2) * sizeof(node));
    // int k=1;
    for(int i=2;i<Nx;i++){
        for(int j=2;j<Ny;j++){
            int im=i-1;
            int ip=i+1;
            int jm=j-1;
            int jp=j+1;
            if(d[i][j]*d[ip][j]<0 || d[i][j]*d[im][j]<0 || d[i][j]*d[i][jp]<0 || d[i][j]*d[i][jm]<0){
                status[i][j]=3;
                //    std::cout << k << std::endl;
                //    k++;
            }
            else
                status[i][j]=1;
        }
    }
    
    for(int i=0;i<Nx+2;i++){
        status[i][0]=3;
        status[i][Ny+1]=3;
    }
    
    for(int j=0;j<Ny+2;j++){
        status[0][j]=3;
        status[Nx+1][j]=3;
    }
    
    for(int i=1;i<Nx+1;i++){
        status[i][1]=1;
        status[i][Ny]=1;
    }
    
    for(int j=1;j<Ny+1;j++){
        status[1][j]=1;
        status[Nx][j]=1;
    }
    
    custom_priority_queue<node> narrow;
    std::vector<intPair > update;
    
    for(int i=2;i<Nx;i++){
        for(int j=2;j<Ny;j++){
            if(status[i][j]==1){
                int im=i-1;
                int ip=i+1;
                int jm=j-1;
                int jp=j+1;
                if(status[ip][j]==3 || status[im][j]==3 ||status[i][jp]==3 || status[i][jm]==3){
                    update.push_back(intPair(i,j));
                    status[i][j]=2;
                }
            }
        }
    }
    
    // std::cout << "totot " << update.size() << "\n";
    
    updateHeap(&narrow,&update,d,mob,status,Tnarrow,dx,dy,Nx,Ny);
    
    
    while(!narrow.empty()){
        
        node Elmin =narrow.top();
        narrow.pop();
        
        int i=Elmin.i;
        int j=Elmin.j;
        
        //std::cout << "totot " << narrow.size() << "\t" << i << "\t" << j << "\n";
        
        
        d[i][j]=copysign(1.0,d[i][j])*Elmin.val;
        status[i][j]=3;
        
        int im=i-1;
        int ip=i+1;
        int jm=j-1;
        int jp=j+1;
        
        // std::cout << "size of status before : " << update.size() << std::endl;
        
        if(status[im][j]<3){
            if(status[im][j]==2)
                narrow.remove(Tnarrow[im][j]);
            update.push_back(intPair(im,j));
            status[im][j]=2;
            //    std::cout << "update im" << std::endl;
        }
        
        if(status[ip][j]<3){
            if(status[ip][j]==2)
                narrow.remove(Tnarrow[ip][j]);
            update.push_back(intPair(ip,j));
            status[ip][j]=2;
            //    std::cout << "update ip" << std::endl;
        }
        
        if(status[i][jm]<3){
            if(status[i][jm]==2)
                narrow.remove(Tnarrow[i][jm]);
            update.push_back(intPair(i,jm));
            status[i][jm]=2;
            //    std::cout << "update jm" << std::endl;
        }
        
        if(status[i][jp]<3){
            if(status[i][jp]==2)
                narrow.remove(Tnarrow[i][jp]);
            update.push_back(intPair(i,jp));
            status[i][jp]=2;
            //   std::cout << "update jp" << std::endl;
        }
        // std::cout << "size of status after : " << update.size() << std::endl;
        
        updateHeap(&narrow,&update,d,mob,status,Tnarrow,dx,dy,Nx,Ny);
        
        
    }
    
    
    /* create the output matrix */
    plhs[0] = mxCreateDoubleMatrix((mwSize)Nx,(mwSize)Ny,mxREAL);
    
    /* get a pointer to the real data in the output matrix */
    outMatrix = mxGetPr(plhs[0]);
    
    for(int i=1;i<Nx+1;i++)
        for(int j=1;j<Ny+1;j++)
            outMatrix[(j-1)+Ny*(i-1)]=d[i][j];
    
    /* create the output matrix */
    plhs[1] = mxCreateDoubleMatrix((mwSize)Nx,(mwSize)Ny,mxREAL);
    
    /* get a pointer to the real data in the output matrix */
    outMatrix2 = mxGetPr(plhs[1]);
    
    for(int i=1;i<Nx+1;i++)
        for(int j=1;j<Ny+1;j++)
            outMatrix2[(j-1)+Ny*(i-1)]=mob[i][j];
    
    
    return;
}
