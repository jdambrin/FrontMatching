#include "mex.h"
#include <math.h>
#include <functional>
#include <queue>
#include <vector>
#include <iostream>
#include <limits>
#include <assert.h>

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
        


void updateHeap(custom_priority_queue<node> *narrow, std::vector<intPair > *update, double **d, double **mob, double **alpha, int **status, node **Tnarrow, double dx, double dy){
    
    for(int k=0;k<(*update).size();k++){
        
        int i=(*update).at(k).first;
        int j=(*update).at(k).second;
        int im=i-1;
        int ip=i+1;
        int jm=j-1;
        int jp=j+1;
        double V,V1,V2,V3,V4,V_H,V_V,Pmax,A,A_H,A_V,A1,A2,A3,A4;
        
        double sig = copysign(1.0, d[i][j]);
        
        double M=mob[i][j];
        
        if(status[im][j]==3){
            V1=d[im][j];
            A1=alpha[im][j];
        }
        else
            V1=inf;
        
        if(status[ip][j]==3){
            V2=d[ip][j];
            A2=alpha[ip][j];
        }
        else
            V2=inf;
        
        if(status[i][jm]==3){
            V3=d[i][jm];
            A3=alpha[i][jm];
        }
        else
            V3=inf;
        
        if(status[i][jp]==3){
            V4=d[i][jp];
            A4=alpha[i][jp];
        }
        else
            V4=inf;
        
        V_H=V1;
        A_H=A1;
        if(V1>V2){
            V_H=V2;
            A_H=A2;
        }
                
        V_V=V3;
        A_V=A3;
        if(V3>V4){
            V_V=V4;
            A_V=A4;
        }
        
        if(V_H==inf){
            V=V_V+sig*dy/M;
            A=A_V;
        }
        else if(V_V==inf){
            V=V_H+sig*dx/M;
            A=A_H;
        }
        else{
            if(V_H<V_V)
                Pmax=pow(V_V-V_H,2)/pow(dx,2);
            else
                Pmax=pow(V_V-V_H,2)/pow(dy,2);
            if(Pmax>(1.0/pow(M,2))){
                if(V_H<V_V){
                    V=V_H+sig*dx/M;
                    A=A_H;
                }
                else{
                    V=V_V+sig*dy/M;
                    A=A_V;
                }
            }
            else{
                double alpha=dx/dy;
                double a=1+pow(alpha,2);
                double b=-2.0*(V_H+pow(alpha,2)*V_V);
                double c=pow(V_H,2)+pow(alpha,2)*pow(V_V,2)-pow(dx,2)/pow(M,2);
                V=(-b+sig*sqrt(fabs(pow(b,2)-4*a*c)))/(2*a); 
                
               // std::cout << a << "\t " << b << "\t" << c << std::endl;
                
                double dxV=(V-V_H)/(dx*dx);
                double dyV=(V-V_V)/(dy*dy);
                
                A=(A_H*dxV+A_V*dyV)/(dxV+dyV);
            }
        }
        
        alpha[i][j]=A;
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
    double *inMatrix3;               /* 1xN input matrix */
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
    
    
    
    inMatrix3 = mxGetPr(prhs[2]);
    
    double **alpha = (double **)malloc((Nx+2)*sizeof(double *));
    for (int i=0; i<Nx+2; i++)
        alpha[i] = (double *)malloc((Ny+2) * sizeof(double));
    
    for(int i=1;i<Nx+1;i++)
        for(int j=1;j<Ny+1;j++)
            alpha[i][j]=inMatrix3[(j-1)+Ny*(i-1)];
    
    for(int i=0;i<Nx+2;i++){
        alpha[i][0]=1;
        alpha[i][Ny+1]=1;
    }
    
    for(int j=0;j<Ny+2;j++){
        alpha[0][j]=1;
        alpha[Nx+1][j]=1;
    }
    
    
//     for(int i=1;i<Nx+1;i++){
//         for(int j=1;j<Ny+1;j++){
//             std::cout << d[i][j] << "\t" ;
//             }
//         std::cout << "\n";
//     }
    
    double dx,dy;
    
    dx=*mxGetPr(prhs[3]);
    dy=*mxGetPr(prhs[4]);
    
    
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
    
    updateHeap(&narrow,&update,d,mob,alpha,status,Tnarrow,dx,dy);
    
    
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
        
        updateHeap(&narrow,&update,d,mob,alpha,status,Tnarrow,dx,dy);
        
        
    }
    
    
    /* create the output matrix */
    plhs[0] = mxCreateDoubleMatrix((mwSize)Nx,(mwSize)Ny,mxREAL);
    plhs[1] = mxCreateDoubleMatrix((mwSize)Nx,(mwSize)Ny,mxREAL);
    
    /* get a pointer to the real data in the output matrix */
    outMatrix = mxGetPr(plhs[0]);
    
    for(int i=1;i<Nx+1;i++)
        for(int j=1;j<Ny+1;j++)
            outMatrix[(j-1)+Ny*(i-1)]=alpha[i][j];
    
     
    outMatrix2 = mxGetPr(plhs[1]);
    
    for(int i=1;i<Nx+1;i++)
        for(int j=1;j<Ny+1;j++)
            outMatrix2[(j-1)+Ny*(i-1)]=d[i][j];
    
    
    return;
}
