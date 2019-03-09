function [Mnew,dnew,Anew]=initBand(M,d,A,dx,dy)
%initialisation de la bande
[Nx,Ny]=size(M);

dnew=zeros(Nx,Ny);
Mnew=zeros(Nx,Ny);
Anew=zeros(Nx,Ny);

for i=1:Nx
    for j=1:Ny
        
        ip=mod(i+1-1,Nx)+1;
        im=mod(i-1-1,Nx)+1;
        jp=mod(j+1-1,Ny)+1;
        jm=mod(j-1-1,Ny)+1;
        type=[0 0 0 0];
        dNc=[];
        MNc=[];
        ANc=[];
        
        if(d(ip,j)*d(i,j)<0)
            theta=d(i,j)/(d(i,j)-d(ip,j));
            Mcand=theta*M(i,j)+(1-theta)*M(ip,j);
            Acand=theta*A(i,j)+(1-theta)*A(ip,j);
            dcand=dx*theta/Mcand;
            dNc=[dNc dcand]; MNc=[MNc Mcand]; ANc=[ANc Acand];
        end
        
        if(d(im,j)*d(i,j)<0)
            theta=d(i,j)/(d(i,j)-d(im,j));
            Mcand=theta*M(i,j)+(1-theta)*M(im,j);
            Acand=theta*A(i,j)+(1-theta)*A(im,j);
            dcand=dx*theta/Mcand;
            dNc=[dNc dcand]; MNc=[MNc Mcand]; ANc=[ANc Acand];
        end
        
        if(d(i,jp)*d(i,j)<0)
            theta=d(i,j)/(d(i,j)-d(i,jp));
            Mcand=theta*M(i,j)+(1-theta)*M(i,jp);
            Acand=theta*A(i,j)+(1-theta)*A(i,jp);
            dcand=dy*theta/Mcand;
            dNc=[dNc dcand]; MNc=[MNc Mcand]; ANc=[ANc Acand];
        end
        
        if(d(i,jm)*d(i,j)<0)
            theta=d(i,j)/(d(i,j)-d(i,jm));
            Mcand=theta*M(i,j)+(1-theta)*M(i,jm);
            Acand=theta*A(i,j)+(1-theta)*A(i,jm);
            dcand=dy*theta/Mcand;
            dNc=[dNc dcand]; MNc=[MNc Mcand]; ANc=[ANc Acand];
        end
        
        if(size(MNc)>0)
            [dnew(i,j),ind]=min(dNc);
            dnew(i,j)=sign(d(i,j))*dnew(i,j);
            Mnew(i,j)=MNc(ind);
            Anew(i,j)=ANc(ind);
        else
            Mnew(i,j)=0;
            Anew(i,j)=0;
            dnew(i,j)=sign(d(i,j));
        end
               
    end
end
end