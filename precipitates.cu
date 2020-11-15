// LSW probability density function
REAL PDF(REAL X)
{
    return (powf(3.,4)*expf(1.)/powf(2.,5./3.))*X*X*expf(-1./(1-2./3.*X))/(powf(X+3.,7./3.)*powf(1.5-X,11./3.));
}

REAL CDF(REAL X, REAL dX, REAL U1)
{
    REAL Xi = 0.;
    REAL Xtemp;
    REAL F = 0.;
    // trapeze integration method
    while(Xi<=X)
    {
        Xtemp = Xi+dX;
        F = F + 0.5*dX*(PDF(Xi)+PDF(Xtemp));
        Xi = Xtemp;
    }
    return F-U1; // g(x)-y
}

void InitializePrecipitates(Precipitate *ppts, Parameters *par)
{
    int Nppts = par[0].Nppts;
    int posi, posj, posk;
    REAL Rtemp;
    int flag;
    int pptsShape = par[0].Fppts;
    REAL nu = par[0].nu;

    // Single precipitate
    if(strcmp(par->dist_name,"Single")==0)
    {
        ppts[0].state = 1;
        ppts[0].form = pptsShape;
        ppts[0].nu = nu;
        ppts[0].i = par[0].Nx/2;
        ppts[0].j = par[0].Ny/2;
        ppts[0].k = par[0].Nz/2;
        ppts[0].R = par[0].Rinit_1;
        
        //sphere integration box
        if (ppts[0].form == 0){
            ppts[0].hx = int(ppts[0].R +2); //CASTER EN INT PARTOUT rajouter en +2 en input 
            ppts[0].hy = int(ppts[0].R +2);
            ppts[0].hz = int(ppts[0].R +2);
        }
        
        //cylinder integration box
        if (ppts[0].form == 1){
            ppts[0].hx = int(ppts[0].R +2);
            ppts[0].hy = int(ppts[0].R +2);
            ppts[0].hz = int(ppts[0].R*ppts[0].nu +2);
        }
        
        
        //elipsoid integration box
        if (ppts[0].form == 2){
            ppts[0].hx = int(ppts[0].R +2);
            ppts[0].hy = int(ppts[0].R*ppts[0].nu +2);
            ppts[0].hz = int(ppts[0].R*ppts[0].nu +2);
        }
        
        
        for(int n=1; n<Nppts; n++)
        {
            ppts[n].state = 0;
            ppts[n].i =0;
            ppts[n].j =0;
            ppts[n].k =0;
            ppts[n].R =0;
        }   
    }
    
    
    //2 ppts 
    
    if(strcmp(par->dist_name,"TwoPpts")==0)
    {
        int R = int(round(par[0].Rinit_2/2.));
        ppts[0].state = 1;
        ppts[0].form = pptsShape;
        ppts[0].nu = nu;
        ppts[0].i = par[0].Nx/2-R;
        ppts[0].j = par[0].Ny/2;
        ppts[0].k = par[0].Nz/2;
        ppts[0].R = par[0].Rinit_1;
        
        ppts[1].state = 1;
        ppts[1].form = pptsShape;
        ppts[1].nu = nu;
        ppts[1].i = par[0].Nx/2+R;
        ppts[1].j = par[0].Ny/2;
        ppts[1].k = par[0].Nz/2;
        ppts[1].R = par[0].Rinit_1-1e-3;
        
        for(int n=0; n<Nppts; n++){
                //sphere integration box
            if (ppts[n].form == 0){
                ppts[n].hx = int(ppts[n].R +2);
                ppts[n].hy = int(ppts[n].R +2);
                ppts[n].hz = int(ppts[n].R +2);
            }
            
            //cylinder integration box
            if (ppts[n].form == 1){
                ppts[n].hx = int(ppts[n].R +2);
                ppts[n].hy = int(ppts[n].R +2);
                ppts[n].hz = int(ppts[n].R*ppts[n].nu +2);
            }
            
            
            //elipsoid integration box
            if (ppts[n].form == 2){
                ppts[n].hx = int(ppts[n].R +2);
                ppts[n].hy = int(ppts[n].R*ppts[n].nu +2);
                ppts[n].hz = int(ppts[n].R*ppts[n].nu +2);
            }
        }
        
        for(int n=2; n<Nppts; n++)
        {
            ppts[n].state = 0;
        }   
    }
    
    // normal distribution with the Box-Muller method
    if(strcmp(par->dist_name,"Normal")==0)
    {
        REAL U1, U2;
        REAL Rmin2 = 16*par[0].Rinit_1*par[0].Rinit_1;
    
        for(int n=0; n<Nppts; n++)
        {
            posi = int(REAL(rand())/REAL(RAND_MAX)*par[0].Nx); 
            posj = int(REAL(rand())/REAL(RAND_MAX)*par[0].Ny);
            posk = int(REAL(rand())/REAL(RAND_MAX)*par[0].Nz);
            
            U1 = REAL(rand())/REAL(RAND_MAX);
            U2 = REAL(rand())/REAL(RAND_MAX);
            Rtemp = par[0].Rinit_1 + (sqrt(-2*log(U1))*cos(2*M_PI*U2))*par[0].Rinit_2;
            
            flag = 1;
            
            // if the precipitate radius is too small, remove it
            if(Rtemp<par[0].Rth)
            {
                flag = 0;
            }

            // if the precipitate is too close to other precipitates, remove it
            else
            {            
                for(int m=0;m<n;m++)
                {
                    if((REAL((posi-ppts[m].i)*(posi-ppts[m].i)+(posj-ppts[m].j)*(posj-ppts[m].j)+(posk-ppts[m].k)*(posk-ppts[m].k))<Rmin2))
                    {
                        flag = 0;
                        break;
                    }
                }
            }
            
            if(flag==1)
            {
                ppts[n].state = 1;
                ppts[n].form = pptsShape;
                ppts[n].nu = nu;
                ppts[n].i = posi;
                ppts[n].j = posj;
                ppts[n].k = posk;
                ppts[n].R = Rtemp;
                
                for(int nn=0; nn<Nppts; nn++){
                        //sphere integration box
                    if (ppts[nn].form == 0){
                        ppts[nn].hx = int(ppts[nn].R +2);
                        ppts[nn].hy = int(ppts[nn].R +2);
                        ppts[nn].hz = int(ppts[nn].R +2);
                    }
                    
                    //cylinder integration box
                    if (ppts[nn].form == 1){
                        ppts[nn].hx = int(ppts[nn].R +2);
                        ppts[nn].hy = int(ppts[nn].R +2);
                        ppts[nn].hz = int(ppts[nn].R*ppts[nn].nu +2);
                    }
                    
                    
                    //elipsoid integration box
                    if (ppts[nn].form == 2){
                        ppts[nn].hx = int(ppts[nn].R +2);
                        ppts[nn].hy = int(ppts[nn].R*ppts[nn].nu +2);
                        ppts[nn].hz = int(ppts[nn].R*ppts[nn].nu +2);
                    }
                }
            }
            else
            {
                n=n-1;
            }
        }        
    }
    
    // LSW distribution with dichotomy method
    if(strcmp(par->dist_name,"LSW")==0)
    {
        REAL U1;
        REAL Rmin2 = 16*par[0].Rinit_1*par[0].Rinit_1;
        REAL Dmax, Dbc;
    
        for(int n=0; n<Nppts; n++)
        {
            posi = int(REAL(rand())/REAL(RAND_MAX)*par[0].Nx); 
            posj = int(REAL(rand())/REAL(RAND_MAX)*par[0].Ny);
            posk = int(REAL(rand())/REAL(RAND_MAX)*par[0].Nz);
            REAL X0=0.001;
            REAL X1=1.499;
            REAL Xf,F0,Ff;
            REAL dX=0.0001;
            U1 = REAL(rand())/REAL(RAND_MAX); // y
            // dichotomy to find Xf for g(Xf) = y
            while(X1 - X0>=dX)
            {    
                Xf = 0.5*(X0+X1);
                F0 = CDF(X0,dX,U1); // g(X0) - y
                Ff = CDF(Xf,dX,U1); // g(Xf) - y
                
                if(F0*Ff<=0.)
                {    
                    X1 = Xf;
                }
                else
                {
                    X0 = Xf;
                }
            }
            
            Rtemp = par[0].Rinit_1*Xf; // Ri = R0 * Xf
            flag = 1;
            
            if(Rtemp<par[0].Rth)
            {
                flag = 0;
            }
            else
            {            
                for(int m=0;m<n;m++)
                {
                    // Checking if precipitates are not overlapping
                    if((REAL((posi-ppts[m].i)*(posi-ppts[m].i)+(posj-ppts[m].j)*(posj-ppts[m].j)+(posk-ppts[m].k)*(posk-ppts[m].k))<Rmin2))
                    {
                        Dmax = sqrt(3.)/2.*(ppts[m].R+Rtemp)+1.;
                        Dbc = sqrt(REAL((posi-ppts[m].i)*(posi-ppts[m].i)+(posj-ppts[m].j)*(posj-ppts[m].j)+(posk-ppts[m].k)*(posk-ppts[m].k)));
                        if(Dbc<=Dmax)
                        {
                            flag = 0;
                            break;
                        }
                    }
                }
            }
            
            if(flag==1)
            {
                ppts[n].state = 1;
                ppts[n].form = pptsShape;
                ppts[n].nu = nu;
                ppts[n].i = posi;
                ppts[n].j = posj;
                ppts[n].k = posk;
                ppts[n].R = Rtemp;
                printf("Precipitate number %d 's radius is %f nm\n", n+1, ppts[n].R*par[0].lc);
                for(int nn=0; nn<Nppts; nn++){
                        //sphere integration box
                    if (ppts[nn].form == 0){
                        ppts[nn].hx = int(ppts[nn].R +2);
                        ppts[nn].hy = int(ppts[nn].R +2);
                        ppts[nn].hz = int(ppts[nn].R +2);
                    }
                    
                    //cylinder integration box
                    if (ppts[nn].form == 1){
                        ppts[nn].hx = int(ppts[nn].R +2);
                        ppts[nn].hy = int(ppts[nn].R +2);
                        ppts[nn].hz = int(ppts[nn].R*ppts[nn].nu +2);
                    }
                    
                    
                    //elipsoid integration box
                    if (ppts[nn].form == 2){
                        ppts[nn].hx = int(ppts[nn].R +2);
                        ppts[nn].hy = int(ppts[nn].R*ppts[nn].nu +2);
                        ppts[nn].hz = int(ppts[nn].R*ppts[nn].nu +2);
                    }
                }
            }
                
            else
            {
                n=n-1;
            }
        }
    }
}


__global__ void PrecipitatesGrowth(REAL *U, REAL *Phi, Precipitate *ppts, Parameters *par) 
{

    int n = blockIdx.x;

    if(ppts[n].state==1)
    {
            int Nx=par[0].Nx;
            int Ny=par[0].Ny;
            int Nz=par[0].Nz;

            int ip=ppts[n].i;
            int jp=ppts[n].j;
            int kp=ppts[n].k;
            
            int i,j,k;
            int ii, jj, kk;
            
            REAL Int=0;
            int Ri=int(ppts[n].R);
            //integration box coordinate whatever ppts.form
            int hx = int(ppts[n].hx);
            int hy = int(ppts[n].hy);
            int hz = int(ppts[n].hz);
            
            int h=Ri+2;
            
            
            int l, m;
            int iiN, jjN, kkN;
                
            //-------//
            // South //
            //-------//      
            j = jp-hy+1;
            for(l=-hx+1;l<hx;l++)
            {
                for(m=-hz+1;m<hz;m++)
                {
                    ii = mod((ip+l),Nx);
                    jj = mod(j,Ny);
                    kk = mod((kp+m),Nz);
                    jjN = mod(j-1,Ny);
                    Int = Int + (U[ii*Ny*Nz+jjN*Nz+kk] - U[ii*Ny*Nz+jj*Nz+kk]);
                }
            }

            //------//
            // East //
            //------//       
            i = ip+hx-1;
            for(l=-hy+1;l<hy;l++)
            {
                for(m=-hz+1;m<hz;m++)
                {
                   ii = mod(i,Nx);
                   jj = mod((jp+l),Ny);
                   kk = mod((kp+m),Nz);
                   iiN = mod(i+1,Nx);
                   Int = Int + (U[iiN*Ny*Nz+jj*Nz+kk] - U[ii*Ny*Nz+jj*Nz+kk]);
                }
            }
                
            //-------//
            // North //
            //-------//
            j = jp+hy-1;
            for(l=-hx+1;l<hx;l++)
            {
                for(m=-hz+1;m<hz;m++)
                {
                    ii = mod((ip+l),Nx);
                    jj = mod(j,Ny);
                    kk = mod((kp+m),Nz);
                    jjN = mod(j+1,Ny);
                    Int = Int + (U[ii*Ny*Nz+jjN*Nz+kk] - U[ii*Ny*Nz+jj*Nz+kk]);
                }
            }
            
            //------//
            // West //
            //------//
            i = ip-hx+1;
            for(l=-hy+1;l<hy;l++)
            {
                for(m=-hz+1;m<hz;m++)
                {
                    ii = mod(i,Nx);
                    jj = mod((jp+l),Ny);
                    kk = mod((kp+m),Nz);
                    iiN = mod(i-1,Nx);
                    Int = Int + (U[iiN*Ny*Nz+jj*Nz+kk] - U[ii*Ny*Nz+jj*Nz+kk]);
                }
            }
              
            //-------//
            // Above //
            //-------//
            k = kp+hz-1;
            for(l=-hx+1;l<hx;l++)
            {
                for(m=-hz+1;m<hz;m++)
                {
                    ii = mod((ip+l),Nx);
                    jj = mod((jp+m),Ny);
                    kk = mod(k,Nz);
                    kkN = mod(k+1,Nz);
                    Int = Int + (U[ii*Ny*Nz+jj*Nz+kkN] - U[ii*Ny*Nz+jj*Nz+kk]);
                }
            }
                
            //-------//
            // Under //
            //-------//
            k = kp-hz+1;
            for(l=-hx+1;l<hx;l++)
            {
                for(m=-hy+1;m<hy;m++)
                {
                    ii = mod((ip+l),Nx);
                    jj = mod((jp+m),Ny);
                    kk = mod(k,Nz);
                    kkN = mod(k-1,Nz);
                    Int = Int + (U[ii*Ny*Nz+jj*Nz+kkN] - U[ii*Ny*Nz+jj*Nz+kk]);
                }
            }
            
            REAL ceq_R = par[0].ceq_m*exp(par[0].d0/ppts[n].R);
            //GT condition with respect to ppts's form need to be changed 
            if (ppts[n].form == 0){
                REAL ceq_R = par[0].ceq_m*exp(par[0].d0/ppts[n].R); // Concentration at the interface matrix/precipitate
            }
            if (ppts[n].form == 1){
                REAL ceq_R = par[0].ceq_m*exp(par[0].d0/ppts[n].R); // Concentration at the interface matrix/precipitate
                //exp(4*(1+ppts[n].nu)*v_at*g/(3*ppts[n].nu*ppts[n].R*4.11*pow(10,-21) *par[0].ceq_p))
            }
            if (ppts[n].form == 2){
                REAL ceq_R = par[0].ceq_m*exp(par[0].d0/ppts[n].R); // Concentration at the interface matrix/precipitate
            }
            REAL v = (par[0].D*(par[0].ceq_p - par[0].ceq_m))/(4.*M_PI*(par[0].ceq_p - ceq_R)*ppts[n].R*ppts[n].R)*Int; //Stefan condition with         GT          generaiser a differente form 
            ppts[n].R += par[0].dt*v;

    }
}


__global__ void PrecipitatesOverlap(REAL *U, REAL *Phi, Precipitate *ppts, Parameters *par) //looking if there is overlapping, if so remove the smallest
{
	//Calculate distance between other ppts 
	int n = blockIdx.x;
	int di ;
	int dj ;
	int dk ;
	int ip ; 
	int jp; 
	int kp;
	REAL dist ;
	int iter ;
	REAL Rnew; 
	if (ppts[n].state ==1) {
		
		for (iter = 0; iter<par[0].Nppts; iter++){ //looking for all other ppts 
		
			if ((iter != n) && (ppts[iter].state ==1)){
			
				di = ppts[n].i - ppts[iter].i;
				dj = ppts[n].j - ppts[iter].j;
				dk = ppts[n].k - ppts[iter].k; 
				dist = sqrt(REAL(di*di+dj*dj+dk*dk)); 
				
				if (dist < (ppts[iter].R + ppts[n].R)*sqrtf(2)){ // if it's too close remose the smallest (integration box ==>Rmax=Rsqrt(2))
					
					if (ppts[iter].R < ppts[n].R){
						//calculation of the new coordinates of the PPT
						ip = round((ppts[iter].R*ppts[iter].i+ppts[n].R*ppts[n].i)/(ppts[iter].R+ppts[n].R));
						jp = round((ppts[iter].R*ppts[iter].j+ppts[n].R*ppts[n].j)/(ppts[iter].R+ppts[n].R));
						kp = round((ppts[iter].R*ppts[iter].k+ppts[n].R*ppts[n].k)/(ppts[iter].R+ppts[n].R));
						//change coordiates 
						ppts[n].i = ip;
						ppts[n].j = jp;
						ppts[n].k = kp;						
						//new radius of the PPT 
						Rnew = pow(ppts[n].R*ppts[n].R*ppts[n].R + ppts[iter].R*ppts[iter].R*ppts[iter].R,1./3.); 
						ppts[n].R = Rnew;
						ppts[iter].state = 0; 
					}
					
				}
				
			}
			
		}
	}
}

__global__ void UpdateFields(REAL *U, REAL *Phi, Precipitate *ppts, Parameters *par) 
{

    int n = blockIdx.x;

    if(ppts[n].state==1)
    {
        int Nx=par[0].Nx;
        int Ny=par[0].Ny;
        int Nz=par[0].Nz;

        int ip=ppts[n].i;
        int jp=ppts[n].j;
        int kp=ppts[n].k;
        
        int i,j,k;
        int ii,jj,kk;
        
        //integration box update
        //Sphere integration box
        if (ppts[n].form == 0){
            ppts[n].hx = int(ppts[n].R +2);
            ppts[n].hy = int(ppts[n].R +2);
            ppts[n].hz = int(ppts[n].R +2);
                    
        }
                    
        //cylinder integration box
        if (ppts[n].form == 1){
            ppts[n].hx = int(ppts[n].R +2);
            ppts[n].hy = int(ppts[n].R +2);
            ppts[n].hz = int(ppts[n].R*ppts[n].nu +2);
        }
                    
                    
        //elipsoid integration box
        if (ppts[n].form == 2){
            ppts[n].hx = int(ppts[n].R +2);
            ppts[n].hy = int(ppts[n].R*ppts[n].nu +2);
            ppts[n].hz = int(ppts[n].R*ppts[n].nu +2);
        }
        
        int h=int(ppts[n].R)+1;
        int hx = int(ppts[n].hx)-1;
        int hy = int(ppts[n].hy)-1;
        int hz = int(ppts[n].hz)-1;
        
        REAL B2 = ppts[n].R*ppts[n].nu*ppts[n].R*ppts[n].nu; 
           
        
        REAL Rppt2 = ppts[n].R*ppts[n].R;
            
        // if the precipitate radius is smaller than a threshold, it disappears
        if(ppts[n].R<par[0].Rth)
        {
            ppts[n].state=0; //switch of precipitate
            int N=0; 
           
            // count number of grid points in ppt
            if (ppts[n].form == 0){ //sphere
                
                for(i=ip-hx;i<ip+hx+1;i++)
                {
                    for(j=jp-hy;j<jp+hy+1;j++)
                    {
                        for(k=kp-hz;k<kp+hz+1;k++)
                        {
                            if((i-ip)*(i-ip)+(j-jp)*(j-jp)+(k-kp)*(k-kp) <= Rppt2)
                            {
                                N = N + 1;
                            }
                        }
                    }
                }
            }
            if (ppts[n].form == 1){ //cylinder
                
                for(i=ip-hx;i<ip+hx+1;i++)
                {
                    for(j=jp-hy;j<jp+hy+1;j++)
                    {
                        for(k=kp-hz;k<kp+hz+1;k++)
                        {
                            if(((i-ip)*(i-ip)+(j-jp)*(j-jp)<= Rppt2) && abs(k-kp)<ppts[n].nu*ppts[n].R) //cynlindrical condition 
                            {
                                N = N + 1;

                            }
                        }
                    }
                }
            }
            if (ppts[n].form == 2){ //elipsoïd
                for(i=ip-hx;i<ip+hx+1;i++)
                {
                    for(j=jp-hy;j<jp+hy+1;j++)
                    {
                        for(k=kp-hz;k<kp+hz+1;k++)
                        {
                            if(REAL((i-ip)*(i-ip))/Rppt2+REAL((j-jp)*(j-jp))/B2+REAL((k-kp)*(k-kp))/B2 <= 1.0) //elipsoïd condition  
                            {
                                N = N + 1;
                            }
                        }
                    }
                }
            }
            
            // redistribute solute
            REAL conc = par[0].ceq_p*4./3.*M_PI*ppts[n].R*ppts[n].R*ppts[n].R/REAL(N); 
            if (ppts[n].form == 0){ //sphere
                
                for(i=ip-hx;i<ip+hx+1;i++)
                {
                    for(j=jp-hy;j<jp+hy+1;j++)
                    {
                        for(k=kp-hz;k<kp+hz+1;k++)
                        {
                            ii = mod(i,Nx);
                            jj = mod(j,Ny);
                            kk = mod(k,Nz);
                            if((i-ip)*(i-ip)+(j-jp)*(j-jp)+(k-kp)*(k-kp) <= Rppt2)
                            {
                                Phi[ii*Ny*Nz+jj*Nz+kk] = 0.0;
                                U[ii*Ny*Nz+jj*Nz+kk] = (par[0].ceq_m - conc)/(par[0].ceq_m - par[0].ceq_p);
                            }
                        }
                    }
                }
            }
            if (ppts[n].form == 1){ //cylinder
                
                for(i=ip-hx;i<ip+hx+1;i++)
                {
                    for(j=jp-hy;j<jp+hy+1;j++)
                    {
                        for(k=kp-hz;k<kp+hz+1;k++)
                        {
                            ii = mod(i,Nx);
                            jj = mod(j,Ny);
                            kk = mod(k,Nz);
                            if(((i-ip)*(i-ip)+(j-jp)*(j-jp)<= Rppt2) && abs(k-kp)<ppts[n].nu*ppts[n].R) //cynlindrical condition 
                            {
                                Phi[ii*Ny*Nz+jj*Nz+kk] = 0.0;
                                U[ii*Ny*Nz+jj*Nz+kk] = (par[0].ceq_m - conc)/(par[0].ceq_m - par[0].ceq_p);

                            }
                        }
                    }
                }
            }
            if (ppts[n].form == 2){ //elipsoïd
                for(i=ip-hx;i<ip+hx+1;i++)
                {
                    for(j=jp-hy;j<jp+hy+1;j++)
                    {
                        for(k=kp-hz;k<kp+hz+1;k++)
                        {
                            ii = mod(i,Nx);
                            jj = mod(j,Ny);
                            kk = mod(k,Nz);
                            if(REAL((i-ip)*(i-ip))/Rppt2+REAL((j-jp)*(j-jp))/B2+REAL((k-kp)*(k-kp))/B2 <= 1.0) //elipsoïd condition 
                            {
                                Phi[ii*Ny*Nz+jj*Nz+kk] = 0.0;
                                U[ii*Ny*Nz+jj*Nz+kk] = (par[0].ceq_m - conc)/(par[0].ceq_m - par[0].ceq_p);
                            }
                        }
                    }
                }
            }
        }
    
        else // set phi=1 and impose concentration on the ppt
        {
            if (ppts[n].form == 0){ //sphere
                // Concentration at matrix/precipitate interface (Gibbs-Thomson)
                REAL U_R = par[0].ceq_m*(exp(par[0].d0/ppts[n].R)-1.)/(par[0].ceq_p - par[0].ceq_m); 
                
                for(i=ip-hx;i<ip+hx+1;i++)
                {
                    for(j=jp-hy;j<jp+hy+1;j++)
                    {
                        for(k=kp-hz;k<kp+hz+1;k++)
                        {
                            ii = mod(i,Nx);
                            jj = mod(j,Ny);
                            kk = mod(k,Nz);
                            Phi[ii*Ny*Nz+jj*Nz+kk] = 0.0;
                            if((i-ip)*(i-ip)+(j-jp)*(j-jp)+(k-kp)*(k-kp) <= Rppt2)
                            {
                                Phi[ii*Ny*Nz+jj*Nz+kk] = 1.0;
                                
                                U[ii*Ny*Nz+jj*Nz+kk] = U_R;
                            }
                        }
                    }
                }
            }
            if (ppts[n].form == 1){ //cylinder
                // Concentration at matrix/precipitate interface (Gibbs-Thomson)
                REAL U_R = par[0].ceq_m*(exp(par[0].d0/ppts[n].R)-1.)/(par[0].ceq_p - par[0].ceq_m); // change it
                
                for(i=ip-hx;i<ip+hx+1;i++)
                {
                    for(j=jp-hy;j<jp+hy+1;j++)
                    {
                        for(k=kp-hz;k<kp+hz+1;k++)
                        {
                            ii = mod(i,Nx);
                            jj = mod(j,Ny);
                            kk = mod(k,Nz);
                            Phi[ii*Ny*Nz+jj*Nz+kk] = 0.0;
                            if(((i-ip)*(i-ip)+(j-jp)*(j-jp)<= Rppt2) && abs(k-kp)<ppts[n].nu*ppts[n].R) //cynlindrical condition 
                            {
                                Phi[ii*Ny*Nz+jj*Nz+kk] = 1.0;
                                
                                U[ii*Ny*Nz+jj*Nz+kk] = U_R;
                            }
                        }
                    }
                }
            }
            if (ppts[n].form == 2){ //elipsoïd
                // Concentration at matrix/precipitate interface (Gibbs-Thomson)
                REAL U_R = par[0].ceq_m*(exp(par[0].d0/ppts[n].R)-1.)/(par[0].ceq_p - par[0].ceq_m); //	 change it
                
                for(i=ip-hx;i<ip+hx+1;i++)
                {
                    for(j=jp-hy;j<jp+hy+1;j++)
                    {
                        for(k=kp-hz;k<kp+hz+1;k++)
                        {
                            ii = mod(i,Nx);
                            jj = mod(j,Ny);
                            kk = mod(k,Nz);
                            Phi[ii*Ny*Nz+jj*Nz+kk] = 0.0;
                            if(REAL((i-ip)*(i-ip))/Rppt2+REAL((j-jp)*(j-jp))/B2+REAL((k-kp)*(k-kp))/B2 <= 1.0) //elipsoïd condition  
                            {
                                Phi[ii*Ny*Nz+jj*Nz+kk] = 1.0;
                                
                                U[ii*Ny*Nz+jj*Nz+kk] = U_R;
                            }
                        }
                    }
                }
            }
        }   
    }
}
