void ReadInputFile(const char *fname, Parameters *par)
{
    FILE *ifp ;
    ifp = fopen(fname, "r"); // input file name "xxxx.in"
    REAL pi=2.*acos(0.);

    char info[200], units[200], line[200];
    double temp;
    double temp1, temp2, temp3;

    // Error if the command used isn't "./Diff3.x xxxx.in 0"
    if (ifp == NULL)  
    {
        fprintf(stderr, "Can't open input data file %s\n",fname);
        exit(1);
    }

    // restart configuration
    for(int i=0;i<3;i++)
    {   
        fscanf(ifp, "%s\n",line);
    }
    fscanf(ifp, "%s %d %s\n", info, &par[0].conf, units);
    fscanf(ifp, "%s %d %s\n", info, &par[0].seed, units);

    // materials parameters
    for(int i=0;i<3;i++)
    {   
        fscanf(ifp, "%s\n",line);
    }
    fscanf(ifp, "%s %lf %s\n", info, &temp, units);
    par[0].D=REAL(temp);
    fscanf(ifp, "%s %lf %s\n", info, &temp, units);
    par[0].d0=REAL(temp);
    fscanf(ifp, "%s %lf %s\n", info, &temp, units);
    par[0].ceq_p=REAL(temp);
    fscanf(ifp, "%s %lf %s\n", info, &temp, units);
    par[0].ceq_m=REAL(temp);
    fscanf(ifp, "%s %lf %s\n", info, &temp, units);
    par[0].cinit_m=REAL(temp);
    


    // read simulation parameters part
    for(int i=0;i<3;i++)
    {
        fscanf(ifp, "%s\n",line);
    }
    fscanf(ifp, "%s %lf %lf %lf %s\n", info, &temp1, &temp2, &temp3, units);
    par[0].Lx=REAL(temp1);
    par[0].Ly=REAL(temp2);
    par[0].Lz=REAL(temp3);
    fscanf(ifp, "%s %d %d %d %s\n", info, &(par[0].Nx), &(par[0].Ny), &(par[0].Nz), units);
    fscanf(ifp, "%s %d %s\n", info, &(par[0].Nsteps), units);
    fscanf(ifp, "%s %lf %s\n", info, &temp, units);
    par[0].dt=REAL(temp);
    fscanf(ifp, "%s %d %s\n", info, &(par[0].NOutput), units);
    fscanf(ifp, "%s %d %s\n", info, &(par[0].every), units);
    
    // read ppt distribution part
    for(int i=0;i<3;i++)
    {
        fscanf(ifp, "%s\n",line);
    }
    fscanf(ifp, "%s %d %s\n", info, &(par[0].Nppts), units);
    fscanf(ifp, "%s %s %s\n", info, par[0].dist_name, units);
    fscanf(ifp, "%s %lf %lf %s\n", info, &temp1, &temp2, units);
    par[0].Rinit_1=REAL(temp1);
    par[0].Rinit_2=REAL(temp2);
    fscanf(ifp, "%s %d %s\n", info, &(par[0].Fppts), units);
    fscanf(ifp, "%s %lf %s\n", info, &(par[0].nu), units);
    
    // compute some variables
    par[0].Npts = par[0].Nx*par[0].Ny*par[0].Nz;
    par[0].dx = par[0].Lx/REAL(par[0].Nx);
    par[0].dy = par[0].Ly/REAL(par[0].Ny);
    par[0].dz = par[0].Lz/REAL(par[0].Nz);

    // If matrix concentration is negative, then it will be set to an equilibrium concentration (Ostwald ripening)
    if(par[0].cinit_m<0.0)
    {
        par[0].cinit_m = par[0].ceq_m*exp(par[0].d0/par[0].Rinit_1);
    }
    par[0].Omega = (par[0].cinit_m - par[0].ceq_m)/(par[0].ceq_p - par[0].ceq_m); //Omega
    
    par[0].Rth = 1.1*par[0].d0/log(par[0].ceq_p/par[0].ceq_m); // Dissolve radius

    // set characteristic values
    par[0].lc=par[0].dx;                     // nm
    par[0].tc=par[0].lc*par[0].lc/par[0].D;  // sec
    
    // normalize values
    par[0].D=par[0].D/(par[0].lc*par[0].lc)*par[0].tc;
    par[0].d0=par[0].d0/par[0].lc;
    par[0].Rinit_1=par[0].Rinit_1/par[0].lc;
    par[0].Rinit_2=par[0].Rinit_2/par[0].lc;
    par[0].Rth=par[0].Rth/par[0].lc;
    par[0].Lx=par[0].Lx/par[0].lc;
    par[0].Ly=par[0].Ly/par[0].lc;
    par[0].Lz=par[0].Lz/par[0].lc;
    par[0].dx=par[0].dx/par[0].lc;
    par[0].dy=par[0].dy/par[0].lc;
    par[0].dz=par[0].dz/par[0].lc;
    
    if(par[0].seed<=0)
    {
        par[0].seed=time(NULL);
    }
}


void DisplayParams(Parameters *par)
{  
    printf("##############################################\n");
    printf("#    ShIFuMo 2020 Geslin Mizrahi Gokelaere   #\n");
    printf("##############################################\n");
    printf("\n");
    printf("#--------------------------------#\n");
    printf("#     SIMULATION PARAMETERS      #\n");
    printf("#--------------------------------#\n"); 
    printf("# Lx    = %12.6g    nm     #\n", par[0].Lx*par[0].lc);
    printf("# Ly    = %12.6g    nm     #\n", par[0].Ly*par[0].lc);
    printf("# Lz    = %12.6g    nm     #\n", par[0].Lz*par[0].lc);
    printf("#--------------------------------#\n");
    printf("# dx    = %12.6g    nm     #\n", par[0].dx*par[0].lc);
    printf("# dy    = %12.6g    nm     #\n", par[0].dy*par[0].lc);
    printf("# dz    = %12.6g    nm     #\n", par[0].dz*par[0].lc);
    printf("#--------------------------------#\n");
    printf("# Nsteps= %12d    steps  #\n", par[0].Nsteps);
    printf("# dt    = %12.6g    sec    #\n", par[0].dt*par[0].tc);
    printf("# Time  = %12.6g    sec    #\n", par[0].dt*par[0].tc*par[0].Nsteps);
    printf("# seed  = %12d           #\n",par[0].seed);
    printf("# Shape  = %12d          #\n",par[0].Fppts);
    printf("# Shape factor  = %12.6f         #\n",par[0].nu);
    printf("#--------------------------------#\n");
    printf("#     CHARACTERISTIC VALUES      #\n");
    printf("#--------------------------------#\n"); 
    printf("# lc    = %12.6g    nm     #\n", par[0].lc);
    printf("# tc    = %12.6g    sec    #\n", par[0].tc);
    printf("#--------------------------------#\n");
}

//-----------------//
// Writing u-field //
//-----------------//
void WriteUVtk(int index, REAL *U, Parameters *par, int iter)
{
    char OutFileName[200];

    sprintf(OutFileName,"u_%06i.vtk",index);
    FILE * OutFile=fopen(OutFileName, "wb");
    
    fprintf(OutFile,"# vtk DataFile Version 2.0\n");
    fprintf(OutFile,"u t=%.8g sec\n",par[0].dt*iter*par[0].tc);
    fprintf(OutFile,"BINARY\n");
    fprintf(OutFile,"DATASET STRUCTURED_POINTS\n");
    fprintf(OutFile,"ORIGIN \t %f %f %f\n",0.,0.,0.);
    fprintf(OutFile,"SPACING \t %f %f %f\n",par[0].dx*par[0].lc,par[0].dy*par[0].lc,par[0].dz*par[0].lc);
    fprintf(OutFile,"DIMENSIONS \t %i %i %i\n",par[0].Nx,par[0].Ny,par[0].Nz);
    fprintf(OutFile,"POINT_DATA \t %i\n",par[0].Nx*par[0].Ny*par[0].Nz);
    fprintf(OutFile,"SCALARS \t volume_scalars float 1\n");
    fprintf(OutFile,"LOOKUP_TABLE \t default\n");
    
    if(LittleEndian())        // swap the bits values before writing
    {
        float temp;
        for(int i=0;i<par[0].Nx;i++)
        {
            for(int j=0;j<par[0].Ny;j++)
            {
                for(int k=0;k<par[0].Nz;k++)
                {
                    temp=ReverseFloat(float(U[i*par[0].Ny*par[0].Nz + par[0].Nz*j + k]));
                    fwrite(&temp, sizeof(float), 1, OutFile);
                }
            }
        }
    }
    else
    {
        fwrite(&U, sizeof(REAL), par[0].Nx*par[0].Ny*par[0].Nz, OutFile);
    }
    
    fclose(OutFile);
    
}

//-------------------//
// Writing phi-field //
//-------------------//
void WritePhiVtk(int index, REAL *Phi, Parameters *par, int iter)
{   
    char OutFileName[200];

    sprintf(OutFileName,"phi_%06i.vtk",index);
    FILE *OutFile=fopen(OutFileName, "wb");
    
    fprintf(OutFile,"# vtk DataFile Version 2.0\n");
    fprintf(OutFile,"phi t=%12.6g sec\n",par[0].dt*iter*par[0].tc);
    fprintf(OutFile,"BINARY\n");
    fprintf(OutFile,"DATASET STRUCTURED_POINTS\n");
    fprintf(OutFile,"ORIGIN \t %f %f %f\n",0.,0.,0.);
    fprintf(OutFile,"SPACING \t %f %f %f\n",par[0].dx*par[0].lc,par[0].dy*par[0].lc,par[0].dz*par[0].lc);
    fprintf(OutFile,"DIMENSIONS \t %i %i %i\n",par[0].Nx,par[0].Ny,par[0].Nz);
    fprintf(OutFile,"POINT_DATA \t %i\n",par[0].Nx*par[0].Ny*par[0].Nz);
    fprintf(OutFile,"SCALARS \t volume_scalars float 1\n");
    fprintf(OutFile,"LOOKUP_TABLE \t default\n");
    
    if(LittleEndian())        // if the system works with little endian, swap the bits values before writing
    {
        float temp;
        for(int i=0;i<par[0].Nx;i++)
        {
            for(int j=0;j<par[0].Ny;j++)
            {
                for(int k=0;k<par[0].Nz;k++)
                {
                    temp=ReverseFloat(float(Phi[i*par[0].Ny*par[0].Nz + par[0].Nz*j + k]));
                    fwrite(&temp, sizeof(float), 1, OutFile);
                }
            }
        }
    }

    fclose(OutFile);
}

//--------------------//
// Writing ppts infos //
//--------------------//
void WritePpts(int index, Precipitate *ppts, Parameters *par, int iter)
{
    char OutFileName[200];

    sprintf(OutFileName,"ppts_%06i.out",index);
    FILE * OutFile=fopen(OutFileName, "wb");
    

    for(int n=0;n<par[0].Nppts;n++)
    {
        fprintf(OutFile, "%d %d %d %d %d %f\n", n, ppts[n].state, ppts[n].i, ppts[n].j, ppts[n].k, ppts[n].R*par[0].lc);
    }
    
    fclose(OutFile);
    
}


//-----------------------------//
// Writing <R(t)> en <U(t)> data //
//-----------------------------//
void WriteRfile(REAL *Rf, REAL *mu, REAL *t, int index)
{
    char OutFileName[200];
    
    sprintf(OutFileName,"R_t.out");
    FILE * OutFile=fopen(OutFileName, "wb");
    
    for(int n=0; n<index;n++)
    {
        fprintf(OutFile, "%d %f %f %f\n", n, Rf[n], t[n], mu[n]);
    }
    
    fclose(OutFile);

}




//----------------------------//
// Writing average properties //
//----------------------------//
void WriteProps(Parameters *par, REAL *U, REAL *Phi, Precipitate *ppts, int iter, int firstTime)
{
    if(firstTime)
    {
        FILE *OutFile=fopen("props.out","w");    
        fprintf(OutFile, "#iter#\t#time_sec#\t#Nppts#\t#R_ave_nm#\t#f_ppt#\t#U_ave#\t#C_ave#\n");
        fclose(OutFile);
    }
    
    int Nppts = par[0].Nppts;
    int Npts = par[0].Npts;
        
    // compute average precipitate radius
    REAL Rave = 0;
    int Ntot=0;

    for(int n=0; n<Nppts; n++)
    {
        if(ppts[n].state == 1)
        {
             Ntot = Ntot+1;
             Rave = Rave + ppts[n].R;
        }
    }
    Rave= Rave/REAL(Ntot);
        
        
    // compute phase fraction of ppts and number of active ppts
    REAL fppt = 0;
    int Nppts_ac = 0;

    for(int n=0; n<Nppts; n++)
    {
        if(ppts[n].state == 1)
        {
            fppt += 4./3.*M_PI*(ppts[n].R*ppts[n].R*ppts[n].R);
            Nppts_ac += 1;
        }
    }
    fppt = fppt/REAL(Npts);
    
    
    // compute total average concentration
    REAL Uave = 0;
    REAL Cave = 0;
        
    for(int n=0; n<Npts;n++)
    {
        Uave = Uave+U[n]*(1.0-Phi[n]); // Concentration in matrix
    }
    
    for(int n=0; n<Nppts;n++)
    {
        if(ppts[n].state==1)
        {
            Uave = Uave + (4./3.)*M_PI*ppts[n].R*ppts[n].R*ppts[n].R;
        }            
    }
    Uave = Uave/REAL(Npts); // average dimension-less concentration in the system
    Cave = par[0].ceq_m + (par[0].ceq_p - par[0].ceq_m)*Uave; // average concentration in the system
    
    FILE *OutFile=fopen("props.out","a");    
    fprintf(OutFile, "%d \t %.8g \t %d \t %.8g \t %.8g \t %.8g \t %.8g  \n", 
                          iter, iter*par[0].dt*par[0].tc, Nppts_ac, Rave*par[0].lc, fppt, Uave, Cave);
    fclose(OutFile);
}
