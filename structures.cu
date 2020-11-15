struct Parameters
{
    int conf;           // type of initial configuration

    REAL D;             // Diffusion coefficient
    REAL d0;            // capillary length (nm)
    REAL ceq_p;         // equilibrium precipitate concentation from Phase Diagram (atomic fraction)
    REAL ceq_m;         // equilibrium matrix concentation from Phase Diagram (atomic fraction)
    REAL cinit_m;       // initial concentation in the matrix (atomic fraction)
	REAL Omega;		    // supersaturation (dimension-less concentration)
    
    int Nsteps;         // number of time-steps
    int Nx, Ny, Nz;     // dimensions in real space
    int Npts;           // number of points
    int seed;           // seed to generate random numbers
    int NOutput;        // number of output files
    int every;          // dump infos only if we got past N time steps since last file 
 
    REAL Lx, Ly, Lz;    // length of the simulation cell
  
    REAL lc;            // characteristic length
    REAL tc;            // characteristic time
        
    REAL dx, dy, dz, dt;// discretisation parameters (dimension-less)
    
    int Nppts;          // number of precipitates in the domain
    char dist_name[200];// name of the distribution for the ppts radii
    REAL Rinit_1;       // parameter 1 for initialization of ppts (average radius)
    REAL Rinit_2;       // parameter 2 for initialization of ppts (std_dev radius)    
    REAL Rth;           // radius under which the precipitate disappears
    int Fppts;          //ppts shapes in the data.in file 
    REAL nu;            //Shape factor
    
};

struct Coord
{
    REAL x;
    REAL y;
	REAL z;
};

struct Precipitate
{
    int state;            // active=1 or not 
    REAL R	;             // radius
    int i;                // precipitate coordinates
    int j;
    int k;
    int form; 	         //Ppts type 0 = sphere, 1 = cylinder, 2 = elïpsoide 
    REAL nu;             //Shape factor 
    int hx;		        //integration box height
    int hy;
    int hz; 
    Coord box[8];       // integration box (not used)
};



