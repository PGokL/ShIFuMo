// --------------------- //
// Initialization on GPU //
// --------------------- //
__global__ void InitializeField(REAL *U,  Parameters *par)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x,
        j = threadIdx.y + blockIdx.y * blockDim.y,
        k = threadIdx.z + blockIdx.z * blockDim.z;

    int Nx=par[0].Nx;
    int Ny=par[0].Ny;
    int Nz=par[0].Nz;

    int pos = Ny*Nz*i + Nz*j + k;

    // set initial composition to 1.0
    if (i<Nx && j<Ny && k<Nz)
    {
        U[pos] = par[0].Omega;
    }
}


// integrate diffusion equation with Euler method
__global__ void Diffusion(REAL *U, REAL *Utemp, REAL *Phi, Parameters *par)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    int j = threadIdx.y + blockIdx.y * blockDim.y;
    int k = threadIdx.z + blockIdx.z * blockDim.z;

    int Nx=par[0].Nx;
    int Ny=par[0].Ny;
    int Nz=par[0].Nz;
    int pos = Ny*Nz*i + Nz*j + k;
    int posN = Ny*Nz*i + Nz*mod((j+1),Ny) + k; //North
    int posS = Ny*Nz*i + Nz*mod((j-1),Ny) + k; //South
    int posE = Ny*Nz*mod(i+1,Nx) + Nz*j + k; //East
    int posW = Ny*Nz*mod(i-1,Nx) + Nz*j + k; //West
    int posA = Ny*Nz*i + Nz*j + mod(k+1,Nz); //Above
    int posU = Ny*Nz*i + Nz*j + mod(k-1,Nz); //Under

    REAL Pref = par[0].D*par[0].dt/(par[0].dx*par[0].dx);

    if(i<Nx && j<Ny && k<Nz)
    {
        Utemp[pos] = U[pos] + Pref*(1.0-Phi[pos])*(U[posN] + U[posS] + U[posE] + U[posW] + U[posA] + U[posU] - 6*U[pos]);
    }

}

