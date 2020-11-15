// defines
#define REAL double
#define BLOCK_SIZE_X 2
#define BLOCK_SIZE_Y 2
#define BLOCK_SIZE_Z 16
#define TIMECODE 1

// includes
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <cuda_runtime.h>
#include <curand_kernel.h>
#include <curand.h>
#include <cuda_profiler_api.h>
#include <iostream>

#include "services.cu"
#include "structures.cu"
#include "inout.cu"
#include "diffusion.cu"
#include "precipitates.cu"


//////////////////
// Program main //
//////////////////

// Host: CPU  --> variables with _h
// Device: GPU --> variables with _d


int main(int argc, char **argv)
{

    //----------------------------------------------------//
    // choose the device (GPU) and display its properties //
    //----------------------------------------------------//
    cudaSetDevice(atoi(argv[2]));
    DisplayDeviceProperties(atoi(argv[2]));


    //---------------------//
    // Read parameter file //
    //---------------------//
    Parameters par_h[1];
    ReadInputFile(argv[1], par_h);
    DisplayParams(par_h);

    Parameters *par_d;    
    cudaMalloc((void**)&par_d, sizeof(Parameters));
    cudaMemcpy(par_d, par_h, sizeof(Parameters), cudaMemcpyHostToDevice);

    //------------//
    // Memory CPU //
    //------------//
    size_t RealSize = par_h[0].Nx * par_h[0].Ny * par_h[0].Nz * sizeof(REAL);
    REAL *U_h =(REAL *) malloc(RealSize);
    REAL *Phi_h =(REAL *) malloc(RealSize);
    Precipitate *ppts_h = (Precipitate *) malloc(sizeof(Precipitate)*par_h[0].Nppts);
    
    //------------//
    // Memory GPU //
    //------------//
    REAL *U_d, *Utemp_d, *Ubuff_d;
    REAL *Phi_d;
    cudaMalloc((void**)&U_d, RealSize);
    cudaMalloc((void**)&Utemp_d, RealSize);
    cudaMalloc((void**)&Phi_d, RealSize);

    Precipitate *ppts_d;
    cudaMalloc((void**)&ppts_d,sizeof(Precipitate)*par_h[0].Nppts);

    // GPU blocks for diffusion
    dim3 DimBlocks(BLOCK_SIZE_X, BLOCK_SIZE_Y, BLOCK_SIZE_Z);
    dim3 NumBlocks(par_h[0].Nx/BLOCK_SIZE_X, par_h[0].Ny/BLOCK_SIZE_Y, par_h[0].Nz/BLOCK_SIZE_Z);
    
    // GPU blocks for ppts
    dim3 DimBlocksPpts(1);
    dim3 NumBlocksPpts(par_h[0].Nppts);


    //---------------------------//
    // set initial configuration //
    //---------------------------//
    
    srand(par_h[0].seed); // initialize random seed
    
    InitializeField<<<NumBlocks,DimBlocks>>>(U_d, par_d);
    cudaMemcpy(U_h, U_d, RealSize, cudaMemcpyDeviceToHost);

    InitializePrecipitates(ppts_h, par_h);
    cudaMemcpy(ppts_d, ppts_h, sizeof(Precipitate)*par_h[0].Nppts, cudaMemcpyHostToDevice);

    UpdateFields<<<NumBlocksPpts,DimBlocksPpts>>>(U_d, Phi_d, ppts_d, par_d);
    cudaMemcpy(Phi_h, Phi_d, RealSize, cudaMemcpyDeviceToHost);
    cudaMemcpy(U_h, U_d, RealSize, cudaMemcpyDeviceToHost);
    cudaMemcpy(ppts_h, ppts_d, sizeof(Precipitate)*par_h[0].Nppts, cudaMemcpyDeviceToHost);


    //----------------------//
    // Output initial state //
    //----------------------//
    int index=0;
    int iter=0;
    WriteUVtk(index, U_h, par_h, iter);
    WritePhiVtk(index, Phi_h, par_h, iter);
    WritePpts(index, ppts_h, par_h, iter);
    WriteProps(par_h, U_h, Phi_h, ppts_h, iter, 1);
    
    
    //-----------------//
    // timer variables //
    //-----------------//
#if(TIMECODE)
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    float timer;
    cudaEventRecord(start, 0);
#endif


    for(iter=1; iter<=par_h[0].Nsteps; iter++) 
    {    

        //-----------//
        // Diffusion //
        //-----------//  
        Diffusion<<<NumBlocks,DimBlocks>>>(U_d, Utemp_d, Phi_d, par_d);
        Ubuff_d=U_d;     U_d=Utemp_d;      Utemp_d=Ubuff_d;
        
        //------------//
        // PPT growth //
        //------------//  
        PrecipitatesGrowth<<<NumBlocksPpts,DimBlocksPpts>>>(U_d, Phi_d, ppts_d, par_d);
        
        //-------------//
        // PPT overlap //
        //-------------// 
        PrecipitatesOverlap<<<NumBlocksPpts,DimBlocksPpts>>>(U_d, Phi_d, ppts_d, par_d);
        
        //------------//
        // PPT upadate //
        //------------//        
        UpdateFields<<<NumBlocksPpts,DimBlocksPpts>>>(U_d, Phi_d, ppts_d, par_d);
        
        
        //--------//
        // Output //
        //--------//  
        
        if((iter%par_h[0].every)==0)
        {
            cudaMemcpy(U_h, U_d, RealSize, cudaMemcpyDeviceToHost);
            cudaMemcpy(Phi_h, Phi_d, RealSize, cudaMemcpyDeviceToHost);
            cudaMemcpy(ppts_h, ppts_d, sizeof(Precipitate)*par_h[0].Nppts, cudaMemcpyDeviceToHost);

            WriteProps(par_h, U_h, Phi_h, ppts_h, iter, 0);            
            
            if((iter%(par_h[0].Nsteps/par_h[0].NOutput))==0)
            {
                index=index+1;
                WriteUVtk(index, U_h, par_h, iter);
                WritePhiVtk(index, Phi_h, par_h, iter);
                WritePpts(index, ppts_h, par_h, iter);
            } 
        }
        

    }
    
    
    cudaError_t err = cudaGetLastError();  // add
    if (err != cudaSuccess) std::cout << "CUDA error: " << cudaGetErrorString(err) << std::endl; // add
    cudaProfilerStop();
    
    #if(TIMECODE)
        cudaEventRecord(stop, 0); 
        cudaEventSynchronize(stop); 
        cudaEventElapsedTime(&timer, start, stop);
        timer=timer/1000.0;
        FILE *OutFile=fopen("timer.out","w");    
        fprintf(OutFile, "#total_time_sec\t #time_perGP_perTS_sec\n");
        fprintf(OutFile, "%.8g \t %.8g\n",timer,(timer/(REAL(par_h[0].Nx*par_h[0].Ny))/(REAL(par_h[0].Nsteps))));
        fclose(OutFile);
        printf("time=%f sec\n",timer);
    #endif
                                    
    return 0;

}

