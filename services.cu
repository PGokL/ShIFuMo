void DisplayDeviceProperties(int Ndev)
{
    cudaDeviceProp deviceProp;
    memset( &deviceProp, 0, sizeof(deviceProp));
    if( cudaSuccess == cudaGetDeviceProperties(&deviceProp, Ndev))
    {
        printf( "#==============================================================");
        printf( "\n#Device Name \t %s ", deviceProp.name );
        printf( "\n#Device Index\t %d ", Ndev );
        printf( "\n#==============================================================");
        printf( "\n#Total Global Memory                  \t %ld KB", (long int)(deviceProp.totalGlobalMem/1024) );
        printf( "\n#Shared memory available per block    \t %ld KB", (long int)(deviceProp.sharedMemPerBlock/1024) );
        printf( "\n#Number of registers per thread block \t %d", deviceProp.regsPerBlock );
        printf( "\n#Warp size in threads             \t %d", deviceProp.warpSize );
        printf( "\n#Memory Pitch                     \t %ld bytes", (long int)(deviceProp.memPitch) );
        printf( "\n#Maximum threads per block        \t %d", deviceProp.maxThreadsPerBlock );
        printf( "\n#Maximum Thread Dimension (block) \t %d * %d * %d", deviceProp.maxThreadsDim[0], deviceProp.maxThreadsDim[1], deviceProp.maxThreadsDim[2] );
        printf( "\n#Maximum Thread Dimension (grid)  \t %d * %d * %d", deviceProp.maxGridSize[0], deviceProp.maxGridSize[1], deviceProp.maxGridSize[2] );
        printf( "\n#Total constant memory            \t %ld bytes", (long int)(deviceProp.totalConstMem) );
        printf( "\n#CUDA ver                         \t %d.%d", deviceProp.major, deviceProp.minor );
        printf( "\n#Clock rate                       \t %d KHz", deviceProp.clockRate );
        printf( "\n#Texture Alignment                \t %ld bytes", (long int)(deviceProp.textureAlignment) );
        printf( "\n#Device Overlap                   \t %s", deviceProp. deviceOverlap?"Allowed":"Not Allowed" );
        printf( "\n#Number of Multi processors       \t %d", deviceProp.multiProcessorCount );
        printf( "\n#==============================================================\n");
    }
    else
    {
        printf( "\n#Could not get properties for device %d.....\n", Ndev);
    }

}

// modulo for PBC
__device__ int mod(int a, int b)
{
    return a<0? b+a: a%b;
}

//void GetMemUsage(int *Array,int Num)
//{
//    int LENMAX=100;
//    char buffer[LENMAX];
//    std::string StrUse = "";
//    FILE* pipe = popen("nvidia-smi -q --display=MEMORY | grep Used ", "r");
//    while(!feof(pipe)) { if(fgets(buffer,LENMAX,pipe) != NULL){ StrUse+=buffer; } }
//    pclose(pipe);
//    for(int dev=0;dev<Num;dev++)
//    {
//        std::istringstream iss(StrUse.substr(StrUse.find(":")+1,StrUse.find("MB")-StrUse.find(":")-1));
//        iss >> Array[dev];
//        StrUse=StrUse.substr(StrUse.find("\n")+1,StrUse.length()-StrUse.find("\n")-1);
//    }
//}

//int GetFreeDevice(int Num)
//{
//    int FreeDev=-1;
//    int MemFree=15;
//    int *Memory_Use = new int[Num] ;
//    
//    // Check utilization of Devices
//    GetMemUsage(Memory_Use,Num);
//    // See if one is free
//    int dev=0;
//    do{
//        if(Memory_Use[dev]<MemFree)
//        { 
//            // Found one...
//            FreeDev=dev; 
//            // Check if it is really free...
//            system("sleep 1s");
//            GetMemUsage(Memory_Use,Num);
//            if(Memory_Use[dev]>MemFree){ FreeDev=-1; }
//            // twice...
//            system("sleep 1s");
//            GetMemUsage(Memory_Use,Num);
//            if(Memory_Use[dev]>MemFree){ FreeDev=-1; }
//        }    
//        dev++;
//    }while(FreeDev==-1 && dev<Num);

//    delete [] Memory_Use;
//    
//    if(FreeDev==-1)
//    {
//        printf("#=======================================\n");
//        system("nvidia-smi -q --display=MEMORY |grep U");
//        printf("#=======================================\n");
//        printf("#NO AVAILABLE GPU: SIMULATION ABORTED...\n");
//        printf("#=======================================\n\n");
//    }
//    return FreeDev;
//}


// test if the machine is working with little or big endian
int LittleEndian(void)
{       
    int num = 1;
    if(*(char *)&num == 1)
    {
        return 1;       //little endian
    }
    else
    {
        return 0;       // big endian
    }
}

float ReverseFloat( const float inFloat )
{
   float retVal;
   char *FloatToConvert = ( char* ) & inFloat;
   char *returnFloat = ( char* ) & retVal;

   // swap the bytes into a temporary buffer
   returnFloat[0] = FloatToConvert[3];
   returnFloat[1] = FloatToConvert[2];
   returnFloat[2] = FloatToConvert[1];
   returnFloat[3] = FloatToConvert[0];
   
   return retVal;
}
