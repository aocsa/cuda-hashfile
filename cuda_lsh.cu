
#include <iostream>
#include <fstream>
#include <ctime>
#include <sstream>
#include <vector>
#include <map>
using namespace std;

#include <curand_kernel.h>
#include <cuda.h>
#include <cstdlib>
#include <thrust/device_ptr.h>

#include "ResultSet.h"
#include "DistanceEvaluator.h"
#include "cuPrintf.cu"

#include "bitonic_kernel.h"

typedef int tObject ; //only point id
typedef ResultSet<tObject> tResult;
#define MAX_FLOAT 3.40282347E+7F
// #define N_THREADS 512
// 
// #define MAX_TOPK 1024
// #define MAX_N_LEVELS 10
// #define LAMBDA_STEP 1000


__global__
void cuda_init_curand(curandState* states) {
    int index = blockIdx.x;
    curand_init(index, index, 0, &states[index]);
}
  

__global__
void initHashFunctionskernel/*<<<nLevels, 1>>>*/(int *hashfunctions, int dim,
												 curandState* states) {
    int index = blockIdx.x;
    curandState state = states[index];
    for (int i = 0; i < dim; i++)
		hashfunctions[blockIdx.x * dim + i] = curand(&state) % 2; 
}


__global__
void hashDatasetkernel(float *dataset, int N, int UpperSize, int dim,
					   int *hashfunctions, int nLevels, float *keys, uint* ids) {

    uint ipoint = blockIdx.x * blockDim.x + threadIdx.x;
    if (ipoint < N) {
        //dataset[ ipoint * dim ] ; // point [i]
        for (int i = 0; i < nLevels; i++) {
            float h = 0.0;
            for (int j = 0; j < dim; j++)
                h +=  ((2 * hashfunctions[ i*dim + j ] - 1) * dataset[ ipoint * dim  + j] );    
            
            h += h*10 + MAX_FLOAT;

            keys[i * UpperSize  + ipoint] = h; 
            ids[i * UpperSize + ipoint] = ipoint;         
        }
    }
    else if (ipoint < UpperSize ) {
        for (int i = 0; i < nLevels; i++) {
            keys[i * UpperSize  + ipoint] = MAX_FLOAT*1000; 
            ids[i * UpperSize + ipoint] = ipoint;         
        }
    }
} 

__device__
int binary_search(float * container, int first, int last, float object)
{
	if( first >= last )
		return first;
	while( first < last )
	{
	   int mid = (first+last)/2;
	   if( object == container[mid ] )
			   return mid;
	   if( object > container[mid ] )
			   first = mid+1;
	   else
			   last  = mid;
	}
	if( object <= container[first] )
		return first;
	return last;
}



__device__
void merge(int* Avalues, int *Bvalues, float* A, float *B, int N, int left, int right) {

    for(int i = 0;i < N; i++)
    {
        if ( A[left]  < B[right] ) {
            A[i] = A[left];
            Avalues[i] = Avalues[left];
            
            left++;
            right++;
        }
        else {
             A[i] = B[right];
             Avalues[i] = Bvalues[left];
             right++;
             left++;
        }
    }
}


// distances initialized with \infty
__device__ float add_query_result(int obj_id, float d,  float *distances, int *result_set, int N)
{
    int pos = 0;
    while(distances[pos] < d && pos < N)
        pos++;

    if ( pos < N ) {
        if (result_set[pos] != obj_id) {
            for (int i = N - 1; i >= pos; i--) {
                distances[i+1] = distances[i];
                result_set[i+1] = result_set[i];
            }
            distances[pos] = d;
            result_set[pos] = obj_id;
        }
    }
    return distances[N - 1];
}


__host__ __device__
float distance(float* A, float* B, int size) {
    float d = 0.0;
    for (int k = 0; k < size; k++) {
        d += fabs(A[k] - B[k]);
    }
    return d;
}


__global__
void TopKMotifHashFileKernel(float* keys, uint* ids, int *hashfunctions,
							 float* dataset, int N, int UpperSize, int dim,
							 int nLevels, int topK, int* out_result_set,
							 float*  best_so_far) {

    uint ipoint = blockIdx.x * blockDim.x + threadIdx.x;
    if (ipoint < N) {
        float *query =  &dataset[ ipoint * dim ]; 
        
        float   distances[MAX_TOPK];
        int     result_set[MAX_TOPK];

        for (int i = 0; i < topK; i++) {
            distances[i] = MAX_FLOAT;
        }
        for (int ilevel = 0; ilevel < nLevels; ilevel++) {
            float h = 0.0;
            for (int j = 0; j < dim; j++)
                h +=  ((2 * hashfunctions[ ilevel*dim + j ] - 1) * query[j] );    
            h += h*10 + MAX_FLOAT;
    
            int pos = binary_search(keys, 0, N, h); 

            int BEGIN = pos - LAMBDA_STEP;
            int END = pos + LAMBDA_STEP;
            
            float bsf = MAX_FLOAT;
            for (int j = BEGIN; j <= END; j++) {
                if (j >= 0 && j < N) {
                    int obj_id = ids[ilevel*UpperSize + j];
                    float d = distance(&dataset[ obj_id * dim ], query, dim);
                    add_query_result(obj_id, d, distances, result_set, topK);
                }  
            }
            float dist = distances[topK - 1];
            float prev_bsf = __int_as_float( atomicMin( (int*) best_so_far, __float_as_int(dist) ) );
            if( dist < prev_bsf ) {
                atomicExch(&out_result_set[0], ipoint);
            }
        }
    }
}


__global__
void AllTopKMotifHashFileKernel(float * keys, uint * ids,  int *hashfunctions,
								float * dataset, int N, int UpperSize, int dim,
								int nLevels, int topK, uint * out_ids, 
								float*  out_distances) {

    uint ipoint = blockIdx.x * blockDim.x + threadIdx.x;
    if (ipoint < N) {
        float *query =  &dataset[ ipoint * dim ]; 
        
        float   distances[MAX_TOPK];
        int     result_set[MAX_TOPK];

        for (int i = 0; i < topK; i++) {
            distances[i] = MAX_FLOAT;
        }
        for (int ilevel = 0; ilevel < nLevels; ilevel++) {
            float h = 0.0;
            for (int j = 0; j < dim; j++)
                h +=  ((2 * hashfunctions[ ilevel*dim + j ] - 1) * query[j] );    
            h += h*10 + MAX_FLOAT;
    
            int pos = binary_search(keys, 0, N, h); 

            int BEGIN = pos - LAMBDA_STEP;
            int END = pos + LAMBDA_STEP;
            
            float bsf = MAX_FLOAT;
            for (int j = BEGIN; j <= END; j++) {
                if (j >= 0 && j < N) {
                    int obj_id = ids[ilevel*UpperSize + j];
                    float d = distance(&dataset[ obj_id * dim ], query, dim);
                    add_query_result(obj_id, d, distances, result_set, topK);
                }  
            }
            out_distances[ipoint] = distances[topK - 1];
            out_ids[ipoint] = ipoint;
        }
    }
    else if (ipoint < UpperSize ) {
        out_distances[ipoint] = MAX_FLOAT;
        out_ids[ipoint] = ipoint;      
    }
}


__global__ 
void NearestQueryHashFileKernel(float *keys, uint *ids, int *hashfunctions,
								float *dataset, int N, int UpperSize, int dim,
								int nLevels, float *query, int topK,
								int *out_result, float *out_distances) {

    for (int i = 0; i < topK; i++) {
        out_distances[i] = MAX_FLOAT;
    }
 
    for (int ilevel = 0; ilevel < nLevels; ilevel++) {
        float h = 0.0;
        for (int j = 0; j < dim; j++) {
            h += ((2 * hashfunctions[ilevel*dim + j ] - 1) * query[j]);    
        }
        h += h*10 + MAX_FLOAT;

        int pos = binary_search(keys, 0, N, h);

        int BEGIN = pos - LAMBDA_STEP;
        int END = pos + LAMBDA_STEP;

        for (int j = BEGIN; j <= END; j++) {
            if (j >= 0 && j < N) {
                int obj_id = ids[ilevel*UpperSize + j];
                float d = distance(&dataset[ obj_id * dim ], query, dim);
                add_query_result(obj_id, d, out_distances, out_result, topK);
            }  
        }
    }
}

__global__
void NearestQueryHashFileKernel_alllevels(float *keys, uint *ids,
										  int *hashfunctions, float *dataset,
										  int N, int UpperSize, int dim,
										  int nLevels, float *query, int topK,
										  int *out_result, float *out_distances) {
    
    uint ilevel = blockIdx.x * blockDim.x + threadIdx.x;

    for (int i = 0; i < topK; i++) {
        out_distances[i] = MAX_FLOAT;
    }
 
    for (int ilevel = 0; ilevel < nLevels; ilevel++) {
        float h = 0.0;
        for (int j = 0; j < dim; j++) {
            h += ((2 * hashfunctions[ilevel*dim + j ] - 1) * query[j]);    
        }
        h += h*10 + MAX_FLOAT;

        int pos = binary_search(keys, 0, N, h);

        int BEGIN = pos - LAMBDA_STEP;
        int END = pos + LAMBDA_STEP;

        for (int j = BEGIN; j <= END; j++) {
            if (j >= 0 && j < N) {
                int obj_id = ids[ilevel*UpperSize + j];
                float d = distance(&dataset[ obj_id * dim ], query, dim);
                int pos = 0;
                while(out_distances[pos] < d && pos < N)
                    pos++;

                if ( pos < N ) {
                    if (out_result[pos] != obj_id) {
                        for (int i = N - 1; i >= pos; i--) {
                            atomicExch(&out_distances[i+1], out_distances[i]);
                            atomicExch(&out_result[i+1], out_result[i]);
                        }
                        atomicExch(&out_distances[pos], d);
                        atomicExch(&out_result[pos], obj_id);
                    }
                }
            }  
        }
    }
     
}


class CUDAHashFile {
    curandState* d_states;  
    int*         d_hashfunctions; 
    float*       d_dataset;
    float*       d_keysContainer;
    uint *       d_idsContainer;
    float*       d_query; 
    int*         d_queryResult;
    float*       d_distances;

    int*         h_result;
    float*       h_distances;

    int N;
    int UpperSize;
    int dim;
    int nLevels;
    MetricEvaluator* metric;
    float* dataset;

public:
    CUDAHashFile(MetricEvaluator* _metric, float *data, int _N,  int _d,  int _nlevels) 
    {
        metric = _metric;
        dataset = data;
        dim = _d;
        nLevels = _nlevels;

        int limits [] = {1024, 2048, 4096, 8192, 16384, 32768, 65536, 131072, 262144, 524288, 1048576};

        N = _N;
        int i=0;
        while (limits[i] < N && i < sizeof(limits) / sizeof(limits[0]) )
            i++;
        
        UpperSize = limits[i];

        int size2 = N * dim * sizeof(float);
        cutilSafeCall( cudaMalloc(&d_states, dim * nLevels * sizeof(curandState)) );
        cutilSafeCall( cudaMalloc((void **)&d_keysContainer,  UpperSize * nLevels * sizeof(float)) );
        cutilSafeCall( cudaMalloc((void **)&d_idsContainer,  UpperSize * nLevels * sizeof(uint)) );
        cutilSafeCall( cudaMalloc((void **)&d_hashfunctions,  N * dim * sizeof(int)) );
        cutilSafeCall( cudaMalloc(&d_query, sizeof(float) * dim));
        cutilSafeCall( cudaMalloc(&d_queryResult, sizeof(int) * MAX_TOPK) );
        cutilSafeCall( cudaMalloc(&d_distances, sizeof(float) * MAX_TOPK) );

        h_result =  new int[MAX_TOPK];
        h_distances = new float[ MAX_TOPK ];

        cutilSafeCall(cudaMalloc(&d_dataset, size2));
        cutilSafeCall(cudaMemcpy(d_dataset, data, size2, cudaMemcpyHostToDevice));  
    
        cuda_init_curand<<<nLevels, 1 >>>(d_states);
        initHashFunctionskernel<<<nLevels, 1>>>(d_hashfunctions, dim, d_states);
    }

    ~CUDAHashFile() {
        cudaFree(d_states);
        cudaFree(d_hashfunctions);
        cudaFree(d_dataset);
        cudaFree(d_keysContainer);
        cudaFree(d_idsContainer);
        cudaFree(d_query);
        cudaFree(d_queryResult);
        cudaFree(d_distances);
     
        delete []h_result;
        delete []h_distances;   
    }

    MetricEvaluator * GetEvaluatorType() {
        return metric;
    }

    virtual void BulkLoad() {

        float *h_InputKey, *d_InputKey,  *h_OutputKeyGPU, *d_OutputKey;
        uint  *d_InputVal, *d_OutputVal, *h_InputVal,  *h_OutputValGPU;

        int  blocks = UpperSize/ N_THREADS;
        int  threads = N_THREADS;
        hashDatasetkernel <<<blocks, threads>>>(d_dataset, N, UpperSize, dim, d_hashfunctions, nLevels, d_keysContainer, d_idsContainer);
       

        /*float *h_ptr  = new float[ UpperSize * nLevels];
        cudaMemcpy(h_ptr, d_keysContainer, sizeof(float) * UpperSize * nLevels, cudaMemcpyDeviceToHost);

        for (int i = 0 ; i < nLevels; i++) {
            for (int j = 0 ; j < UpperSize; j++){
                printf("%d \t %f\n", j, h_ptr[i*UpperSize + j] );
            }
            printf("\n\n");
        }   

        printf("sorting...\n\n");*/

        for (int i = 0 ; i < nLevels; i++) {
        
            d_InputKey = d_keysContainer;
            d_InputVal = d_idsContainer;

            d_OutputKey = d_InputKey;
            d_OutputVal = d_InputVal;
            uint offset = i;

            bitonicSort(
                d_OutputKey,
                d_OutputVal,
                d_InputKey,
                d_InputVal,
                1,
                offset,
                UpperSize,
                1
            );  
        }
        /*
        cudaMemcpy(h_ptr, d_OutputKey, sizeof(float) * UpperSize * nLevels, cudaMemcpyDeviceToHost);

        for (int i = 0 ; i < nLevels; i++) {
            for (int j = 0 ; j < UpperSize; j++){
                printf("%d \t %f\n", j, h_ptr[i*UpperSize + j] );
            }
            printf("\n\n");
        }*/
    }


    virtual  tResult TopKMotif(int k) {
        int  blocks = (N + N_THREADS - 1)/ N_THREADS;
        int  threads = N_THREADS;

        float bsf = MAX_FLOAT;
        float* d_distance_bsf;
        cudaMalloc( &d_distance_bsf, sizeof(float) );
        cudaMemcpy( d_distance_bsf, &bsf, sizeof(float), cudaMemcpyHostToDevice );

        TopKMotifHashFileKernel<<<blocks, threads>>>( d_keysContainer, d_idsContainer,  d_hashfunctions, d_dataset, N, UpperSize, dim, nLevels, k, d_queryResult, d_distance_bsf);
        cudaMemcpy(h_result, d_queryResult, sizeof(int), cudaMemcpyDeviceToHost);
        cudaMemcpy(&bsf, d_distance_bsf, sizeof(float), cudaMemcpyDeviceToHost );

        int motif_id = h_result[0];
        cout << "motif_id: " << motif_id << endl;
        return NearestQuery(&dataset[motif_id*dim], k);
    }


    void AllTopKMotifs(int nMotifs, int topK) {        
        uint  *d_ids;
        uint* h_ids = new uint[nMotifs];
        
        float* d_distance_bsf;
        cutilSafeCall( cudaMalloc(&d_distance_bsf, UpperSize * sizeof(float) ));
        cutilSafeCall( cudaMalloc(&d_ids,  UpperSize * sizeof(uint)) );

        int  blocks = UpperSize / N_THREADS;
        int  threads = N_THREADS;
        AllTopKMotifHashFileKernel<<<blocks, threads>>>( d_keysContainer, d_idsContainer,  d_hashfunctions, d_dataset, N, UpperSize, dim, nLevels, topK, d_ids, d_distance_bsf);
        
        float *h_InputKey, *d_InputKey,  *h_OutputKeyGPU, *d_OutputKey;
        uint  *d_InputVal, *d_OutputVal, *h_InputVal,  *h_OutputValGPU;
        d_InputKey = d_distance_bsf;
        d_InputVal = d_ids;
        d_OutputKey = d_InputKey;
        d_OutputVal = d_InputVal;
        cout << "sorting " << endl;
        bitonicSort(
            d_OutputKey,
            d_OutputVal,
            d_InputKey,
            d_InputVal,
            1,
            0,
            UpperSize,
            1
        );  

        cudaMemcpy(h_ids, d_ids, nMotifs * sizeof(int), cudaMemcpyDeviceToHost );
        
        for (int i = 0; i < nMotifs; i++) {
            int motif_id = h_ids[i];
            cout << "motif_id: " << motif_id << endl;
            tResult result = NearestQuery(&dataset[motif_id*dim], topK);
            printf("Motif found %d NNs at distance %0.6lf. NNs are:\n", result.GetNumOfEntries(), result.GetMaximumDistance());
            for(int j = 0; j < result.GetNumOfEntries(); j++){
                printf("%09d\tdist:%0.6lf\n", result[j].GetObject(), result[j].GetDistance());
            }
        }
    }

     virtual tResult NearestQuery(float* query, int k) {
        tResult resultSet;
        resultSet.SetQueryInfo(0, tResult::KNEARESTQUERY, k, MAX_FLOAT, true);
        cudaMemcpy(d_query, query, sizeof(float) * dim, cudaMemcpyHostToDevice);

        cudaPrintfInit();
        NearestQueryHashFileKernel<<<1, 1>>>( d_keysContainer, d_idsContainer,  d_hashfunctions, d_dataset, N, UpperSize, dim, nLevels, d_query, k, d_queryResult, d_distances); // UpperSize to make id * N  + offset
        cudaPrintfDisplay();
        cudaPrintfEnd();

        cudaMemcpy(h_result, d_queryResult, sizeof(int) * k, cudaMemcpyDeviceToHost);
        cudaMemcpy(h_distances, d_distances, sizeof(float) * k, cudaMemcpyDeviceToHost);

        for (int i = 0; i < k; i++) {
            int obj_id = h_result[i];
            float dist = h_distances[i];
            resultSet.AddPair(obj_id, dist);
        }
        return resultSet;
    }
};


void load_data(std::string datafilepath, float * data, int N, int dim){
    std::ifstream ifs (datafilepath.c_str()); 
    std::string line;
    float value;

    for (int i = 0; i < N; i++) {
        getline(ifs, line);
        std::stringstream ss (line);
        for (int j = 0; j < dim; j++) {
            ss >> value;
            data[i * dim + j] = value;
        }
    }
}

void usage(char *programName){
    printf("Usage: %s #pts_in_data_set #queries dimension", programName);
	printf("successProbability #neighbors data_set_file queries_file ");
	printf("distaceType[L1=0|L2=1] nLevels\n");
}


int main(int nargs, char **args){
    if (nargs < 9) {
        usage(args[0]);
        exit(1);
    }
    int nPoints = atoi(args[1]);
    int nQueries = atoi(args[2]);
    int dimension = atoi(args[3]);
    double successProbability = atof(args[4]);
    int nNN = atoi(args[5]);

    std::string db = args[6];
    std::string db_test = args[7];
    float* dataset = new float[nPoints*dimension];
    float* queries = new float[nQueries*dimension];

    printf("Loading data: \n");
	clock_t begin = clock();
    load_data(db, dataset, nPoints, dimension);
    load_data(db_test, queries, nQueries, dimension);
    clock_t end = clock();
    printf("Exec time loadData (ms): %d\n", ( end - begin ));

    int distanceType = atof(args[8]);
    int nLevels =  atoi(args[9]);

    MetricEvaluator* metric;
    if (distanceType == 0)
        metric= new L1Distance();
    else
        metric= new L2Distance();

    metric->dim = dimension; 
   
	begin = clock();
    CUDAHashFile index(metric, dataset, nPoints, dimension, nLevels);
    end = clock();
    printf("Exec time init_hash_functions (ms): %d\n", ( end - begin ));

	begin = clock();
    index.BulkLoad();

    end = clock();
    printf("Exec time dev_lsh_bulkload (ms): %d\n", ( end - begin ));
     
    /*double meanQueryTime =  0.0f;
    for (int i = 0; i < nQueries; i++) {
        clock_t begin = clock();
        tResult result = index.NearestQuery( &queries[i * dimension], nNN);
        clock_t end = clock();
        double time     = (end - begin)/1000.0;
        printf("Total time for R-NN query at radius %0.6lf (radius no. %d):\t%0.6lf\n", (double)nNN, 0, time);
        meanQueryTime += time;

        if(result.GetNumOfEntries() == 0)
            printf("Query point %d: no NNs found.\n", i);
        else {
            printf("Query point %d: found %d NNs at distance %0.6lf (radius no. %d). NNs are:\n", i, result.GetNumOfEntries(), result.GetMaximumDistance(), 0);
            for(int j = 0; j < result.GetNumOfEntries(); j++){
                printf("%09d\tdist:%0.6lf\n", result[j].GetObject(), result[j].GetDistance());
            }
        }
    }
    if (nQueries > 0){
        meanQueryTime = meanQueryTime / nQueries;
        printf("Mean query time: %0.6lf\n", (double)meanQueryTime);
    }
    printf("Number of Distance Computations:%d\n", index.GetEvaluatorType()->GetDistanceCount());*/

    index.AllTopKMotifs(5, nNN);

    /*begin = clock();
    tResult result = index.TopKMotif(nNN);
    end = clock();
    double time     = (end - begin)/1000.0;

    printf("Total time for [Top-%d MOTIF] \t %0.6lf\n", nNN, time);

    if(result.GetNumOfEntries() == 0)
        printf("No Motif found.\n");
    else {
        printf("Motif found %d NNs at distance %0.6lf. NNs are:\n", result.GetNumOfEntries(), result.GetMaximumDistance());
        for(int j = 0; j < result.GetNumOfEntries(); j++){
            printf("%09d\tdist:%0.6lf\n", result[j].GetObject(), result[j].GetDistance());
        }
    }*/
    delete metric;
    delete []dataset;
    delete []queries;
    return 0;
}
