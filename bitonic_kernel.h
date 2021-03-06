
#include <cuda.h>

#include <stdio.h> 
#include <stdlib.h> 
#include <assert.h> 

#include "cutil.h"
#include "cutil_inline_runtime.h"

#define SHARED_SIZE_LIMIT 1024U
typedef unsigned int uint;

//Map to single instructions on G8x / G9x / G100
#define    UMUL(a, b) __umul24((a), (b))
#define UMAD(a, b, c) ( UMUL((a), (b)) + (c) )


__device__ inline void Comparator(
								  float& keyA,
								  uint& valA,
								  float& keyB,
								  uint& valB,
								  uint dir
								  ){
									  uint t;
									  if( (keyA > keyB) == dir ){
											t = keyA; keyA = keyB; keyB = t;
											t = valA; valA = valB; valB = t;
									  }
}


////////////////////////////////////////////////////////////////////////////////
// Monolithic bitonic sort kernel for short arrays fitting into shared memory
////////////////////////////////////////////////////////////////////////////////
__global__ void bitonicSortShared(
								  float *d_DstKey,
								  uint *d_DstVal,
								  float *d_SrcKey,
								  uint *d_SrcVal,
								  uint offset,
								  uint arrayLength,
								  uint dir
								  ){

								  	  d_SrcKey = &d_SrcKey[offset * arrayLength];
								  	  d_SrcVal = &d_SrcVal[offset * arrayLength];
									  
									  d_DstKey = &d_DstKey[offset * arrayLength];
								  	  d_DstVal = &d_DstVal[offset * arrayLength];
									  
									  //Shared memory storage for one or more short vectors
									  __shared__ float s_key[SHARED_SIZE_LIMIT];
									  __shared__ uint s_val[SHARED_SIZE_LIMIT];

									  //Offset to the beginning of subbatch and load data
									  d_SrcKey += blockIdx.x * SHARED_SIZE_LIMIT + threadIdx.x;
									  d_SrcVal += blockIdx.x * SHARED_SIZE_LIMIT + threadIdx.x;
									  d_DstKey += blockIdx.x * SHARED_SIZE_LIMIT + threadIdx.x;
									  d_DstVal += blockIdx.x * SHARED_SIZE_LIMIT + threadIdx.x;
									  s_key[threadIdx.x +                       0] = d_SrcKey[                      0];
									  s_val[threadIdx.x +                       0] = d_SrcVal[                      0];
									  s_key[threadIdx.x + (SHARED_SIZE_LIMIT / 2)] = d_SrcKey[(SHARED_SIZE_LIMIT / 2)];
									  s_val[threadIdx.x + (SHARED_SIZE_LIMIT / 2)] = d_SrcVal[(SHARED_SIZE_LIMIT / 2)];

									  for(uint size = 2; size < arrayLength; size <<= 1){
										  //Bitonic merge
										  uint ddd = dir ^ ( (threadIdx.x & (size / 2)) != 0 );
										  for(uint stride = size / 2; stride > 0; stride >>= 1){
											  __syncthreads();
											  uint pos = 2 * threadIdx.x - (threadIdx.x & (stride - 1));
											  Comparator(
												  s_key[pos +      0], s_val[pos +      0],
												  s_key[pos + stride], s_val[pos + stride],
												  ddd
												  );
										  }
									  }

									  //ddd == dir for the last bitonic merge step
									  {
										  for(uint stride = arrayLength / 2; stride > 0; stride >>= 1){
											  __syncthreads();
											  uint pos = 2 * threadIdx.x - (threadIdx.x & (stride - 1));
											  Comparator(
												  s_key[pos +      0], s_val[pos +      0],
												  s_key[pos + stride], s_val[pos + stride],
												  dir
												  );
										  }
									  }

									  __syncthreads();
									  d_DstKey[                      0] = s_key[threadIdx.x +                       0];
									  d_DstVal[                      0] = s_val[threadIdx.x +                       0];
									  d_DstKey[(SHARED_SIZE_LIMIT / 2)] = s_key[threadIdx.x + (SHARED_SIZE_LIMIT / 2)];
									  d_DstVal[(SHARED_SIZE_LIMIT / 2)] = s_val[threadIdx.x + (SHARED_SIZE_LIMIT / 2)];
}


__global__ void bitonicSortShared1(
								   float *d_DstKey,
								   uint *d_DstVal,
								   float *d_SrcKey,
								   uint *d_SrcVal,
								   uint offset,
								   uint arrayLength
								   ){

								  	  d_SrcKey = &d_SrcKey[offset * arrayLength];
								  	  d_SrcVal = &d_SrcVal[offset * arrayLength];
									  
									  d_DstKey = &d_DstKey[offset * arrayLength];
								  	  d_DstVal = &d_DstVal[offset * arrayLength];

									   //Shared memory storage for current subarray
									   __shared__ float s_key[SHARED_SIZE_LIMIT];
									   __shared__ uint s_val[SHARED_SIZE_LIMIT];

									   //Offset to the beginning of subarray and load data
									   d_SrcKey += blockIdx.x * SHARED_SIZE_LIMIT + threadIdx.x;
									   d_SrcVal += blockIdx.x * SHARED_SIZE_LIMIT + threadIdx.x;
									   d_DstKey += blockIdx.x * SHARED_SIZE_LIMIT + threadIdx.x;
									   d_DstVal += blockIdx.x * SHARED_SIZE_LIMIT + threadIdx.x;
									   s_key[threadIdx.x +                       0] = d_SrcKey[                      0];
									   s_val[threadIdx.x +                       0] = d_SrcVal[                      0];
									   s_key[threadIdx.x + (SHARED_SIZE_LIMIT / 2)] = d_SrcKey[(SHARED_SIZE_LIMIT / 2)];
									   s_val[threadIdx.x + (SHARED_SIZE_LIMIT / 2)] = d_SrcVal[(SHARED_SIZE_LIMIT / 2)];

									   for(uint size = 2; size < SHARED_SIZE_LIMIT; size <<= 1){
										   //Bitonic merge
										   uint ddd = (threadIdx.x & (size / 2)) != 0;
										   for(uint stride = size / 2; stride > 0; stride >>= 1){
											   __syncthreads();
											   uint pos = 2 * threadIdx.x - (threadIdx.x & (stride - 1));
											   Comparator(
												   s_key[pos +      0], s_val[pos +      0],
												   s_key[pos + stride], s_val[pos + stride],
												   ddd
												   );
										   }
									   }

									   //Odd / even arrays of SHARED_SIZE_LIMIT elements
									   //sorted in opposite directions
									   uint ddd = blockIdx.x & 1;
									   {
										   for(uint stride = SHARED_SIZE_LIMIT / 2; stride > 0; stride >>= 1){
											   __syncthreads();
											   uint pos = 2 * threadIdx.x - (threadIdx.x & (stride - 1));
											   Comparator(
												   s_key[pos +      0], s_val[pos +      0],
												   s_key[pos + stride], s_val[pos + stride],
												   ddd
												   );
										   }
									   }


									   __syncthreads();
									   d_DstKey[                      0] = s_key[threadIdx.x +                       0];
									   d_DstVal[                      0] = s_val[threadIdx.x +                       0];
									   d_DstKey[(SHARED_SIZE_LIMIT / 2)] = s_key[threadIdx.x + (SHARED_SIZE_LIMIT / 2)];
									   d_DstVal[(SHARED_SIZE_LIMIT / 2)] = s_val[threadIdx.x + (SHARED_SIZE_LIMIT / 2)];
}

//Bitonic merge iteration for stride >= SHARED_SIZE_LIMIT
__global__ void bitonicMergeGlobal(
								   float *d_DstKey,
								   uint *d_DstVal,
								   float *d_SrcKey,
								   uint *d_SrcVal,
								   uint offset,
								   uint arrayLength,
								   uint size,
								   uint stride,
								   uint dir
								   ){
								   	  d_SrcKey = &d_SrcKey[offset * arrayLength];
								  	  d_SrcVal = &d_SrcVal[offset * arrayLength];
									  
									  d_DstKey = &d_DstKey[offset * arrayLength];
								  	  d_DstVal = &d_DstVal[offset * arrayLength];

									   uint global_comparatorI = blockIdx.x * blockDim.x + threadIdx.x;
									   uint        comparatorI = global_comparatorI & (arrayLength / 2 - 1);

									   //Bitonic merge
									   uint ddd = dir ^ ( (comparatorI & (size / 2)) != 0 );
									   uint pos = 2 * global_comparatorI - (global_comparatorI & (stride - 1));

									   float keyA = d_SrcKey[pos +      0];
									   uint valA = d_SrcVal[pos +      0];
									   float keyB = d_SrcKey[pos + stride];
									   uint valB = d_SrcVal[pos + stride];

									   Comparator(
										   keyA, valA,
										   keyB, valB,
										   ddd
										   );

									   d_DstKey[pos +      0] = keyA;
									   d_DstVal[pos +      0] = valA;
									   d_DstKey[pos + stride] = keyB;
									   d_DstVal[pos + stride] = valB;
}

//Combined bitonic merge steps for
//size > SHARED_SIZE_LIMIT and stride = [1 .. SHARED_SIZE_LIMIT / 2]
__global__ void bitonicMergeShared(
								   float *d_DstKey,
								   uint *d_DstVal,
								   float *d_SrcKey,
								   uint *d_SrcVal,
								   uint offset,
								   uint arrayLength,
								   uint size,
								   uint dir
								   ){
								   	  d_SrcKey = &d_SrcKey[offset * arrayLength];
								  	  d_SrcVal = &d_SrcVal[offset * arrayLength];
									  
									  d_DstKey = &d_DstKey[offset * arrayLength];
								  	  d_DstVal = &d_DstVal[offset * arrayLength];
								  	  
									   //Shared memory storage for current subarray
									   __shared__ float s_key[SHARED_SIZE_LIMIT];
									   __shared__ uint s_val[SHARED_SIZE_LIMIT];

									   d_SrcKey += blockIdx.x * SHARED_SIZE_LIMIT + threadIdx.x;
									   d_SrcVal += blockIdx.x * SHARED_SIZE_LIMIT + threadIdx.x;
									   d_DstKey += blockIdx.x * SHARED_SIZE_LIMIT + threadIdx.x;
									   d_DstVal += blockIdx.x * SHARED_SIZE_LIMIT + threadIdx.x;
									   s_key[threadIdx.x +                       0] = d_SrcKey[                      0];
									   s_val[threadIdx.x +                       0] = d_SrcVal[                      0];
									   s_key[threadIdx.x + (SHARED_SIZE_LIMIT / 2)] = d_SrcKey[(SHARED_SIZE_LIMIT / 2)];
									   s_val[threadIdx.x + (SHARED_SIZE_LIMIT / 2)] = d_SrcVal[(SHARED_SIZE_LIMIT / 2)];

									   //Bitonic merge
									   uint comparatorI = UMAD(blockIdx.x, blockDim.x, threadIdx.x) & ((arrayLength / 2) - 1);
									   uint ddd = dir ^ ( (comparatorI & (size / 2)) != 0 );
									   for(uint stride = SHARED_SIZE_LIMIT / 2; stride > 0; stride >>= 1){
										   __syncthreads();
										   uint pos = 2 * threadIdx.x - (threadIdx.x & (stride - 1));
										   Comparator(
											   s_key[pos +      0], s_val[pos +      0],
											   s_key[pos + stride], s_val[pos + stride],
											   ddd
											   );
									   }

									   __syncthreads();
									   d_DstKey[                      0] = s_key[threadIdx.x +                       0];
									   d_DstVal[                      0] = s_val[threadIdx.x +                       0];
									   d_DstKey[(SHARED_SIZE_LIMIT / 2)] = s_key[threadIdx.x + (SHARED_SIZE_LIMIT / 2)];
									   d_DstVal[(SHARED_SIZE_LIMIT / 2)] = s_val[threadIdx.x + (SHARED_SIZE_LIMIT / 2)];
}



uint factorRadix2(uint *log2L, uint L){
	if(!L){
		*log2L = 0;
		return 0;
	}else{
		for(*log2L = 0; (L & 1) == 0; L >>= 1, *log2L++);
		return L;
	}
}
 uint bitonicSort(
			float *d_DstKey,
			uint *d_DstVal,
			float *d_SrcKey,
			uint *d_SrcVal,
			uint batchSize,
			uint offset,
			uint arrayLength,
			uint dir
			){
				//Nothing to sort
				if(arrayLength < 2)
					return 0;

				//Only power-of-two array lengths are supported by this implementation
				uint log2L;
				uint factorizationRemainder = factorRadix2(&log2L, arrayLength);
				assert( factorizationRemainder == 1 );

				dir = (dir != 0);

				uint  blockCount = batchSize * arrayLength / SHARED_SIZE_LIMIT;
				uint threadCount = SHARED_SIZE_LIMIT / 2;

				if(arrayLength <= SHARED_SIZE_LIMIT){
					assert( (batchSize * arrayLength) % SHARED_SIZE_LIMIT == 0 );
					bitonicSortShared<<<blockCount, threadCount>>>(d_DstKey, d_DstVal, d_SrcKey, d_SrcVal, offset, arrayLength, dir);
				}else{
					bitonicSortShared1<<<blockCount, threadCount>>>(d_DstKey, d_DstVal, d_SrcKey, d_SrcVal, offset, arrayLength);

					for(uint size = 2 * SHARED_SIZE_LIMIT; size <= arrayLength; size <<= 1)
						for(unsigned stride = size / 2; stride > 0; stride >>= 1)
							if(stride >= SHARED_SIZE_LIMIT){
								bitonicMergeGlobal<<<(batchSize * arrayLength) / 512, 256>>>(d_DstKey, d_DstVal, d_DstKey, d_DstVal, offset, arrayLength, size, stride, dir);
							}else{
								bitonicMergeShared<<<blockCount, threadCount>>>(d_DstKey, d_DstVal, d_DstKey, d_DstVal, offset, arrayLength, size, dir);
								break;
							}
				}
				return threadCount;   
}


#define shrLog printf

 ////////////////////////////////////////////////////////////////////////////////
 // Validate sorted keys array (check for integrity and proper order)
 ////////////////////////////////////////////////////////////////////////////////
 uint validateSortedKeys(
	 uint *resKey,
	 uint *srcKey,
	 uint batchSize,
	 uint arrayLength,
	 uint numValues,
	 uint dir
	 ){
		 uint *srcHist;
		 uint *resHist;
		 if(arrayLength < 2){
			 shrLog("validateSortedKeys(): arrayLength too short, exiting...\n");
			 return 1;
		 }

		 shrLog("...inspecting keys array: ");

		 srcHist = (uint *)malloc(numValues * sizeof(uint));
		 resHist = (uint *)malloc(numValues * sizeof(uint));

		 int flag = 1;
		 for(uint j = 0; j < batchSize; j++, srcKey += arrayLength, resKey += arrayLength){
			 //Build histograms for keys arrays
			 memset(srcHist, 0, numValues * sizeof(uint));
			 memset(resHist, 0, numValues * sizeof(uint));
			 for(uint i = 0; i < arrayLength; i++){
				 if( srcKey[i] < numValues && resKey[i] < numValues ){
					 srcHist[srcKey[i]]++;
					 resHist[resKey[i]]++;
				 }else{
					 flag = 0;
					 break;
				 }
			 }

			 if(!flag){
				 shrLog("***Set %u source/result key arrays are not limited properly***\n", j);
				 goto brk;
			 }

			 //Compare the histograms
			 for(uint i = 0; i < numValues; i++)
				 if(srcHist[i] != resHist[i]){
					 flag = 0;
					 break;
				 }

				 if(!flag){
					 shrLog("***Set %u source/result keys histograms do not match***\n", j);
					 goto brk;
				 }

				 if(dir){
					 //Ascending order
					 for(uint i = 0; i < arrayLength - 1; i++)
						 if(resKey[i + 1] < resKey[i]){
							 flag = 0;
							 break;
						 }
				 }else{
					 //Descending order
					 for(uint i = 0; i < arrayLength - 1; i++)
						 if(resKey[i + 1] > resKey[i]){
							 flag = 0;
							 break;
						 }
				 }

				 if(!flag){
					 shrLog("***Set %u result key array is not ordered properly***\n", j);
					 goto brk;
				 }
		 }

brk:
		 free(resHist);
		 free(srcHist);

		 if(flag) shrLog("OK\n");
		 return flag;
 }



 int validateValues(
	 uint *resKey,
	 uint *resVal,
	 uint *srcKey,
	 uint batchSize,
	 uint arrayLength
	 ){
		 int correctFlag = 1, stableFlag = 1;

		 shrLog("...inspecting keys and values array: ");
		 for(uint i = 0; i < batchSize; i++, resKey += arrayLength, resVal += arrayLength){
			 for(uint j = 0; j < arrayLength; j++){
				 if(resKey[j] != srcKey[resVal[j]])
					 correctFlag = 0;

				 if( (j < arrayLength - 1) && (resKey[j] == resKey[j + 1]) && (resVal[j] > resVal[j + 1]) )
					 stableFlag = 0;
			 }
		 }

		 shrLog(correctFlag ? "OK\n" : "***corrupted!!!***\n");
		 shrLog(stableFlag ? "...stability property: stable!\n" : "...stability property: NOT stable\n");

		 return correctFlag;
 }
