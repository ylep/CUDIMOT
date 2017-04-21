#define SIZE_PART 12800 	// Number of voxels to compute in each subpart of the volume. 
				// It uses about 500MB allocated in global memory
				// Could increase it, but almost the same performance
				// Enough to exploit parallelism and low memory requirements for old/laptop GPUs

#define VOXELS_BLOCK 8		// Number of voxels per Block. 

// 12800/8  = 1600 blocks




