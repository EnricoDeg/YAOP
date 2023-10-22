#ifdef __cplusplus
extern "C" {
#endif


// Constructor
void TKE_Init(int nproma, int nlevs, int nblocks,
              int block_size, int start_index, int end_index);

// Destructor
void TKE_Finalize();

// Calculation
void TKE_Calc(int start_block, int end_block, double * tke, int *dolic_c);


#ifdef __cplusplus
}
#endif
