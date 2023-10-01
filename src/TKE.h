#ifdef __cplusplus
extern "C" {
#endif


// Constructor
void TKE_Init(int nproma, int nlevs, int nblocks);

// Destructor
void TKE_Finalize();

// Calculation
void TKE_Calc(double * temperature);


#ifdef __cplusplus
}
#endif
