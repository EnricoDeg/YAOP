#ifdef __cplusplus
extern "C" {
#endif


// Constructor
void TKE_Init(int nproma, int nlevs, int nblocks);

// Calculation
void TKE_Calc(double * temperature);


#ifdef __cplusplus
}
#endif
