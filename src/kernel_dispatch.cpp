#include <cstdlib>
#include <iostream>

#include "kernel_common.h"

extern "C" {
    int CliPPStatus(int No_mutation, int* r, int *n, int* minor, int* total, double ploidy,
	      double* Lambda_list, int Lambda_num, double alpha, double rho, double gamma, int Run_limit, double precision,
	      int control_large, int least_mut, double post_th, double least_diff,
	      double* coef_1d, double* wcut_1d, double purity, char* preliminary){

#ifdef USE_CUDA
        if(std::getenv("CLIPP_FORCE_CPU") == nullptr){
            int cuda_status = CliPPCUDA(No_mutation, r, n, minor, total, ploidy, Lambda_list, Lambda_num, alpha, rho, gamma, Run_limit, precision, control_large, least_mut, post_th, least_diff, coef_1d, wcut_1d, purity, preliminary);
            if(cuda_status == kCliPPOk) return kCliPPOk;
            if(cuda_status == kCudaUnavailable){
                std::cerr << "CUDA unavailable. Falling back to CPU CliPP implementation." << std::endl;
            }else{
                std::cerr << "CUDA CliPP implementation failed; not falling back automatically." << std::endl;
                return cuda_status;
            }
        }
#endif
	return CliPPCPP(No_mutation, r, n, minor, total, ploidy, Lambda_list, Lambda_num, alpha, rho, gamma, Run_limit, precision, control_large, least_mut, post_th, least_diff, coef_1d, wcut_1d, purity, preliminary);
    }

    void CliPP(int No_mutation, int* r, int *n, int* minor, int* total, double ploidy,
	      double* Lambda_list, int Lambda_num, double alpha, double rho, double gamma, int Run_limit, double precision,
	      int control_large, int least_mut, double post_th, double least_diff,
	      double* coef_1d, double* wcut_1d, double purity, char* preliminary){

        const int rc = CliPPStatus(No_mutation, r, n, minor, total, ploidy, Lambda_list, Lambda_num, alpha, rho, gamma, Run_limit, precision, control_large, least_mut, post_th, least_diff, coef_1d, wcut_1d, purity, preliminary);
        if(rc != kCliPPOk){
            std::cerr << "CliPP failed with status " << rc << std::endl;
        }
    }
}
