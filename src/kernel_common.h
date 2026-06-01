#pragma once

#include <string>

constexpr int kCliPPOk = 0;
constexpr int kCliPPError = 1;
constexpr int kCudaUnavailable = 2;
constexpr int kCudaFailedAfterWork = 3;

std::string lambda_to_string(double Lambda);

int validate_clipp_inputs(
    int No_mutation,
    const int* c_r,
    const int* c_n,
    const int* c_minor,
    const int* c_total,
    const double* Lambda_list,
    int Lambda_num,
    double ploidy,
    double purity,
    double alpha,
    double rho,
    double gamma,
    int Run_limit,
    double precision,
    int control_large,
    int least_mut,
    double post_th,
    double least_diff,
    const double* c_coef_1d,
    const double* c_wcut_1d,
    const char* preliminary);

int CliPPCPP(
    int No_mutation,
    int* r,
    int* n,
    int* minor,
    int* total,
    double ploidy,
    double* Lambda_list,
    int Lambda_num,
    double alpha,
    double rho,
    double gamma,
    int Run_limit,
    double precision,
    int control_large,
    int least_mut,
    double post_th,
    double least_diff,
    double* coef_1d,
    double* wcut_1d,
    double purity,
    char* preliminary);

#ifdef USE_CUDA
int CliPPCUDA(
    int No_mutation,
    int* r,
    int* n,
    int* minor,
    int* total,
    double ploidy,
    double* Lambda_list,
    int Lambda_num,
    double alpha,
    double rho,
    double gamma,
    int Run_limit,
    double precision,
    int control_large,
    int least_mut,
    double post_th,
    double least_diff,
    double* coef_1d,
    double* wcut_1d,
    double purity,
    char* preliminary);
#endif
