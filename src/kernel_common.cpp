#include <cmath>
#include <iostream>
#include <string>

#include "kernel_common.h"

std::string lambda_to_string(double Lambda)
{
    std::string s = std::to_string(Lambda);
    const size_t last_nonzero = s.find_last_not_of('0');

    if(last_nonzero == std::string::npos){
        return "0";
    }

    s.erase(last_nonzero + 1);

    if(!s.empty() && s.back() == '.'){
        s.pop_back();
    }

    if(s.empty() || s == "-0"){
        s = "0";
    }

    return s;
}

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
    const char* preliminary)
{
    if(No_mutation <= 0){
        std::cerr << "No_mutation must be positive." << std::endl;
        return kCliPPError;
    }

    if(Lambda_num <= 0){
        std::cerr << "Lambda_num must be positive." << std::endl;
        return kCliPPError;
    }

    if(!c_r || !c_n || !c_minor || !c_total || !Lambda_list || !c_coef_1d || !c_wcut_1d || !preliminary){
        std::cerr << "CliPP received a null input pointer." << std::endl;
        return kCliPPError;
    }

    if(!std::isfinite(ploidy) || !std::isfinite(purity)){
        std::cerr << "ploidy and purity must be finite." << std::endl;
        return kCliPPError;
    }

    if(alpha <= 0.0 || rho <= 0.0 || gamma <= 1.0 ||
       !std::isfinite(alpha) || !std::isfinite(rho) || !std::isfinite(gamma)){
        std::cerr << "Invalid optimization parameters: require finite alpha > 0, rho > 0, gamma > 1." << std::endl;
        return kCliPPError;
    }

    if(Run_limit <= 0){
        std::cerr << "Run_limit must be positive." << std::endl;
        return kCliPPError;
    }

    if(precision <= 0.0 || !std::isfinite(precision)){
        std::cerr << "precision must be positive and finite." << std::endl;
        return kCliPPError;
    }

    if(control_large <= 0){
        std::cerr << "control_large must be positive." << std::endl;
        return kCliPPError;
    }

    if(least_mut < 0){
        std::cerr << "least_mut must be nonnegative." << std::endl;
        return kCliPPError;
    }

    if(least_mut > No_mutation){
        std::cerr << "least_mut must be less than or equal to No_mutation." << std::endl;
        return kCliPPError;
    }

    if(post_th < 0.0 || !std::isfinite(post_th)){
        std::cerr << "post_th must be nonnegative and finite." << std::endl;
        return kCliPPError;
    }

    if(least_diff < 0.0 || !std::isfinite(least_diff)){
        std::cerr << "least_diff must be nonnegative and finite." << std::endl;
        return kCliPPError;
    }

    for(int l = 0; l < Lambda_num; ++l){
        if(!std::isfinite(Lambda_list[l]) || Lambda_list[l] < 0.0){
            std::cerr << "Invalid Lambda_list[" << l << "]." << std::endl;
            return kCliPPError;
        }
    }

    for(int i = 0; i < No_mutation; ++i){
        if(c_r[i] < 0){
            std::cerr << "Invalid r[" << i << "] < 0." << std::endl;
            return kCliPPError;
        }

        if(c_n[i] <= 0){
            std::cerr << "Invalid depth n[" << i << "] <= 0." << std::endl;
            return kCliPPError;
        }

        if(c_r[i] > c_n[i]){
            std::cerr << "Invalid r[" << i << "] > n[" << i << "]." << std::endl;
            return kCliPPError;
        }

        if(c_minor[i] <= 0){
            std::cerr << "Invalid minor[" << i << "] <= 0." << std::endl;
            return kCliPPError;
        }

        if(c_total[i] <= 0){
            std::cerr << "Invalid total[" << i << "] <= 0." << std::endl;
            return kCliPPError;
        }

        if(c_minor[i] > c_total[i]){
            std::cerr << "Invalid minor[" << i << "] > total[" << i << "]." << std::endl;
            return kCliPPError;
        }
    }

    return kCliPPOk;
}
