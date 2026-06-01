#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <fstream>
#include <algorithm>
#include <cstdlib>
#include <cstdint>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <atomic>
#include <unordered_map>
#include <mutex>
//#include <chrono>
#include <Eigen/Dense>
#include "kernel_common.h"
#ifdef _OPENMP 
#include <omp.h>
#endif
#ifdef USE_CUDA
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <nvrtc.h>
#endif

using Eigen::MatrixXd;
using Eigen::MatrixXi;

using namespace std;

class CudaPackedPairDiff {
public:
    CudaPackedPairDiff(const std::vector<double>& eta, const std::vector<long long>& pair_start)
        : eta_(eta), pair_start_(pair_start)
    {
    }

    double value(int a, int b) const
    {
        if(a == b) return 0.0;
        const long long p = pair_index(a, b);
        auto it = overrides_.find(p);
        if(it != overrides_.end()) return it->second;
        return eta_[p];
    }

    void set(int a, int b, double value)
    {
        if(a >= b) return;
        overrides_[pair_start_[a] + b - a - 1] = value;
    }

private:
    long long pair_index(int a, int b) const
    {
        if(a > b) std::swap(a, b);
        return pair_start_[a] + b - a - 1;
    }

    const std::vector<double>& eta_;
    const std::vector<long long>& pair_start_;
    mutable std::unordered_map<long long, double> overrides_;
};

#ifdef USE_CUDA
namespace {

constexpr int kCudaBlockSize = 256;
constexpr size_t kCudaMemoryOverheadBytes = 256ULL * 1024ULL * 1024ULL;

std::runtime_error make_cuda_error(const char* kind, const char* call, const char* file, int line, const std::string& message)
{
    std::ostringstream oss;
    oss << kind << " error at " << file << ":" << line << " while calling " << call << ": " << message;
    return std::runtime_error(oss.str());
}

void check_cuda_runtime(cudaError_t err, const char* call, const char* file, int line)
{
    if(err != cudaSuccess){
        throw make_cuda_error("CUDA runtime", call, file, line, cudaGetErrorString(err));
    }
}

void check_cuda_driver(CUresult err, const char* call, const char* file, int line)
{
    if(err != CUDA_SUCCESS){
        const char* name = nullptr;
        const char* message = nullptr;
        cuGetErrorName(err, &name);
        cuGetErrorString(err, &message);
        std::ostringstream oss;
        if(name) oss << name;
        if(message) oss << " (" << message << ")";
        throw make_cuda_error("CUDA driver", call, file, line, oss.str());
    }
}

void check_nvrtc(nvrtcResult err, const char* call, const char* file, int line)
{
    if(err != NVRTC_SUCCESS){
        throw make_cuda_error("NVRTC", call, file, line, nvrtcGetErrorString(err));
    }
}

#define CLIPP_CUDA_CHECK(call) check_cuda_runtime((call), #call, __FILE__, __LINE__)
#define CLIPP_CUDA_DRIVER_CHECK(call) check_cuda_driver((call), #call, __FILE__, __LINE__)
#define CLIPP_NVRTC_CHECK(call) check_nvrtc((call), #call, __FILE__, __LINE__)

#ifndef CLIPP_VERBOSE
#define CLIPP_VERBOSE 0
#endif

size_t estimate_cuda_lambda_bytes(int No_mutation)
{
    const long long pair_count = static_cast<long long>(No_mutation) * (No_mutation - 1) / 2;
    const size_t edge_bytes_per_pair = 2ULL * sizeof(double) + sizeof(unsigned char);
    if(pair_count > 0 && static_cast<unsigned long long>(pair_count) > std::numeric_limits<size_t>::max() / edge_bytes_per_pair){
        return std::numeric_limits<size_t>::max();
    }

    const size_t edge_bytes =
        static_cast<size_t>(pair_count) * edge_bytes_per_pair;
    const size_t node_double_count =
        2ULL   // w_new, w_old
      + 5ULL   // A, B, linear, Minv, dot_terms
      + 1ULL;  // row_residual_max
    const size_t node_bytes = static_cast<size_t>(No_mutation) * sizeof(double) * node_double_count
        + static_cast<size_t>(No_mutation) * sizeof(int);
    const size_t reduction_bytes = static_cast<size_t>(std::max<long long>(1, (No_mutation + 2LL * kCudaBlockSize - 1) / (2LL * kCudaBlockSize))) * sizeof(double) * 2ULL;

    if(edge_bytes > std::numeric_limits<size_t>::max() - node_bytes ||
       edge_bytes + node_bytes > std::numeric_limits<size_t>::max() - reduction_bytes ||
       edge_bytes + node_bytes + reduction_bytes > std::numeric_limits<size_t>::max() - kCudaMemoryOverheadBytes){
        return std::numeric_limits<size_t>::max();
    }

    return edge_bytes + node_bytes + reduction_bytes + kCudaMemoryOverheadBytes;
}

class CudaStream {
public:
    CudaStream()
    {
        CLIPP_CUDA_CHECK(cudaStreamCreateWithFlags(&stream_, cudaStreamNonBlocking));
    }

    ~CudaStream()
    {
        if(stream_){
            cudaStreamDestroy(stream_);
        }
    }

    CudaStream(const CudaStream&) = delete;
    CudaStream& operator=(const CudaStream&) = delete;

    cudaStream_t get() const { return stream_; }

private:
    cudaStream_t stream_ = nullptr;
};

template<typename T>
class GpuBuffer {
public:
    GpuBuffer() = default;
    explicit GpuBuffer(size_t count) { reset(count); }
    ~GpuBuffer() { release(); }

    GpuBuffer(const GpuBuffer&) = delete;
    GpuBuffer& operator=(const GpuBuffer&) = delete;

    GpuBuffer(GpuBuffer&& other) noexcept : ptr_(other.ptr_), count_(other.count_)
    {
        other.ptr_ = nullptr;
        other.count_ = 0;
    }

    GpuBuffer& operator=(GpuBuffer&& other) noexcept
    {
        if(this != &other){
            release();
            ptr_ = other.ptr_;
            count_ = other.count_;
            other.ptr_ = nullptr;
            other.count_ = 0;
        }
        return *this;
    }

    void reset(size_t count)
    {
        release();
        if(count > std::numeric_limits<size_t>::max() / sizeof(T)){
            throw std::overflow_error("GpuBuffer allocation size overflow");
        }
        count_ = count;
        if(count_ > 0){
            CLIPP_CUDA_CHECK(cudaMalloc(reinterpret_cast<void**>(&ptr_), count_ * sizeof(T)));
        }
    }

    void release()
    {
        if(ptr_){
            cudaFree(ptr_);
            ptr_ = nullptr;
        }
        count_ = 0;
    }

    T* get() const { return ptr_; }
    size_t count() const { return count_; }

    void swap(GpuBuffer& other)
    {
        std::swap(ptr_, other.ptr_);
        std::swap(count_, other.count_);
    }

private:
    T* ptr_ = nullptr;
    size_t count_ = 0;
};

template<typename T>
void copy_to_device(GpuBuffer<T>& dst, const std::vector<T>& src, cudaStream_t stream = nullptr)
{
    dst.reset(src.size());
    if(!src.empty()){
        CLIPP_CUDA_CHECK(cudaMemcpyAsync(dst.get(), src.data(), src.size() * sizeof(T), cudaMemcpyHostToDevice, stream));
    }
}

template<typename T>
void copy_to_host(std::vector<T>& dst, const GpuBuffer<T>& src, cudaStream_t stream = nullptr)
{
    dst.resize(src.count());
    if(src.count() > 0){
        CLIPP_CUDA_CHECK(cudaMemcpyAsync(dst.data(), src.get(), src.count() * sizeof(T), cudaMemcpyDeviceToHost, stream));
        CLIPP_CUDA_CHECK(cudaStreamSynchronize(stream));
    }
}

struct CudaStaticInputs {
    GpuBuffer<double> d_n;
    GpuBuffer<double> d_minor;
    GpuBuffer<double> d_total;
    GpuBuffer<double> d_theta_hat;
    GpuBuffer<double> d_phi_hat;
    GpuBuffer<double> d_coef;
    GpuBuffer<double> d_wcut;

    void reset(
        const std::vector<double>& n,
        const std::vector<double>& minor,
        const std::vector<double>& total,
        const std::vector<double>& theta_hat,
        const std::vector<double>& phi_hat,
        const std::vector<double>& coef,
        const std::vector<double>& wcut,
        cudaStream_t stream)
    {
        copy_to_device(d_n, n, stream);
        copy_to_device(d_minor, minor, stream);
        copy_to_device(d_total, total, stream);
        copy_to_device(d_theta_hat, theta_hat, stream);
        copy_to_device(d_phi_hat, phi_hat, stream);
        copy_to_device(d_coef, coef, stream);
        copy_to_device(d_wcut, wcut, stream);
        CLIPP_CUDA_CHECK(cudaStreamSynchronize(stream));
    }
};

const char* cuda_solver_source()
{
    // The included file is a C++ raw string literal consumed by NVRTC.
    return
#include "kernel_cuda_kernels.inc"
        ;
}

class CudaProgram {
public:
    CudaProgram()
    {
        cudaDeviceProp prop{};
        CLIPP_CUDA_CHECK(cudaGetDeviceProperties(&prop, 0));
        compile(prop.major, prop.minor);
    }

    ~CudaProgram()
    {
        if(module_) cuModuleUnload(module_);
    }

    CudaProgram(const CudaProgram&) = delete;
    CudaProgram& operator=(const CudaProgram&) = delete;

    void launch(const char* kernel_name, dim3 grid, dim3 block, void** args, size_t shared_bytes = 0, cudaStream_t stream = nullptr)
    {
        if(grid.x == 0 || grid.y == 0 || grid.z == 0 || block.x == 0 || block.y == 0 || block.z == 0){
            throw std::runtime_error(std::string("Invalid CUDA launch dimensions for ") + kernel_name);
        }

        CUfunction fn = get_function(kernel_name);
        CLIPP_CUDA_DRIVER_CHECK(cuLaunchKernel(
            fn,
            grid.x, grid.y, grid.z,
            block.x, block.y, block.z,
            static_cast<unsigned int>(shared_bytes),
            reinterpret_cast<CUstream>(stream),
            args,
            nullptr));
#ifdef CLIPP_CUDA_DEBUG_SYNC
        if(stream){
            CLIPP_CUDA_CHECK(cudaStreamSynchronize(stream));
        }else{
            CLIPP_CUDA_DRIVER_CHECK(cuCtxSynchronize());
        }
#endif
    }

private:
    CUfunction get_function(const char* kernel_name)
    {
        std::lock_guard<std::mutex> lock(function_mutex_);
        auto it = functions_.find(kernel_name);
        if(it != functions_.end()) return it->second;

        CUfunction fn;
        CLIPP_CUDA_DRIVER_CHECK(cuModuleGetFunction(&fn, module_, kernel_name));
        functions_.emplace(kernel_name, fn);
        return fn;
    }

    void compile(int major, int minor)
    {
        nvrtcProgram program;
        CLIPP_NVRTC_CHECK(nvrtcCreateProgram(
            &program,
            cuda_solver_source(),
            "clipp_cuda_solver.cu",
            0,
            nullptr,
            nullptr));

        std::string arch = "--gpu-architecture=compute_" + std::to_string(major) + std::to_string(minor);
        const char* options[] = {
            arch.c_str(),
            "--std=c++17"
        };

        nvrtcResult compile_result = nvrtcCompileProgram(program, 2, options);

        size_t log_size = 0;
        CLIPP_NVRTC_CHECK(nvrtcGetProgramLogSize(program, &log_size));
        std::string log;
        if(log_size > 1){
            log.resize(log_size);
            CLIPP_NVRTC_CHECK(nvrtcGetProgramLog(program, &log[0]));
        }

        if(compile_result != NVRTC_SUCCESS){
            nvrtcDestroyProgram(&program);
            std::ostringstream oss;
            oss << "failed to compile CUDA kernels";
            if(!log.empty()) oss << ":\n" << log;
            throw std::runtime_error(oss.str());
        }

        if(!log.empty()){
            std::string trimmed = log;
            trimmed.erase(std::remove(trimmed.begin(), trimmed.end(), '\0'), trimmed.end());
            if(!trimmed.empty()){
                std::cerr << "NVRTC log:\n" << trimmed << std::endl;
            }
        }

        size_t ptx_size = 0;
        CLIPP_NVRTC_CHECK(nvrtcGetPTXSize(program, &ptx_size));
        std::string ptx(ptx_size, '\0');
        CLIPP_NVRTC_CHECK(nvrtcGetPTX(program, &ptx[0]));
        CLIPP_NVRTC_CHECK(nvrtcDestroyProgram(&program));

        CLIPP_CUDA_DRIVER_CHECK(cuModuleLoadDataEx(&module_, ptx.data(), 0, nullptr, nullptr));
    }

    CUmodule module_ = nullptr;
    std::unordered_map<std::string, CUfunction> functions_;
    std::mutex function_mutex_;
};

CudaProgram& cuda_program()
{
    static CudaProgram program;
    return program;
}

unsigned int blocks_for(long long count, int block_size = kCudaBlockSize)
{
    if(count <= 0) return 0;
    return static_cast<unsigned int>((count + block_size - 1) / block_size);
}

unsigned int reduction_blocks_for(long long count, int block_size = kCudaBlockSize)
{
    if(count <= 0) return 0;
    return static_cast<unsigned int>((count + (2LL * block_size) - 1) / (2LL * block_size));
}

struct CliPPCudaWorkspace {
    GpuBuffer<double> d_w_new;
    GpuBuffer<double> d_w_old;
    GpuBuffer<double> d_A;
    GpuBuffer<double> d_B;
    GpuBuffer<double> d_linear;
    GpuBuffer<double> d_Minv;
    GpuBuffer<double> d_dot_terms;
    GpuBuffer<double> d_eta;
    GpuBuffer<double> d_tau;
    GpuBuffer<double> d_row_residual_max;
    GpuBuffer<double> d_reduce_tmp_a;
    GpuBuffer<double> d_reduce_tmp_b;
    GpuBuffer<int> d_problematic_flags;
    GpuBuffer<unsigned char> d_fused;

    void reset(int No_mutation, long long pair_count)
    {
        d_w_new.reset(No_mutation);
        d_w_old.reset(No_mutation);
        d_A.reset(No_mutation);
        d_B.reset(No_mutation);
        d_linear.reset(No_mutation);
        d_Minv.reset(No_mutation);
        d_dot_terms.reset(No_mutation);
        d_eta.reset(static_cast<size_t>(pair_count));
        d_tau.reset(static_cast<size_t>(pair_count));
        d_row_residual_max.reset(No_mutation);
        d_problematic_flags.reset(No_mutation);
        d_fused.reset(static_cast<size_t>(pair_count));

        const long long max_reduce_input = No_mutation;
        const size_t reduce_temp_count = static_cast<size_t>(std::max<long long>(1, reduction_blocks_for(max_reduce_input)));
        d_reduce_tmp_a.reset(reduce_temp_count);
        d_reduce_tmp_b.reset(reduce_temp_count);
    }
};

double reduce_device_values(CudaProgram& program, const char* kernel_name, double* input, long long count, GpuBuffer<double>& tmp_a, GpuBuffer<double>& tmp_b, cudaStream_t stream)
{
    if(count <= 0) return 0.0;
    if(count == 1){
        double result = 0.0;
        CLIPP_CUDA_CHECK(cudaMemcpyAsync(&result, input, sizeof(double), cudaMemcpyDeviceToHost, stream));
        CLIPP_CUDA_CHECK(cudaStreamSynchronize(stream));
        return result;
    }

    long long current_count = count;
    double* current = input;
    bool use_a = true;

    while(current_count > 1){
        unsigned int grid = reduction_blocks_for(current_count);
        double* output = use_a ? tmp_a.get() : tmp_b.get();
        void* args[] = { &current, &output, &current_count };
        program.launch(kernel_name, dim3(grid), dim3(kCudaBlockSize), args, kCudaBlockSize * sizeof(double), stream);
        current = output;
        current_count = grid;
        use_a = !use_a;
    }

    double result = 0.0;
    CLIPP_CUDA_CHECK(cudaMemcpyAsync(&result, current, sizeof(double), cudaMemcpyDeviceToHost, stream));
    CLIPP_CUDA_CHECK(cudaStreamSynchronize(stream));
    return result;
}

double reduce_device_sum(CudaProgram& program, double* input, long long count, GpuBuffer<double>& tmp_a, GpuBuffer<double>& tmp_b, cudaStream_t stream)
{
    return reduce_device_values(program, "reduce_sum_kernel", input, count, tmp_a, tmp_b, stream);
}

double reduce_device_max(CudaProgram& program, double* input, long long count, GpuBuffer<double>& tmp_a, GpuBuffer<double>& tmp_b, cudaStream_t stream)
{
    if(count <= 0) return 0.0;
    return reduce_device_values(program, "reduce_max_kernel", input, count, tmp_a, tmp_b, stream);
}

constexpr int kCliPPNeedFullEta = 4;

void build_cpu_order_class_labels_from_fused(
    int No_mutation,
    const std::vector<unsigned char>& fused,
    const std::vector<long long>& pair_start,
    std::vector<int>& class_label,
    std::vector<int>& group_size)
{
    class_label.assign(No_mutation, -1);
    class_label[0] = 0;

    group_size.clear();
    group_size.push_back(1);

    int labl = 1;
    for(int i = 0; i < No_mutation; i++){
        for(int j = 0; j < i; j++){
            const long long pair_index = pair_start[j] + i - j - 1;
            if(fused[pair_index] != 0){
                class_label[i] = class_label[j];
                group_size[class_label[j]] = group_size[class_label[j]] + 1;
                break;
            }
        }
        if(class_label[i] == -1){
            class_label[i] = labl;
            labl += 1;
            group_size.push_back(1);
        }
    }
}

int cpu_order_minimum_group_size(const std::vector<int>& group_size, int No_mutation)
{
    int temp_size = No_mutation;
    for(unsigned int ii = 0; ii < group_size.size(); ii++){
        if(group_size[ii] > 0 && temp_size > group_size[ii]) temp_size = group_size[ii];
    }
    return temp_size;
}

int finalize_cuda_postprocess_from_labels(
    int No_mutation,
    double Lambda,
    int k,
    double least_diff,
    const std::vector<double>& n,
    const std::vector<double>& phi_hat,
    std::vector<int>& class_label,
    const std::vector<int>& problematic_snvs,
    std::string& preliminary_folder)
{
    int i, j, count;
    double temp, temp1;

    std::vector<int> labels = class_label;
    std::sort(labels.begin(), labels.end());
    vector<int>::iterator ip;
    ip = std::unique(labels.begin(), labels.end());
    labels.resize(std::distance(labels.begin(), ip));

    std::vector<double> phi_out(labels.size());
    std::fill(phi_out.begin(), phi_out.end(), 0);

    for(unsigned int ii = 0; ii < labels.size(); ii++){
        std::vector<int> ind;
        ind.reserve(No_mutation);

        for(i = 0; i < No_mutation; i++) {
            if(class_label[i] == labels[ii]) ind.push_back(i);
        }

        temp = 0;
        temp1 = 0;
        for(i = 0; i < int(ind.size()); i++) {
            class_label[ind[i]] = ii;
            temp += phi_hat[ind[i]] * n[ind[i]];
            temp1 += n[ind[i]];
        }

        phi_out[ii] = temp / temp1;
    }

    if(labels.size() > 1){
        std::vector<double> sort_phi = phi_out;
        std::sort(sort_phi.begin(), sort_phi.end());
        std::vector<double> phi_diff;

        for(i = 1; i < int(sort_phi.size()); i++) {
            phi_diff.push_back(sort_phi[i] - sort_phi[i -1]);
        }

        double min_val = 10000000000;
        int min_ind = 0;

        for(i = 0; i < int(phi_diff.size()); i++){
            if (min_val > phi_diff[i]){
                min_val = phi_diff[i];
                min_ind = i;
            }
        }

        std::vector<double> min_val_vector;
        min_val_vector.push_back(min_val);

        while(min_val < least_diff){
            std::vector<int> combine_ind, combine_to_ind;
            for(i = 0; i < int(phi_out.size()); i++) {
                if(phi_out[i] == sort_phi[min_ind]) combine_ind.push_back(i);
                if(phi_out[i] == sort_phi[min_ind + 1]) combine_to_ind.push_back(i);
            }

            for(i = 0; i < int(class_label.size()); i++) {
                if(class_label[i] == combine_ind[0]) class_label[i] = combine_to_ind[0];
            }

            labels.clear();
            labels = class_label;
            std::sort(labels.begin(), labels.end());
            ip = std::unique(labels.begin(), labels.end());
            labels.resize(std::distance(labels.begin(), ip));

            phi_out.clear();
            for(i = 0; i < int(labels.size()); i++) {
                phi_out.push_back(0.0);
            }

            for(unsigned int ii = 0; ii < labels.size(); ii++){
                std::vector<int> ind;
                ind.reserve(No_mutation);

                for(i = 0; i < No_mutation; i++) {
                    if(class_label[i] == labels[ii]) ind.push_back(i);
                }
                temp = 0;
                temp1 = 0;

                for(i = 0; i < int(ind.size()); i++) {
                    class_label[ind[i]] = ii;
                    temp += phi_hat[ind[i]] * n[ind[i]];
                    temp1 += n[ind[i]];
                }
                phi_out[ii] = temp / temp1;
            }

            if(labels.size() == 1) break;
            else{
                sort_phi = phi_out;
                std::sort(sort_phi.begin(), sort_phi.end());
                phi_diff.clear();
                for(i = 1; i < int(sort_phi.size()); i++) {
                    phi_diff.push_back(sort_phi[i] - sort_phi[i -1]);
                }

                min_val = 10000000000;
                for(i = 0; i < int(phi_diff.size()); i++){
                    if (min_val > phi_diff[i]){
                        min_val = phi_diff[i];
                        min_ind = i;
                    }
                }

                min_val_vector.push_back(min_val);

                if(min_val_vector.size() > 4){
                    int unique_min_val = 1;
                    for(double value:min_val_vector){
                        if(value != min_val) unique_min_val += 1;
                    }
                    if(unique_min_val == 1){
                        std::cout << "Lambda: " << Lambda << "\titeration: " << k << "\tfailed" << std::endl;
                        return -1;
                    }
                }
            }
        }
    }

    ofstream phi_file, label_file, summary_file;

    std::string lambda_str = lambda_to_string(Lambda);

    std::string phi_file_path = preliminary_folder + "/lam" + lambda_str + "_phi.txt";
    std::string label_file_path = preliminary_folder + "/lam" + lambda_str + "_label.txt";
    std::string summary_file_path = preliminary_folder + "/lam" + lambda_str + "_summary_table.txt";

    phi_file.open(phi_file_path);
    label_file.open(label_file_path);
    summary_file.open(summary_file_path);

    if(!phi_file.is_open()){
        std::cerr << "Cannot open file " << phi_file_path  << std::endl;
        return 1;
    }

    if(!label_file.is_open()){
        std::cerr << "Cannot open file " << label_file_path << std::endl;
        return 1;
    }
    if(!summary_file.is_open()){
        std::cerr << "Cannot open file " << summary_file_path << std::endl;
        return 1;
    }

    for(i = 0; i < int(class_label.size()); i++){
        double res = 0;
        for(j = 0; j < int(phi_out.size()); j++){
            if(class_label[i] == j){
                res = phi_out[j];
                break;
            }
        }
        phi_file << round(res * 1000.0 ) / 1000.0 << "\n";
        label_file << class_label[i] << "\n";
    }

    for(i = 0; i < int(phi_out.size()); i++){
        count = 0;
        double res = phi_out[i];
        for(j = 0; j < int(class_label.size()); j++){
            if(class_label[j] == i){
                count += 1;
            }
        }
        summary_file << i << "\t" << count << "\t" << round(res * 1000.0 ) / 1000.0 << "\n";
    }

    phi_file.close();
    label_file.close();
    summary_file.close();

    if (problematic_snvs.size()>0){
        std::string problematic_snvs_file_path = preliminary_folder + "/lam" + lambda_str + "_problematic_snvs.txt";
        ofstream problematic_snvs_file;

        problematic_snvs_file.open(problematic_snvs_file_path);

        if(!problematic_snvs_file.is_open()){
            std::cerr << "Cannot open file " << problematic_snvs_file_path  << std::endl;
            return 1;
        }

        for(i = 0; i < int(problematic_snvs.size()); i++){
            problematic_snvs_file << problematic_snvs[i] << "\n";
        }
    }

    return 0;
}

int postprocess_cuda_result_from_flags(
    int No_mutation,
    double Lambda,
    int k,
    int least_mut,
    double least_diff,
    const std::vector<double>& n,
    const std::vector<double>& phi_hat,
    const std::vector<unsigned char>& fused,
    const std::vector<long long>& pair_start,
    const std::vector<int>& problematic_snv_flags,
    std::string& preliminary_folder)
{
    std::vector<int> problematic_snvs;
    for(int i = 0; i < No_mutation; i++){
        if(problematic_snv_flags[i]) problematic_snvs.push_back(i);
    }

    std::vector<int> class_label;
    std::vector<int> group_size;
    build_cpu_order_class_labels_from_fused(No_mutation, fused, pair_start, class_label, group_size);

    if(double(cpu_order_minimum_group_size(group_size, No_mutation)) < least_mut){
        return kCliPPNeedFullEta;
    }

    return finalize_cuda_postprocess_from_labels(
        No_mutation,
        Lambda,
        k,
        least_diff,
        n,
        phi_hat,
        class_label,
        problematic_snvs,
        preliminary_folder);
}

int postprocess_cuda_result(
    int No_mutation,
    double Lambda,
    int k,
    int least_mut,
    double post_th,
    double least_diff,
    const std::vector<double>& n,
    const std::vector<double>& phi_hat,
    std::vector<double>& eta_new,
    const std::vector<long long>& pair_start,
    const std::vector<int>& problematic_snv_flags,
    std::string& preliminary_folder)
{
    int i, j, count;
    double temp, temp1;

    std::vector<int> problematic_snvs;
    for(i = 0; i < No_mutation; i++){
        if(problematic_snv_flags[i]) problematic_snvs.push_back(i);
    }

    for(long long pair_index = 0; pair_index < static_cast<long long>(eta_new.size()); pair_index++){
        if(fabs(eta_new[pair_index]) <= post_th) eta_new[pair_index] = 0.0;
    }

    std::vector<int> class_label(No_mutation);
    std::fill(class_label.begin(), class_label.end(), -1);
    class_label[0] = 0;

    std::vector<int> group_size;
    group_size.push_back(1);

    int labl = 1;

    for(i = 0; i < No_mutation; i++){
        for(j = 0; j < i; j++){
            const long long pair_index = pair_start[j] + i - j - 1;
            if(eta_new[pair_index] == 0.0){
                class_label[i] = class_label[j];
                group_size[class_label[j]] = group_size[class_label[j]] + 1;
                break;
            }
        }
        if(class_label[i] == -1){
            class_label[i] = labl;
            labl += 1;
            group_size.push_back(1);
        }
    }

    int temp_size = No_mutation;
    for(unsigned int ii = 0; ii < group_size.size(); ii++){
        if(group_size[ii] > 0 && temp_size > group_size[ii]) temp_size = group_size[ii];
    }

    std::vector<int> tmp_grp;
    for(unsigned int ii = 0; ii < group_size.size(); ii++){
        if(group_size[ii] == temp_size) {
            tmp_grp.push_back(ii);
        }
    }

    int refine = 0;

    if(double(temp_size) < least_mut) refine = 1;
    CudaPackedPairDiff diff(eta_new, pair_start);
    MatrixXd tmp_diff;
    if(refine == 1){
        tmp_diff.resize(No_mutation, 1);
    }
    while(refine == 1){
        refine = 0;

        std::vector<int> tmp_col;
        for(i = 0; i < No_mutation; i++){
            if(class_label[i] == tmp_grp[0]){
                tmp_col.push_back(i);
            }
        }

        for(unsigned int ii = 0; ii < tmp_col.size(); ii++){
            if(tmp_col[ii] != 0 && tmp_col[ii] != No_mutation - 1){
                for(int jj = 0; jj < tmp_col[ii]; jj++){
                    tmp_diff(jj, 0) = fabs(diff.value(jj, tmp_col[ii]));
                }
                tmp_diff(tmp_col[ii], 0) = 100.0;

                for(int jj = tmp_col[ii] + 1; jj < No_mutation; jj++){
                    tmp_diff(jj, 0) = fabs(diff.value(tmp_col[ii], jj));
                }

                for(unsigned int jj = 0; jj < tmp_col.size(); jj++){
                    tmp_diff(tmp_col[jj], 0) = tmp_diff(tmp_col[jj], 0) + 100.0;
                }

                for(unsigned int jj = 0; jj < tmp_col.size(); jj++){
                    diff.set(jj, tmp_col[ii], tmp_diff(jj, 0));
                }

                for(int jj = tmp_col[ii] + 1; jj < No_mutation; jj++){
                    diff.set(tmp_col[ii], jj, tmp_diff(jj, 0));
                }
            }else {
                if(tmp_col[ii] == 0){
                    tmp_diff(0, 0) = 100.0;
                    for(int jj = 1; jj < No_mutation; jj++){
                        tmp_diff(jj, 0) = fabs(diff.value(0, jj));
                    }
                    for(unsigned int jj = 0; jj < tmp_col.size(); jj++){
                        tmp_diff(tmp_col[jj], 0) = tmp_diff(tmp_col[jj], 0) + 100.0;
                    }
                    for(int jj = 1; jj < No_mutation; jj++){
                        diff.set(0, jj, tmp_diff(jj, 0));
                    }
                }else{
                    for(int jj = 0; jj < No_mutation - 1; jj++){
                        tmp_diff(jj, 0) = fabs(diff.value(jj, No_mutation - 1));
                    }
                    tmp_diff(No_mutation - 1, 0) = 100.0;
                    for(unsigned int jj = 0; jj < tmp_col.size(); jj++){
                        tmp_diff(tmp_col[jj], 0) = tmp_diff(tmp_col[jj], 0) + 100.0;
                    }
                    for(int jj = 0; jj < No_mutation - 1; jj++){
                        diff.set(jj, No_mutation - 1, tmp_diff(jj, 0));
                    }
                }
            }

            int ind = 0;
            temp = 100000000.0;
            for(int jj= 0; jj < No_mutation; jj++){
                if(tmp_diff(jj, 0) < temp){
                    temp = tmp_diff(jj, 0);
                    ind = jj;
                }
            }
            group_size[class_label[tmp_col[ii]]] = group_size[class_label[tmp_col[ii]]] - 1;
            class_label[tmp_col[ii]] = class_label[ind];
            group_size[class_label[tmp_col[ii]]] = group_size[class_label[tmp_col[ii]]] + 1;
        }

        temp_size = No_mutation;
        for(unsigned int ii = 0; ii < group_size.size(); ii++){
            if(group_size[ii] > 0 && temp_size > group_size[ii]) temp_size = group_size[ii];
        }

        tmp_grp.clear();
        for(unsigned int ii = 0; ii < group_size.size(); ii++){
            if(group_size[ii] == temp_size) {
                tmp_grp.push_back(ii);
            }
        }
        refine = 0;
        if (temp_size < least_mut) refine = 1;
    }

    std::vector<int> labels = class_label;
    std::sort(labels.begin(), labels.end());
    vector<int>::iterator ip;
    ip = std::unique(labels.begin(), labels.end());
    labels.resize(std::distance(labels.begin(), ip));

    std::vector<double> phi_out(labels.size());
    std::fill(phi_out.begin(), phi_out.end(), 0);

    for(unsigned int ii = 0; ii < labels.size(); ii++){
        std::vector<int> ind;
        ind.reserve(No_mutation);

        for(i = 0; i < No_mutation; i++) {
            if(class_label[i] == labels[ii]) ind.push_back(i);
        }

        temp = 0;
        temp1 = 0;
        for(i = 0; i < int(ind.size()); i++) {
            class_label[ind[i]] = ii;
            temp += phi_hat[ind[i]] * n[ind[i]];
            temp1 += n[ind[i]];
        }

        phi_out[ii] = temp / temp1;
    }

    if(labels.size() > 1){
        std::vector<double> sort_phi = phi_out;
        std::sort(sort_phi.begin(), sort_phi.end());
        std::vector<double> phi_diff;

        for(i = 1; i < int(sort_phi.size()); i++) {
            phi_diff.push_back(sort_phi[i] - sort_phi[i -1]);
        }

        double min_val = 10000000000;
        int min_ind = 0;

        for(i = 0; i < int(phi_diff.size()); i++){
            if (min_val > phi_diff[i]){
                min_val = phi_diff[i];
                min_ind = i;
            }
        }

        std::vector<double> min_val_vector;
        min_val_vector.push_back(min_val);

        while(min_val < least_diff){
            std::vector<int> combine_ind, combine_to_ind;
            for(i = 0; i < int(phi_out.size()); i++) {
                if(phi_out[i] == sort_phi[min_ind]) combine_ind.push_back(i);
                if(phi_out[i] == sort_phi[min_ind + 1]) combine_to_ind.push_back(i);
            }

            for(i = 0; i < int(class_label.size()); i++) {
                if(class_label[i] == combine_ind[0]) class_label[i] = combine_to_ind[0];
            }

            labels.clear();
            labels = class_label;
            std::sort(labels.begin(), labels.end());
            ip = std::unique(labels.begin(), labels.end());
            labels.resize(std::distance(labels.begin(), ip));

            phi_out.clear();
            for(i = 0; i < int(labels.size()); i++) {
                phi_out.push_back(0.0);
            }

            for(unsigned int ii = 0; ii < labels.size(); ii++){
                std::vector<int> ind;
                ind.reserve(No_mutation);

                for(i = 0; i < No_mutation; i++) {
                    if(class_label[i] == labels[ii]) ind.push_back(i);
                }
                temp = 0;
                temp1 = 0;

                for(i = 0; i < int(ind.size()); i++) {
                    class_label[ind[i]] = ii;
                    temp += phi_hat[ind[i]] * n[ind[i]];
                    temp1 += n[ind[i]];
                }
                phi_out[ii] = temp / temp1;
            }

            if(labels.size() == 1) break;
            else{
                sort_phi = phi_out;
                std::sort(sort_phi.begin(), sort_phi.end());
                phi_diff.clear();
                for(i = 1; i < int(sort_phi.size()); i++) {
                    phi_diff.push_back(sort_phi[i] - sort_phi[i -1]);
                }

                min_val = 10000000000;
                for(i = 0; i < int(phi_diff.size()); i++){
                    if (min_val > phi_diff[i]){
                        min_val = phi_diff[i];
                        min_ind = i;
                    }
                }

                min_val_vector.push_back(min_val);

                if(min_val_vector.size() > 4){
                    int unique_min_val = 1;
                    for(double value:min_val_vector){
                        if(value != min_val) unique_min_val += 1;
                    }
                    if(unique_min_val == 1){
                        std::cout << "Lambda: " << Lambda << "\titeration: " << k << "\tfailed" << std::endl;
                        return -1;
                    }
                }
            }
        }
    }

    ofstream phi_file, label_file, summary_file;

    std::string lambda_str = lambda_to_string(Lambda);

    std::string phi_file_path = preliminary_folder + "/lam" + lambda_str + "_phi.txt";
    std::string label_file_path = preliminary_folder + "/lam" + lambda_str + "_label.txt";
    std::string summary_file_path = preliminary_folder + "/lam" + lambda_str + "_summary_table.txt";

    phi_file.open(phi_file_path);
    label_file.open(label_file_path);
    summary_file.open(summary_file_path);

    if(!phi_file.is_open()){
        std::cerr << "Cannot open file " << phi_file_path  << std::endl;
        return 1;
    }

    if(!label_file.is_open()){
        std::cerr << "Cannot open file " << label_file_path << std::endl;
        return 1;
    }
    if(!summary_file.is_open()){
        std::cerr << "Cannot open file " << summary_file_path << std::endl;
        return 1;
    }

    for(i = 0; i < int(class_label.size()); i++){
        double res = 0;
        for(j = 0; j < int(phi_out.size()); j++){
            if(class_label[i] == j){
                res = phi_out[j];
                break;
            }
        }
        phi_file << round(res * 1000.0 ) / 1000.0 << "\n";
        label_file << class_label[i] << "\n";
    }

    for(i = 0; i < int(phi_out.size()); i++){
        count = 0;
        double res = phi_out[i];
        for(j = 0; j < int(class_label.size()); j++){
            if(class_label[j] == i){
                count += 1;
            }
        }
        summary_file << i << "\t" << count << "\t" << round(res * 1000.0 ) / 1000.0 << "\n";
    }

    phi_file.close();
    label_file.close();
    summary_file.close();

    if (problematic_snvs.size()>0){
        std::string problematic_snvs_file_path = preliminary_folder + "/lam" + lambda_str + "_problematic_snvs.txt";
        ofstream problematic_snvs_file;

        problematic_snvs_file.open(problematic_snvs_file_path);

        if(!problematic_snvs_file.is_open()){
            std::cerr << "Cannot open file " << problematic_snvs_file_path  << std::endl;
            return 1;
        }

        for(i = 0; i < int(problematic_snvs.size()); i++){
            problematic_snvs_file << problematic_snvs[i] << "\n";
        }
    }

    return 0;
}

int CliPPIndividualCUDA(
    int No_mutation,
    const std::vector<double>& n,
    const CudaStaticInputs& static_inputs,
    CliPPCudaWorkspace& workspace,
    double Lambda,
    double alpha,
    double rho,
    double gamma,
    int Run_limit,
    double precision,
    int control_large,
    int least_mut,
    double post_th,
    double least_diff,
    std::string& preliminary_folder,
    const std::vector<double>& phi_hat,
    const std::vector<long long>& pair_start,
    double scale_parameter,
    cudaStream_t stream)
{
    const long long pair_count = static_cast<long long>(No_mutation) * (No_mutation - 1) / 2;
    CudaProgram& program = cuda_program();

    if(pair_count > 0 && static_cast<unsigned long long>(pair_count) > std::numeric_limits<size_t>::max() / sizeof(double) / 2ULL){
        std::cerr << "CUDA memory size overflow for pair arrays." << std::endl;
        return kCudaUnavailable;
    }

    GpuBuffer<double>& d_w_new = workspace.d_w_new;
    GpuBuffer<double>& d_w_old = workspace.d_w_old;
    GpuBuffer<double>& d_A = workspace.d_A;
    GpuBuffer<double>& d_B = workspace.d_B;
    GpuBuffer<double>& d_linear = workspace.d_linear;
    GpuBuffer<double>& d_Minv = workspace.d_Minv;
    GpuBuffer<double>& d_dot_terms = workspace.d_dot_terms;
    GpuBuffer<double>& d_eta = workspace.d_eta;
    GpuBuffer<double>& d_tau = workspace.d_tau;
    GpuBuffer<double>& d_row_residual_max = workspace.d_row_residual_max;
    GpuBuffer<double>& d_reduce_tmp_a = workspace.d_reduce_tmp_a;
    GpuBuffer<double>& d_reduce_tmp_b = workspace.d_reduce_tmp_b;
    GpuBuffer<int>& d_problematic_flags = workspace.d_problematic_flags;
    GpuBuffer<unsigned char>& d_fused = workspace.d_fused;

    const GpuBuffer<double>& d_n = static_inputs.d_n;
    const GpuBuffer<double>& d_minor = static_inputs.d_minor;
    const GpuBuffer<double>& d_total = static_inputs.d_total;
    const GpuBuffer<double>& d_theta_hat = static_inputs.d_theta_hat;
    const GpuBuffer<double>& d_phi_hat = static_inputs.d_phi_hat;
    const GpuBuffer<double>& d_coef = static_inputs.d_coef;
    const GpuBuffer<double>& d_wcut = static_inputs.d_wcut;

    if(No_mutation > 0){
        int n_mut = No_mutation;
        int block = kCudaBlockSize;
        unsigned int grid = blocks_for(No_mutation);
        double* phi_ptr = d_phi_hat.get();
        double* w_ptr = d_w_new.get();
        void* args[] = { &n_mut, &phi_ptr, &scale_parameter, &control_large, &w_ptr };
        program.launch("init_w_kernel", dim3(grid), dim3(block), args, 0, stream);
    }

    if(pair_count > 0){
        int n_mut = No_mutation;
        double* w_ptr = d_w_new.get();
        double* eta_ptr = d_eta.get();
        void* init_eta_args[] = { &n_mut, &w_ptr, &eta_ptr };
        program.launch("init_eta_kernel", dim3(No_mutation - 1), dim3(kCudaBlockSize), init_eta_args, 0, stream);

        long long count_ll = pair_count;
        double one = 1.0;
        double* tau_ptr = d_tau.get();
        void* fill_tau_args[] = { &count_ll, &tau_ptr, &one };
        program.launch("fill_double_kernel", dim3(blocks_for(pair_count)), dim3(kCudaBlockSize), fill_tau_args, 0, stream);
    }

    if(No_mutation > 0){
        CLIPP_CUDA_CHECK(cudaMemsetAsync(
            d_problematic_flags.get(),
            0,
            static_cast<size_t>(No_mutation) * sizeof(int),
            stream));
    }

    double residual = 100.0;
    double trace_g = 0.0;
    int k = 0;

    while(residual > precision && k < Run_limit){
        k += 1;
#if CLIPP_VERBOSE
        std::cout << "\rLambda: " << Lambda << "\titeration: " << k << "\tresidual: " << residual << std::flush;
#endif

        d_w_old.swap(d_w_new);

        {
            int n_mut = No_mutation;
            double* w_old = d_w_old.get();
            double* n_ptr = d_n.get();
            double* minor_ptr = d_minor.get();
            double* total_ptr = d_total.get();
            double* theta_hat_ptr = d_theta_hat.get();
            double* coef_ptr = d_coef.get();
            double* wcut_ptr = d_wcut.get();
            double* A_ptr = d_A.get();
            double* B_ptr = d_B.get();
            int* flags_ptr = d_problematic_flags.get();
            void* args[] = { &n_mut, &w_old, &n_ptr, &minor_ptr, &total_ptr, &theta_hat_ptr, &coef_ptr, &wcut_ptr, &A_ptr, &B_ptr, &flags_ptr };
            program.launch("compute_ab_kernel", dim3(blocks_for(No_mutation)), dim3(kCudaBlockSize), args, 0, stream);
        }

        {
            int n_mut = No_mutation;
            double* eta_old = d_eta.get();
            double* tau_old = d_tau.get();
            double* A_ptr = d_A.get();
            double* B_ptr = d_B.get();
            double* linear_ptr = d_linear.get();
            double* Minv_ptr = d_Minv.get();
            double* dot_ptr = d_dot_terms.get();
            void* args[] = { &n_mut, &alpha, &eta_old, &tau_old, &A_ptr, &B_ptr, &linear_ptr, &Minv_ptr, &dot_ptr };
            program.launch("compute_linear_minv_kernel", dim3(No_mutation), dim3(kCudaBlockSize), args, kCudaBlockSize * sizeof(double), stream);
        }

        double sum_minv = reduce_device_sum(program, d_Minv.get(), No_mutation, d_reduce_tmp_a, d_reduce_tmp_b, stream);
        trace_g = -alpha * sum_minv;
        if(!std::isfinite(trace_g)){
            std::cout << "Lambda: " << Lambda << "\titeration: " << k << "\tnon-finite trace_g" << std::endl;
            return -1;
        }

        double minv_linear_dot = reduce_device_sum(program, d_dot_terms.get(), No_mutation, d_reduce_tmp_a, d_reduce_tmp_b, stream);
        const double rank1_denom = 1.0 + trace_g;
        if(!std::isfinite(rank1_denom) || std::fabs(rank1_denom) < 1e-12){
            std::cerr << "Invalid rank1 denominator at Lambda: " << Lambda
                      << ", iteration: " << k
                      << ", trace_g: " << trace_g << std::endl;
            return -1;
        }
        double rank1_scale = alpha / rank1_denom;

        {
            int n_mut = No_mutation;
            double* Minv_ptr = d_Minv.get();
            double* linear_ptr = d_linear.get();
            double* w_new = d_w_new.get();
            void* args[] = { &n_mut, &Minv_ptr, &linear_ptr, &minv_linear_dot, &rank1_scale, &control_large, &w_new };
            program.launch("update_w_kernel", dim3(blocks_for(No_mutation)), dim3(kCudaBlockSize), args, 0, stream);
        }

        if(pair_count > 0){
            double mcp_denom = 1.0 - 1.0 / ((gamma - 1.0) * alpha);
            if(!std::isfinite(mcp_denom) || std::fabs(mcp_denom) < 1e-12){
                std::cerr << "Invalid MCP denominator at Lambda: " << Lambda
                          << ", iteration: " << k
                          << ", alpha: " << alpha
                          << ", gamma: " << gamma << std::endl;
                return -1;
            }

            int n_mut = No_mutation;
            double* w_new = d_w_new.get();
            double* eta_ptr = d_eta.get();
            double* tau_ptr = d_tau.get();
            double* residual_rows = d_row_residual_max.get();
            void* args[] = { &n_mut, &alpha, &Lambda, &gamma, &mcp_denom, &w_new, &eta_ptr, &tau_ptr, &residual_rows };
            program.launch("update_eta_tau_kernel", dim3(No_mutation), dim3(kCudaBlockSize), args, kCudaBlockSize * sizeof(double), stream);
            residual = reduce_device_max(program, d_row_residual_max.get(), No_mutation, d_reduce_tmp_a, d_reduce_tmp_b, stream);
        }else{
            residual = -100000.0;
        }

        alpha = alpha * rho;
    }

#if CLIPP_VERBOSE
    std::cout << std::endl;
#endif

    std::vector<int> problematic_flags;
    copy_to_host(problematic_flags, d_problematic_flags, stream);

    if(pair_count > 0){
        long long count_ll = pair_count;
        double* eta_ptr = d_eta.get();
        unsigned char* fused_ptr = d_fused.get();
        void* args[] = { &count_ll, &post_th, &eta_ptr, &fused_ptr };
        program.launch("threshold_eta_flag_kernel", dim3(blocks_for(pair_count)), dim3(kCudaBlockSize), args, 0, stream);

        std::vector<unsigned char> fused_host;
        copy_to_host(fused_host, d_fused, stream);

        const int flag_postprocess_status = postprocess_cuda_result_from_flags(
            No_mutation,
            Lambda,
            k,
            least_mut,
            least_diff,
            n,
            phi_hat,
            fused_host,
            pair_start,
            problematic_flags,
            preliminary_folder);

        if(flag_postprocess_status != kCliPPNeedFullEta){
            return flag_postprocess_status;
        }
    }

    std::vector<double> eta_host;
    if(pair_count > 0){
        copy_to_host(eta_host, d_eta, stream);
    }

    return postprocess_cuda_result(
        No_mutation,
        Lambda,
        k,
        least_mut,
        post_th,
        least_diff,
        n,
        phi_hat,
        eta_host,
        pair_start,
        problematic_flags,
        preliminary_folder);
}

} // namespace

extern "C" int CliPPWarmupCUDA()
{
#ifdef USE_CUDA
    try{
        int device_count = 0;
        CLIPP_CUDA_CHECK(cudaGetDeviceCount(&device_count));
        if(device_count <= 0){
            return kCudaUnavailable;
        }

        CLIPP_CUDA_CHECK(cudaSetDevice(0));
        CLIPP_CUDA_DRIVER_CHECK(cuInit(0));
        CLIPP_CUDA_CHECK(cudaFree(nullptr));

        CudaProgram& program = cuda_program();
        GpuBuffer<double> d_values(1);
        long long count = 1;
        double value = 0.0;
        double* values = d_values.get();
        void* args[] = { &count, &values, &value };
        program.launch("fill_double_kernel", dim3(1), dim3(kCudaBlockSize), args);
        CLIPP_CUDA_CHECK(cudaDeviceSynchronize());
        return kCliPPOk;
    }catch(const std::exception& err){
#if CLIPP_VERBOSE
        std::cerr << "CUDA warmup failed: " << err.what() << std::endl;
#endif
        return kCudaUnavailable;
    }
#else
    return kCudaUnavailable;
#endif
}

int CliPPCUDA(int No_mutation, int* c_r, int *c_n, int *c_minor, int *c_total, double ploidy, double* Lambda_list, int Lambda_num, double alpha, double rho, double gamma, int Run_limit, double precision, int control_large, int least_mut, double post_th, double least_diff, double* c_coef_1d, double* c_wcut_1d, double purity, char* preliminary)
{
    bool cuda_work_started = false;
    try{
        int validation_status = validate_clipp_inputs(
            No_mutation,
            c_r,
            c_n,
            c_minor,
            c_total,
            Lambda_list,
            Lambda_num,
            ploidy,
            purity,
            alpha,
            rho,
            gamma,
            Run_limit,
            precision,
            control_large,
            least_mut,
            post_th,
            least_diff,
            c_coef_1d,
            c_wcut_1d,
            preliminary);
        if(validation_status != kCliPPOk){
            return validation_status;
        }

        int device_count = 0;
        CLIPP_CUDA_CHECK(cudaGetDeviceCount(&device_count));
        if(device_count <= 0){
            std::cerr << "CUDA backend requested, but no CUDA device is available." << std::endl;
            return kCudaUnavailable;
        }

        CLIPP_CUDA_CHECK(cudaSetDevice(0));
        CLIPP_CUDA_DRIVER_CHECK(cuInit(0));
        CLIPP_CUDA_CHECK(cudaFree(nullptr));

        std::string preliminary_folder(preliminary);

#ifdef _OPENMP
        Eigen::setNbThreads(1);
#endif

        std::vector<double> n(No_mutation);
        std::vector<double> minor_(No_mutation);
        std::vector<double> total(No_mutation);
        std::vector<double> theta_hat(No_mutation);
        std::vector<double> phi_hat(No_mutation);

        for(int i = 0; i < No_mutation; i++){
            n[i] = double(c_n[i]);
            minor_[i] = double(c_minor[i]);
            total[i] = double(c_total[i]);
            theta_hat[i] = double(c_r[i]) / n[i];
            phi_hat[i] = theta_hat[i] * ((ploidy - purity * ploidy + purity * total[i]) / minor_[i]);
        }

        double phi_hat_max = 0.0;
        if(No_mutation > 0){
            phi_hat_max = *std::max_element(phi_hat.begin(), phi_hat.end());
        }
        double scale_parameter = max(1.0, phi_hat_max);

        std::vector<double> wcut_1d(No_mutation * 2);
        std::vector<double> coef_1d(No_mutation * 6);

        for(int i = 0; i < No_mutation * 2; i++){
            wcut_1d[i] = c_wcut_1d[i];
        }

        for(int i = 0; i < No_mutation * 6; i++){
            coef_1d[i] = c_coef_1d[i];
        }

        std::vector<long long> pair_start(No_mutation, 0);
        long long p_count = 0;
        for(int i = 0; i < No_mutation - 1; i++){
            pair_start[i] = p_count;
            p_count += No_mutation - i - 1;
        }
        if(No_mutation > 0) pair_start[No_mutation - 1] = p_count;

        (void)cuda_program();

        CudaStream static_stream;
        CudaStaticInputs static_inputs;
        static_inputs.reset(
            n,
            minor_,
            total,
            theta_hat,
            phi_hat,
            coef_1d,
            wcut_1d,
            static_stream.get());

        size_t free_bytes = 0;
        size_t total_bytes = 0;
        CLIPP_CUDA_CHECK(cudaMemGetInfo(&free_bytes, &total_bytes));

        const size_t per_lambda_bytes = estimate_cuda_lambda_bytes(No_mutation);
        int memory_parallelism = 1;
        if(per_lambda_bytes != std::numeric_limits<size_t>::max() && per_lambda_bytes > 0){
            memory_parallelism = static_cast<int>(
                std::max<size_t>(1, static_cast<size_t>(0.70 * static_cast<double>(free_bytes)) / per_lambda_bytes));
        }
        memory_parallelism = std::max(1, std::min(memory_parallelism, Lambda_num));

        int lambda_parallelism = 1;
        if(const char* requested_parallelism = std::getenv("CLIPP_CUDA_LAMBDA_PARALLELISM")){
            const int requested = std::atoi(requested_parallelism);
            if(requested > 0){
                lambda_parallelism = std::max(1, std::min({requested, Lambda_num, memory_parallelism}));
            }
        }

#if CLIPP_VERBOSE
        std::cerr << "CUDA lambda parallelism: " << lambda_parallelism
                  << " of " << Lambda_num << " lambdas." << std::endl;
#endif

        std::vector<CliPPCudaWorkspace> workspaces(lambda_parallelism);
        for(int i = 0; i < lambda_parallelism; ++i){
            workspaces[i].reset(No_mutation, p_count);
        }

        for(int batch_start = 0; batch_start < Lambda_num; batch_start += lambda_parallelism){
            const int batch_count = std::min(lambda_parallelism, Lambda_num - batch_start);
            std::atomic<int> batch_status{kCliPPOk};
            cuda_work_started = true;

#ifdef _OPENMP
#pragma omp parallel for num_threads(batch_count) schedule(static)
#endif
            for(int local_index = 0; local_index < batch_count; ++local_index){
                if(batch_status.load(std::memory_order_relaxed) != kCliPPOk) continue;

                const int lambda_index = batch_start + local_index;
                const double Lambda = Lambda_list[lambda_index];

                try{
                    CLIPP_CUDA_CHECK(cudaSetDevice(0));
                    CudaStream stream;
                    int rc = CliPPIndividualCUDA(
                        No_mutation,
                        n,
                        static_inputs,
                        workspaces[local_index],
                        Lambda,
                        alpha,
                        rho,
                        gamma,
                        Run_limit,
                        precision,
                        control_large,
                        least_mut,
                        post_th,
                        least_diff,
                        preliminary_folder,
                        phi_hat,
                        pair_start,
                        scale_parameter,
                        stream.get());
                    CLIPP_CUDA_CHECK(cudaStreamSynchronize(stream.get()));
                    if(rc != kCliPPOk){
                        batch_status.store(rc, std::memory_order_relaxed);
                    }
                }catch(const std::exception& err){
#ifdef _OPENMP
#pragma omp critical
#endif
                    {
                        std::cerr << "CUDA lambda " << lambda_index << " failed: " << err.what() << std::endl;
                    }
                    batch_status.store(kCudaFailedAfterWork, std::memory_order_relaxed);
                }
            }

            const int rc = batch_status.load(std::memory_order_relaxed);
            if(rc != kCliPPOk){
                return batch_start == 0 ? rc : kCudaFailedAfterWork;
            }
        }

        CLIPP_CUDA_CHECK(cudaDeviceSynchronize());
        return kCliPPOk;
    }catch(const std::exception& err){
        std::cerr << "CUDA backend failed: " << err.what() << std::endl;
        return cuda_work_started ? kCudaFailedAfterWork : kCudaUnavailable;
    }
}
#endif
