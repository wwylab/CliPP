R"CUDA(
#define CLIPP_DBL_MAX 1.79769313486231570814527423731704357e+308

__device__ __forceinline__ long long pair_start_ll(int i, int n)
{
    return (long long)i * (long long)n - ((long long)i * (long long)(i + 1)) / 2LL;
}

__device__ __forceinline__ long long pair_index_ll(int i, int j, int n)
{
    return pair_start_ll(i, n) + (long long)(j - i - 1);
}

__device__ __forceinline__ double expit_dev(double x)
{
    return 1.0 / (1.0 + exp(-x));
}

__device__ __forceinline__ double logit_dev(double x)
{
    return log(x / (1.0 - x));
}

__device__ __forceinline__ double st_dev(double x, double lam)
{
    const double val = fabs(x) - lam;
    if(x == 0.0 || val <= 0.0) return 0.0;
    return x > 0.0 ? val : -val;
}

extern "C" __global__ void fill_double_kernel(long long count, double* values, double value)
{
    const long long idx = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if(idx < count) values[idx] = value;
}


extern "C" __global__ void init_w_kernel(
    int n_mut,
    const double* __restrict__ phi_hat,
    double scale_parameter,
    int control_large,
    double* __restrict__ w_new)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if(i >= n_mut) return;

    double temp = phi_hat[i] / scale_parameter;
    const double upper = expit_dev((double)control_large);
    const double lower = expit_dev(-(double)control_large);

    if(temp > upper) temp = upper;
    if(temp < lower) temp = lower;

    temp = logit_dev(temp);

    if(temp > control_large) temp = control_large;
    if(temp < -control_large) temp = -control_large;

    w_new[i] = temp;
}

extern "C" __global__ void init_eta_kernel(
    int n_mut,
    const double* __restrict__ w_new,
    double* __restrict__ eta_new)
{
    const int i = blockIdx.x;
    if(i >= n_mut - 1) return;

    const long long start = pair_start_ll(i, n_mut);
    const int row_len = n_mut - i - 1;

    for(int offset = threadIdx.x; offset < row_len; offset += blockDim.x){
        const int j = i + 1 + offset;
        eta_new[start + offset] = w_new[i] - w_new[j];
    }
}

extern "C" __global__ void compute_ab_kernel(
    int n_mut,
    const double* __restrict__ w_old,
    const double* __restrict__ n,
    const double* __restrict__ minor,
    const double* __restrict__ total,
    const double* __restrict__ theta_hat,
    const double* __restrict__ coef,
    const double* __restrict__ wcut,
    double* __restrict__ A,
    double* __restrict__ B,
    int* __restrict__ problematic_flags)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if(i >= n_mut) return;

    const double wi = w_old[i];
    const double exp_w = exp(wi);

    double theta = exp_w * minor[i] / (2.0 + exp_w * total[i]);

    double tag1, tag2, tag3, tag4;

    if(wi <= wcut[i * 2 + 0]){
        tag1 = 1.0;
        tag3 = 0.0;
    }else{
        tag1 = 0.0;
        tag3 = 1.0;
    }

    if(wi >= wcut[i * 2 + 1]){
        tag2 = 1.0;
        tag4 = 0.0;
    }else{
        tag2 = 0.0;
        tag4 = 1.0;
    }

    if(theta >= 1.0){
        theta = 0.99;
        problematic_flags[i] = 1;
    }

    const double sqrt_n = sqrt(n[i]);
    const double sqrt_theta = sqrt(theta * (1.0 - theta));

    A[i] =
        sqrt_n *
        (tag1 * coef[i * 6 + 1] + tag2 * coef[i * 6 + 5] + tag3 * tag4 * coef[i * 6 + 3] - theta_hat[i]) /
        sqrt_theta;

    B[i] =
        sqrt_n *
        (tag1 * coef[i * 6 + 0] + tag2 * coef[i * 6 + 4] + tag3 * tag4 * coef[i * 6 + 2]) /
        sqrt_theta;
}


extern "C" __global__ void compute_linear_minv_kernel(
    int n_mut,
    double alpha,
    const double* __restrict__ eta_old,
    const double* __restrict__ tau_old,
    const double* __restrict__ A,
    const double* __restrict__ B,
    double* __restrict__ linear,
    double* __restrict__ Minv,
    double* __restrict__ dot_terms)
{
    extern __shared__ double shared[];

    const int i = blockIdx.x;
    if(i >= n_mut) return;

    double local_sum = 0.0;

    for(int j = threadIdx.x; j < n_mut; j += blockDim.x){
        if(j < i){
            const long long p = pair_index_ll(j, i, n_mut);
            const double pair_value = alpha * eta_old[p] + tau_old[p];
            local_sum -= pair_value;
        }else if(j > i){
            const long long p = pair_index_ll(i, j, n_mut);
            const double pair_value = alpha * eta_old[p] + tau_old[p];
            local_sum += pair_value;
        }
    }

    shared[threadIdx.x] = local_sum;
    __syncthreads();

    for(int stride = blockDim.x / 2; stride > 0; stride >>= 1){
        if(threadIdx.x < stride) shared[threadIdx.x] += shared[threadIdx.x + stride];
        __syncthreads();
    }

    if(threadIdx.x == 0){
        const double bi = B[i];
        const double lin = shared[0] - bi * A[i];
        const double mi = 1.0 / (bi * bi + (double)n_mut * alpha);

        linear[i] = lin;
        Minv[i] = mi;
        dot_terms[i] = mi * lin;
    }
}


extern "C" __global__ void update_w_kernel(
    int n_mut,
    const double* __restrict__ Minv,
    const double* __restrict__ linear,
    double minv_linear_dot,
    double rank1_scale,
    int control_large,
    double* __restrict__ w_new)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if(i >= n_mut) return;

    double wi = Minv[i] * linear[i] + rank1_scale * Minv[i] * minv_linear_dot;

    if(wi > control_large) wi = control_large;
    if(wi < -control_large) wi = -control_large;

    w_new[i] = wi;
}

extern "C" __global__ void update_eta_tau_kernel(
    int n_mut,
    double alpha,
    double Lambda,
    double gamma,
    double mcp_denom,
    const double* __restrict__ w_new,
    double* __restrict__ eta,
    double* __restrict__ tau,
    double* __restrict__ row_residual_max)
{
    extern __shared__ double shared[];

    const int i = blockIdx.x;
    if(i >= n_mut - 1){
        if(threadIdx.x == 0) row_residual_max[i] = -CLIPP_DBL_MAX;
        return;
    }

    const long long start = pair_start_ll(i, n_mut);
    const int row_len = n_mut - i - 1;
    double local_max = -CLIPP_DBL_MAX;

    for(int offset = threadIdx.x; offset < row_len; offset += blockDim.x){
        const int j = i + 1 + offset;
        const long long p = start + offset;
        const double pair_diff = w_new[i] - w_new[j];
        const double tau_old = tau[p];
        const double temp = pair_diff - tau_old / alpha;

        double tag1, tag2, tag3, tag4;

        if(fabs(temp) > gamma * Lambda){
            tag1 = 1.0;
            tag3 = 0.0;
        }else{
            tag1 = 0.0;
            tag3 = 1.0;
        }

        if(fabs(temp) < (Lambda + Lambda / alpha)){
            tag2 = 1.0;
            tag4 = 0.0;
        }else{
            tag2 = 0.0;
            tag4 = 1.0;
        }

        const double eta_value =
            temp * tag1
            + st_dev(temp, Lambda / alpha) * tag2
            + st_dev(temp, gamma * Lambda / ((gamma - 1.0) * alpha)) /
                  mcp_denom * tag3 * tag4;

        eta[p] = eta_value;

        const double residual_i = pair_diff - eta_value;
        tau[p] = tau_old - alpha * residual_i;

        if(residual_i > local_max) local_max = residual_i;
    }

    shared[threadIdx.x] = local_max;
    __syncthreads();

    for(int stride = blockDim.x / 2; stride > 0; stride >>= 1){
        if(threadIdx.x < stride){
            const double other = shared[threadIdx.x + stride];
            if(other > shared[threadIdx.x]) shared[threadIdx.x] = other;
        }
        __syncthreads();
    }

    if(threadIdx.x == 0) row_residual_max[i] = shared[0];
}


extern "C" __global__ void threshold_eta_flag_kernel(
    long long pair_count,
    double post_th,
    double* __restrict__ eta_new,
    unsigned char* __restrict__ fused)
{
    const long long p = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if(p >= pair_count) return;

    if(fabs(eta_new[p]) <= post_th){
        eta_new[p] = 0.0;
        fused[p] = 1;
    }else{
        fused[p] = 0;
    }
}

extern "C" __global__ void reduce_sum_kernel(
    const double* __restrict__ input,
    double* __restrict__ output,
    long long n)
{
    extern __shared__ double shared[];
    const unsigned int tid = threadIdx.x;
    const long long grid_stride = (long long)gridDim.x * blockDim.x * 2LL;
    double sum = 0.0;

    for(long long idx = (long long)blockIdx.x * blockDim.x * 2LL + tid; idx < n; idx += grid_stride){
        sum += input[idx];
        const long long idx2 = idx + blockDim.x;
        if(idx2 < n) sum += input[idx2];
    }

    shared[tid] = sum;
    __syncthreads();

    for(unsigned int stride = blockDim.x / 2; stride > 0; stride >>= 1){
        if(tid < stride) shared[tid] += shared[tid + stride];
        __syncthreads();
    }

    if(tid == 0) output[blockIdx.x] = shared[0];
}

extern "C" __global__ void reduce_max_kernel(
    const double* __restrict__ input,
    double* __restrict__ output,
    long long n)
{
    extern __shared__ double shared[];
    const unsigned int tid = threadIdx.x;
    const long long grid_stride = (long long)gridDim.x * blockDim.x * 2LL;
    double value = -CLIPP_DBL_MAX;

    for(long long idx = (long long)blockIdx.x * blockDim.x * 2LL + tid; idx < n; idx += grid_stride){
        const double v1 = input[idx];
        if(v1 > value) value = v1;
        const long long idx2 = idx + blockDim.x;
        if(idx2 < n){
            const double v2 = input[idx2];
            if(v2 > value) value = v2;
        }
    }

    shared[tid] = value;
    __syncthreads();

    for(unsigned int stride = blockDim.x / 2; stride > 0; stride >>= 1){
        if(tid < stride){
            const double other = shared[tid + stride];
            if(other > shared[tid]) shared[tid] = other;
        }
        __syncthreads();
    }

    if(tid == 0) output[blockIdx.x] = shared[0];
}
)CUDA"
