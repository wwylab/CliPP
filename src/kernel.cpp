#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <fstream>
#include <algorithm>
//#include <chrono>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#ifdef _OPENMP 
#include <omp.h>
#endif

using Eigen::MatrixXd;
using Eigen::MatrixXi;
using Eigen::VectorXi;
using Eigen::VectorXd;
using Eigen::ArrayXi;
using Eigen::SparseMatrix;

using namespace std;

#define expit(a) (1/(1+exp(-(a))))
#define logit(a) (log((a)/(1-(a))))

double ST(double x, double lam){
    double val = fabs(x) - lam;
    if (x == 0 || val <= 0) return (0.0);
    if (x > 0) return (val);
    if (x < 0) return (-val);
    return (0.0);
}

int  CliPPIndividual(int, MatrixXd &, MatrixXd &, MatrixXd &, MatrixXd &, double, double, double, double, double, int, double, int, int, double, double, MatrixXd &, MatrixXd &, double, std::string &, MatrixXd &, MatrixXd &, MatrixXi &, SparseMatrix<double, Eigen::RowMajor> &, double);



int CliPPCPP(int No_mutation, int* c_r, int *c_n, int *c_minor, int *c_total, double ploidy, double* Lambda_list, int Lambda_num, double alpha, double rho, double gamma, int Run_limit, double precision, int control_large, int least_mut, double post_th, double least_diff, double* c_coef_1d, double* c_wcut_1d, double purity, char* preliminary){

  std::string preliminary_folder(preliminary);
  //std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
  
#ifdef _OPENMP 
    int num_threads = 1;
#pragma omp parallel 
    {
#pragma omp single
	num_threads = omp_get_num_threads();
    }

    omp_set_num_threads(num_threads);
    Eigen::setNbThreads(num_threads);
#endif
    
    int p_i, p_j, p_count, p_starting; 

    //double Lambda = Lambda_list[6];
    
    MatrixXd r(No_mutation, 1);
    MatrixXd n(No_mutation, 1);
    MatrixXd minor_(No_mutation, 1);
    MatrixXd total(No_mutation, 1);
    MatrixXd theta_hat(No_mutation, 1);
    MatrixXd phi_hat(No_mutation, 1);
    
    MatrixXi row1(No_mutation * (No_mutation - 1) / 2, 1);
    MatrixXi row2(No_mutation * (No_mutation - 1) / 2, 1);
    MatrixXi ids(No_mutation * (No_mutation - 1) / 2, 2);
    MatrixXd wcut_1d(No_mutation * 2, 1);
    MatrixXd coef_1d(No_mutation * 6, 1);
    
    for(p_i = 0; p_i < No_mutation; p_i++){
	r(p_i, 0) = double(c_r[p_i]);
	n(p_i, 0) = double(c_n[p_i]);
	minor_(p_i, 0) = double(c_minor[p_i]);
	total(p_i, 0) = double(c_total[p_i]);
    }

    theta_hat = r.array() / n.array();
    phi_hat = theta_hat.array() * ((ploidy - purity * ploidy + purity * total.array()) / minor_.array());
    
    double phi_hat_max = phi_hat.maxCoeff();
    double scale_parameter = max(1.0, phi_hat_max);
    
    for(p_i = 0; p_i < No_mutation * 2; p_i ++){
	wcut_1d(p_i, 0) = c_wcut_1d[p_i];
    }

    for(p_i = 0; p_i < No_mutation * 6; p_i ++){
	coef_1d(p_i, 0) = c_coef_1d[p_i];
    }

    p_count = 0;
    for(p_i = 0; p_i < No_mutation; p_i++){
	for(p_j = 0; p_j < No_mutation; p_j++){
	    if(p_j > p_i){
		ids(p_count, 0) = p_i;
		ids(p_count, 1) = p_j;
		p_count += 1;
	    }
	}
    }
    
    p_starting = 0;
    for(p_i = 0; p_i < No_mutation - 1; p_i++){
	for(p_j = 0; p_j < No_mutation - p_i - 1; p_j++){
	    row1(p_starting + p_j, 0) = p_i;
	    row2(p_starting + p_j, 0) = p_i + 1 + p_j;
	}
	p_starting = p_starting + No_mutation - p_i - 1;
    }
    
    typedef Eigen::Triplet<double> T;
    std::vector<T> triplet_list;
    triplet_list.reserve(No_mutation * (No_mutation - 1));

    for(p_i = 0; p_i < No_mutation * (No_mutation - 1) / 2; p_i++){
	triplet_list.push_back(T(row1(p_i, 0), p_i, 1.0));
    }
    
    for(p_i = 0; p_i < No_mutation * (No_mutation - 1) / 2; p_i++){
	triplet_list.push_back(T(row2(p_i, 0), p_i, -1.0));
    }
    SparseMatrix<double, Eigen::RowMajor> DELTA(No_mutation, No_mutation * (No_mutation - 1) / 2);
    DELTA.setFromTriplets(triplet_list.begin(), triplet_list.end());

#ifdef _OPENMP 
#pragma omp parallel for
#endif
    for(int lambda_index = 0; lambda_index < Lambda_num; lambda_index ++){
	double Lambda = Lambda_list[lambda_index];
	//if(Lambda != 0.01) continue;
    
	CliPPIndividual(No_mutation, r, n, minor_, total, ploidy, Lambda, alpha, rho, gamma, Run_limit, precision, control_large, least_mut, post_th, least_diff, coef_1d, wcut_1d, purity, preliminary_folder, theta_hat, phi_hat, ids, DELTA, scale_parameter);
    }

    //std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    //std::cout << "Time cost = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;

    //ofstream time_cost_file;
    
    //std::string time_file_path = preliminary_folder + "/kernel_timing.txt";
    //time_cost_file.open(time_file_path);

    //if(!time_cost_file.is_open()){
    //std::cerr << "Cannot open file " << time_file_path  << std::endl;
    //	return 1;
    //} 
    //time_cost_file << "Time cost = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;;
    
    //time_cost_file.close();
    
    return 0;
}

int  CliPPIndividual(int No_mutation, MatrixXd &r, MatrixXd &n, MatrixXd &minor_, MatrixXd &total, double ploidy, double Lambda, double alpha, double rho, double gamma, int Run_limit, double precision, int control_large, int least_mut, double post_th, double least_diff, MatrixXd &coef_1d, MatrixXd &wcut_1d, double purity, std::string &preliminary_folder, MatrixXd &theta_hat, MatrixXd &phi_hat, MatrixXi &ids, SparseMatrix<double, Eigen::RowMajor> &DELTA, double scale_parameter)
{
    
    int i, j, k, count;
    
    MatrixXd theta(No_mutation, 1);
    MatrixXd w_new(No_mutation, 1);
    MatrixXd w_old(No_mutation, 1);
    
    MatrixXd eta_new(No_mutation * (No_mutation - 1)/2, 1);
    MatrixXd eta_old(No_mutation * (No_mutation - 1)/2, 1);

    MatrixXd tau_new = MatrixXd::Ones(No_mutation * (No_mutation - 1)/2, 1);
    MatrixXd tau_old(No_mutation, 1);
    
    double temp, temp1, temp2;
#ifdef _OPENMP 
#pragma omp parallel for private(temp, temp1, temp2)
#endif
    for(i = 0; i < No_mutation; i++){
	temp = phi_hat(i, 0) / scale_parameter;
	temp1 = expit((double)control_large);
	if(temp > temp1) temp = temp1;
	temp2 = expit( - (double)control_large);
	if(temp < temp2) temp = temp2;
	
	temp = logit(temp);
	if(temp > control_large) temp = control_large;
	if(temp < -control_large) temp = -control_large;
	w_new(i, 0) = temp;
    }

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(i = 0; i < No_mutation * (No_mutation -1)/2; i++){
	eta_new(i, 0) = w_new(ids(i, 0), 0) - w_new(ids(i, 1), 0);    
    }
    
    MatrixXd A(No_mutation, 1);
    MatrixXd B(No_mutation, 1);
    MatrixXd linear(No_mutation, 1);
    MatrixXd Minv(No_mutation, 1);
    MatrixXd Minv_diag(No_mutation, No_mutation);
    MatrixXd Minv_outer(No_mutation, No_mutation);
    MatrixXd delt(No_mutation * (No_mutation - 1) /2, 1);
    MatrixXd inverted(No_mutation, No_mutation);
    
    double tag1, tag2, tag3, tag4, max_val = -10000.0;
    double residual = 100.0;
    double trace_g = 0.0;
    
    k = 0;
    std::vector<int> problematic_snvs;
    
    while(residual > precision && k < Run_limit){
	k += 1;
	std::cout << "\rLambda: " << Lambda << "\titeration: " << k << "\tresidual: " << residual << std::flush;
	//std::cout << "Lambda: " << Lambda << "iteration: " << k << "residual: " << residual << std::endl;
	
	w_old = w_new;
	tau_old = tau_new;
	eta_old = eta_new;
#ifdef _OPENMP
#pragma omp parallel for private(tag1, tag2, tag3, tag4)
#endif
	for(i = 0; i < No_mutation; i++){
	    theta(i, 0) = exp(w_old(i, 0)) * minor_(i, 0) / ( 2 + exp(w_old(i, 0)) * total(i, 0));
	    if(w_old(i, 0) <= wcut_1d(i * 2 + 0, 0)) {
		tag1 = 1.0; tag3 = 0.0;
		
	    }else{
		tag1 = 0.0;
		tag3 = 1.0;
	    }

	    if(w_old(i, 0) >= wcut_1d(i * 2 + 1, 0)) {
		tag2 = 1.0; tag4 = 0.0;
	    }else{
		tag2 = 0.0; tag4 = 1.0;
	    }

	    if(theta(i, 0) >= 1.0) {
	      theta(i, 0) = 0.99;
	      #pragma omp critical
	      {
		problematic_snvs.push_back(i);
	      }
	    }
	    
	    
	    A(i, 0) = sqrt(n(i, 0)) * (tag1 * coef_1d(i * 6 + 1, 0) + tag2 * coef_1d(i * 6 + 5, 0) + tag3 * tag4 * coef_1d(i * 6 + 3, 0) - theta_hat(i, 0)) / (sqrt(theta(i, 0) * (1 - theta(i, 0))));
	    B(i, 0) = sqrt(n(i, 0)) * (tag1 * coef_1d(i * 6 + 0, 0) + tag2 * coef_1d(i * 6 + 4, 0) + tag3 * tag4 * coef_1d(i * 6 + 2, 0)) / (sqrt(theta(i, 0) * (1 - theta(i, 0))));
	}

	//std::cout << "problematic_snvs size:\t" << problematic_snvs.size() << std::endl;

	
        linear = DELTA * (alpha * eta_old + tau_new) - B.cwiseProduct(A);

	Minv = 1.0 / ((B.cwiseProduct(B)).array() + double(No_mutation) * alpha);
	Minv_diag = Minv.asDiagonal();

	trace_g = -alpha * Minv.sum();
	if(isnan(trace_g)){
	    std::cout << "Lambda: " << Lambda << "\titeration: " << k << "\tNan" << std::endl;
	    return -1;
	}
	
	
	Minv_outer = Minv * (Minv.transpose());

	inverted = Minv_diag.array() - 1.0 / (1.0 + trace_g) * (-alpha) * Minv_outer.array();
	w_new = inverted * linear;

#ifdef _OPENMP
#pragma omp parallel for
#endif
	for(i = 0; i < No_mutation; i++){
	    if(w_new(i, 0) > control_large) w_new(i, 0) = control_large;
	    if(w_new(i, 0) < -control_large) w_new(i, 0) = -control_large;
	}
	
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for(i = 0; i < No_mutation * (No_mutation - 1) /2; i++){
	    delt(i, 0) = w_new(ids(i, 0), 0) - w_new(ids(i, 1), 0) - 1.0/alpha * tau_old(i, 0);
	}

#ifdef _OPENMP
#pragma omp parallel for private(tag1, tag2, tag3, tag4, temp)
#endif
	for(i = 0; i < No_mutation * (No_mutation - 1) / 2; i++){

	    temp = delt(i, 0);
	    if(fabs(temp) > gamma * Lambda) {
		tag1 = 1.0; tag3 = 0.0;
	    }else{
		tag1 = 0.0; tag3 = 1.0;
	    }

	    if (fabs(temp) < (Lambda + Lambda / alpha)) {
		tag2 = 1.0; tag4 = 0.0;
	    }else{
		tag2 = 0.0; tag4 = 1.0;
	    }
	    eta_new(i, 0) = temp * tag1 + ST(temp, Lambda / alpha) * tag2 + ST(temp, gamma * Lambda / ((gamma - 1.0) * alpha)) / (1.0 - 1.0 / ((gamma - 1.0) * alpha)) * tag3 *tag4;

	    tau_new(i, 0) = tau_old(i, 0) - alpha * (w_new(ids(i, 0), 0) - w_new(ids(i, 1), 0) - eta_new(i, 0));
	}
	alpha = alpha * rho;
	max_val = -100000.0;

#ifdef _OPENMP	
#pragma omp parallel for private(temp) reduction(max:max_val)
#endif
	for(i = 0; i < No_mutation * (No_mutation - 1) / 2; i++){
	    temp = w_new(ids(i, 0), 0) - w_new(ids(i, 1), 0) - eta_new(i, 0);
	    if (temp > max_val){
		max_val = temp;
	    }
	}

	residual = max_val;
	
	//std::cout << residual << std::endl;
    }

    sort( problematic_snvs.begin(), problematic_snvs.end() );
    problematic_snvs.erase( unique( problematic_snvs.begin(), problematic_snvs.end() ), problematic_snvs.end() );
    
    std::cout << std::endl;
    
    MatrixXd diff(No_mutation, No_mutation);

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(i = 0; i < No_mutation * (No_mutation - 1) / 2; i++){
	if(fabs(eta_new(i, 0)) <= post_th) eta_new(i, 0) = 0.0;
    } 
#ifdef _OPENMP
#pragma omp parallel for
#endif
    
    for(i = 0; i < No_mutation * (No_mutation - 1) / 2; i++){
	diff(ids(i, 0), ids(i, 1)) = eta_new(i, 0);
    }
    
    std::vector<int> class_label(No_mutation);
    std::fill(class_label.begin(), class_label.end(), -1);
    class_label[0] = 0;
    
    std::vector<int> group_size;
    group_size.push_back(1);
    
    int labl = 1;

    for(i = 0; i < No_mutation; i++){
	for(j = 0; j < i; j++){
	    if(diff(j, i) == 0.0){
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
    MatrixXd tmp_diff(No_mutation, 1);
    
    if(double(temp_size) < least_mut) refine = 1;
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
		    tmp_diff(jj, 0) = fabs(diff(jj, tmp_col[ii]));
		}
		tmp_diff(tmp_col[ii], 0) = 100.0;

		for(int jj = tmp_col[ii] + 1; jj < No_mutation; jj++){
		    tmp_diff(jj, 0) = fabs(diff(tmp_col[ii], jj));
		}

		for(unsigned int jj = 0; jj < tmp_col.size(); jj++){
		    tmp_diff(tmp_col[jj], 0) = tmp_diff(tmp_col[jj], 0) + 100.0;
		}

		for(unsigned int jj = 0; jj < tmp_col.size(); jj++){
		    diff(jj, tmp_col[ii]) = tmp_diff(jj, 0);
		}

		for(int jj = tmp_col[ii] + 1; jj < No_mutation; jj++){
		    diff(tmp_col[ii], jj) = tmp_diff(jj, 0);
		}
	    }else {
		if(tmp_col[ii] == 0){
		    tmp_diff(0, 0) = 100.0;
		    for(int jj = 1; jj < No_mutation; jj++){
			tmp_diff(jj, 0) = diff(0, jj);
		    }
		    for(unsigned int jj = 0; jj < tmp_col.size(); jj++){
			tmp_diff(tmp_col[jj], 0) = tmp_diff(tmp_col[jj], 0) + 100.0;
		    }
		    for(int jj = 1; jj < No_mutation; jj++){
			diff(0, jj) = tmp_diff(jj, 0);
		    }
		}else{
		    for(int jj = 0; jj < No_mutation - 1; jj++){
			tmp_diff(jj, 0) = diff(jj, No_mutation - 1);
		    }
		    tmp_diff(No_mutation - 1, 0) = 100.0;
		    for(unsigned int jj = 0; jj < tmp_col.size(); jj++){
			tmp_diff(tmp_col[jj], 0) = tmp_diff(tmp_col[jj], 0) + 100.0;
		    }
		    for(int jj = 0; jj < No_mutation - 1; jj++){
			diff(jj, No_mutation - 1) = tmp_diff(jj, 0);
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
	    temp += phi_hat(ind[i], 0) * n(ind[i], 0);
	    temp1 += n(ind[i], 0);
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
		    temp += phi_hat(ind[i], 0) * n(ind[i], 0);
		    temp1 += n(ind[i], 0);
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

    // write to file
    
    ofstream phi_file, label_file, summary_file;

    std::string lambda_str = std::to_string(Lambda);
    while (lambda_str[lambda_str.size() - 1] == '0' || lambda_str[lambda_str.size() - 1] == '.')
	lambda_str.resize(lambda_str.size() - 1);
    
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


extern "C" {
    void CliPP(int No_mutation, int* r, int *n, int* minor, int* total, double ploidy,
	      double* Lambda_list, int Lambda_num, double alpha, double rho, double gamma, int Run_limit, double precision,
	      int control_large, int least_mut, double post_th, double least_diff,
	      double* coef_1d, double* wcut_1d, double purity, char* preliminary){

	CliPPCPP(No_mutation, r, n, minor, total, ploidy, Lambda_list, Lambda_num, alpha, rho, gamma, Run_limit, precision, control_large, least_mut, post_th, least_diff, coef_1d, wcut_1d, purity, preliminary);
    }
}
