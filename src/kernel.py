'''
CliP the python version, tested with 3.3.1 and 3.4.5. Should work with all 3.X and 2.7.6+, but not tested
Author: Kaixian Yu, Yujie Jiang, Yuxin Tang
email: yujiejiang679@gmail.com
04/02/2021
'''
import numpy as np
import scipy as sp
import scipy.sparse
from scipy.special import logit
from scipy.special import expit

# the Soft Thresholding function
def ST(x, lam):
	val = np.abs(x) - lam
	val = np.sign(x)*(val > 0) * val
	return val

# The CliP function
def CliP(r, n, minor, total, ploidy, Lambda, alpha, rho, gamma, Run_limit, precision,
         control_large, least_mut, post_th, least_diff, coef, wcut, purity):

	No_mutation = len(r)
	NO_MUTATION = len(r)
	# VAF
	theta_hat = r / n
	phi_hat = theta_hat * ((ploidy - purity * ploidy + purity * total) / minor)
	# constrain phi_hat in (0,1)
	scale_parameter = np.max([1, np.max(phi_hat)])
	phi_new = phi_hat / scale_parameter
	# control_large is used to avoid having large result
	phi_new[phi_new > expit(control_large)] = expit(control_large)
	phi_new[phi_new < expit(-control_large)] = expit(-control_large)
	w_new = logit(phi_new)
	w_new[w_new > control_large] = control_large
	w_new[w_new < -control_large] = -control_large
	k = 0  # iterator
	diff = np.subtract.outer(w_new, w_new)
	ids = np.triu_indices(diff.shape[1], 1)
	eta_new = diff[ids]
	tau_new = np.ones((int(No_mutation * (No_mutation - 1) / 2), 1))
	col_id = np.append(np.array(range(int(No_mutation * (No_mutation - 1) / 2))),
					   np.array(range(int(No_mutation * (No_mutation - 1) / 2))))
	row1 = np.zeros(int(No_mutation * (No_mutation - 1) / 2))
	row2 = np.zeros(int(No_mutation * (No_mutation - 1) / 2))
	starting = 0
	for i in range(No_mutation - 1):
		row1[starting:(starting + No_mutation - i - 1)] = i
		row2[starting:(starting + No_mutation - i - 1)] = np.array(range(No_mutation))[(i + 1):]
		starting = starting + No_mutation - i - 1
	row_id = np.append(row1, row2)
	vals = np.append(np.ones(int(No_mutation * (No_mutation - 1) / 2)),
					 -np.ones(int(No_mutation * (No_mutation - 1) / 2)))
	DELTA = sp.sparse.coo_matrix((vals, (row_id, col_id)),
								 shape=(No_mutation, int(No_mutation * (No_mutation - 1) / 2))).tocsr()
	residual = 100

	while residual > precision and k < Run_limit:
		print('\r', k, ',', residual, end="")
		k = k + 1
		w_old = w_new
		tau_old = tau_new
		eta_old = eta_new
		theta = np.exp(w_old) * minor / (2 + np.exp(w_old) * total)

		A = np.sqrt(n) * (
					(w_old <= wcut[:, 0]) * coef[:, 1] + (w_old >= wcut[:, 1]) * coef[:, 5] + (w_old > wcut[:, 0]) * (
						w_old < wcut[:, 1]) * coef[:, 3] - theta_hat) / np.sqrt(theta * (1 - theta))
		B = np.sqrt(n) * (
					(w_old <= wcut[:, 0]) * coef[:, 0] + (w_old >= wcut[:, 1]) * coef[:, 4] + (w_old > wcut[:, 0]) * (
						w_old < wcut[:, 1]) * coef[:, 2]) / np.sqrt(theta * (1 - theta))

		linear = (DELTA * np.matrix((alpha * eta_old + tau_new.T).T)).flatten() - (B * A)

		Minv = 1 / (B ** 2 + No_mutation * alpha)
		Minv_diag = np.diag(Minv)

		trace_g = -alpha * np.sum(Minv)

		Minv_outer = np.outer(Minv,Minv)
		inverted = Minv_diag - (1 / (1 + trace_g) * (-alpha) * Minv_outer)
		w_new    = np.matmul(inverted, linear.T)

		w_new = np.array(w_new).ravel()
		w_new[w_new > control_large] = control_large
		w_new[w_new < -control_large] = -control_large
		diff = np.subtract.outer(w_new, w_new)
		delt = (diff[ids] - 1 / alpha * tau_old.T).ravel()
		eta_new = delt * (np.abs(delt) > gamma * Lambda) + ST(delt, Lambda / alpha) * (
					np.abs(delt) < (Lambda + Lambda / alpha)) + ST(delt, gamma * Lambda / ((gamma - 1) * alpha)) / (
							  1 - 1 / ((gamma - 1) * alpha)) * (np.abs(delt) <= (gamma * Lambda)) * (
							  np.abs(delt) >= (Lambda + Lambda / alpha))
		tau_new = tau_old - np.array([alpha * (diff[ids] - eta_new)]).T
		alpha = alpha * rho
		residual = np.max(diff[ids] - eta_new)
	# assign mutations based on the distance matrix
	eta_new[np.where(np.abs(eta_new) <= post_th)] = 0
	diff[ids] = eta_new
	class_label = -np.ones(No_mutation)
	class_label[0] = 0
	group_size = [1]
	labl = 1

	for i in range(1, No_mutation):
		for j in range(i):
			if diff[j, i] == 0:
				class_label[i] = class_label[j]
				group_size[int(class_label[j])] += 1
				break
		if class_label[i] == -1:
			class_label[i] = labl
			labl += 1
			group_size.append(1)


	# quality control
	tmp_size = np.min(np.array(group_size)[np.array(group_size) > 0])
	tmp_grp = np.where(group_size == tmp_size)
	refine = False
	if tmp_size < least_mut:
		refine = True
	while refine:
		refine = False
		tmp_col = np.where(class_label == tmp_grp[0][0])[0]
		for i in range(len(tmp_col)):
			if tmp_col[i] != 0 and tmp_col[i] != No_mutation - 1:
				tmp_diff = np.abs(np.append(np.append(diff[0:tmp_col[i], tmp_col[i]].T.ravel(), 100),
											diff[tmp_col[i], (tmp_col[i] + 1):No_mutation].ravel()))
				tmp_diff[tmp_col] += 100
				diff[0:tmp_col[i], tmp_col[i]] = tmp_diff[0:tmp_col[i]]
				diff[tmp_col[i], (tmp_col[i] + 1):No_mutation] = tmp_diff[(tmp_col[i] + 1):No_mutation]
			elif tmp_col[i] == 0:
				tmp_diff = np.append(100, diff[0, 1:No_mutation])
				tmp_diff[tmp_col] += 100
				diff[0, 1:No_mutation] = tmp_diff[1:No_mutation]
			else:
				tmp_diff = np.append(diff[0:(No_mutation - 1), No_mutation - 1], 100)
				tmp_diff[tmp_col] += 100
				diff[0:(No_mutation - 1), No_mutation - 1] = tmp_diff[0:(No_mutation - 1)]
			ind = tmp_diff.argmin()
			group_size[class_label.astype(np.int64, copy=False)[tmp_col[i]]] -= 1
			class_label[tmp_col[i]] = class_label[ind]
			group_size[class_label.astype(np.int64, copy=False)[tmp_col[i]]] += 1
		tmp_size = np.min(np.array(group_size)[np.array(group_size) > 0])
		tmp_grp = np.where(group_size == tmp_size)
		refine = False
		if tmp_size < least_mut:
			refine = True
	labels = np.unique(class_label)
	phi_out = np.zeros(len(labels))
	for i in range(len(labels)):
		ind = np.where(class_label == labels[i])[0]
		class_label[ind] = i
		phi_out[i] = np.sum(phi_hat[ind] * n[ind]) / np.sum(n[ind])
	if len(labels) > 1:
		sort_phi = np.sort(phi_out)
		phi_diff = sort_phi[1:] - sort_phi[:-1]
		min_val = phi_diff.min()
		min_ind = phi_diff.argmin()
		while min_val < least_diff:
			combine_ind = np.where(phi_out == sort_phi[min_ind])[0]
			combine_to_ind = np.where(phi_out == sort_phi[min_ind + 1])[0]
			class_label[class_label == combine_ind] = combine_to_ind
			labels = np.unique(class_label)
			phi_out = np.zeros(len(labels))
			for i in range(len(labels)):
				ind = np.where(class_label == labels[i])[0]
				class_label[ind] = i
				phi_out[i] = np.sum(phi_hat[ind] * n[ind]) / np.sum(n[ind])
			if len(labels) == 1:
				break
			else:
				sort_phi = np.sort(phi_out)
				phi_diff = sort_phi[1:] - sort_phi[:-1]
				min_val = phi_diff.min()
				min_ind = phi_diff.argmin()
	phi_res = np.zeros(No_mutation)
	for lab in range(len(phi_out)):
		phi_res[class_label == lab] = phi_out[lab]
	return {'phi': phi_res, 'label': class_label}
