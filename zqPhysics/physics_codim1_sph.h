/***********************************************************/
/**	\file       Physics codim 1 sph operator
	\brief		Physics codim 1 sph operator
	\author		Zhiqi Li
	\date		1/9/2022
*/
/***********************************************************/
#ifndef __PHYSICS_COMDIM1_SPH_H__
#define __PHYSICS_COMDIM1_SPH__

#include <iostream>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <cstdlib>
#include <string.h>
#include<zqBasicMath/math_vector.h>
#include<zqBasicMath/math_matrix.h>
#include<zqPhysics/physics_sph_kernel.h>
#include<ResearchP_config.h>
#ifdef  RESEARCHP_ENABLE_CUDA
#include <zqBasicUtils/utils_cuda.h>
#include <zqBasicUtils/utils_array.h>
#endif //  RESEARCHP_ENABLE_CUDA
#include<zqBasicMath/math_MLS_regression.h>
namespace zq{ namespace physics{
	
	template<int d>
#ifdef  RESEARCHP_ENABLE_CUDA
	__host__ __device__
#endif
	Vec<real, d - 1> projectPlane(
		const Vec<real, d>& u,
		const Mat<real, d>& e
	) {
		Vec<real, d - 1> t;
		for (int i = 0; i < d - 1; i++) {
			t[i] = dot(u, e.col(i));
		}
		return t;
	}

	template<int d>
#ifdef  RESEARCHP_ENABLE_CUDA
	__host__ __device__
#endif
	Vec<real, d> planeVector(
		const Vec<real, d - 1>& u,
		const Mat<real, d>& e
	) {
		Vec<real, d> t;
		for (int i = 0; i < d - 1; i++) {
			t += u[i] * e.col(i);
		}
		return t;
	}

	template<int d>
#ifdef  RESEARCHP_ENABLE_CUDA
	__host__ __device__
#endif
	Mat<real,d > frameFromNormal( Vec<real,d> normal) {
		Mat<real,d> new_e;
		if constexpr (d == 2) {
			new_e.colSet(0, orthogonalVector<real, d>(normal).Normalize());
			new_e.colSet(1, normal.Normalize());
			return new_e;
		}
		else if constexpr (d == 3) {
			new_e.colSet(0, orthogonalVector<real,d>(normal).Normalize());
			new_e.colSet(2, normal.Normalize());
			new_e.colSet(1, cross(new_e.col(0),new_e.col(2)));
			return new_e;
		}
	}
	/// Calculate the volume of a bubble
	template<int d,int side>
	real calculateVolume(
		Array<Vec<real,d>,side> x,
		Array<Mat<real,d>, side> e,
		Array<real,side> s
	) {
		real dd = d;
		Array<real> temp(x.size());
		Array<Vec<real, d>> local_x = x;
		Array<Mat<real, d>> local_e = e;
		Array<real> local_s = s;
		zq::utils::Calc_Each(
			[&] (const int idx)->real {
				return (1.0 / dd) * (local_s[idx]) * dot(local_x[idx], local_e[idx].col(d - 1));
			}
			,
			temp
		);
#ifdef  RESEARCHP_ENABLE_CUDA
		return thrust::reduce(temp.begin(), temp.end());
#else
		real sum = 0;
		for (int i = 0; i < temp.size(); i++) sum += temp[i];
		return sum;
#endif
	}

	template<typename Fint,int d>
#ifdef  RESEARCHP_ENABLE_CUDA
	__host__ __device__
#endif
	real codim1WeightSPH(
		int i,
		int nbs_num,
		int* nbs,
		Vec<real, d>* x,
		Mat<real, d>* e,
		real* vol,
		real* h,
		const Fint& phi,
		SPHKernel* sph_kernel,
		SPH_Kernel_Type kernel_type = SPH_Kernel_Type::SPIKY_SPH_KERNEL
	) {
		decltype(phi(0, 0)) phi_weight;
		if constexpr (std::is_same<decltype(phi(0)), real>::value) {
			phi_weight = 0;
		}
		for (int k = 0; k < nbs_num; k++) {
			int j = nbs[k];
			Vec<real,d> r_ij = x[i] - x[j];
			Vec<real,d-1> coord = projectPlane(r_ij, e[i]);
			phi_weight += vol[j] / h[j] * phi(j) * sph_kernel->Weight<d - 1>(coord.Length(), kernel_type);
		}
		return phi_weight;
	}

	template<int d>
#ifdef  RESEARCHP_ENABLE_CUDA
	__host__ __device__
#endif
	real codim1HeightSPH(
		int i, 
		int nbs_num,
		int* nbs, 
		Vec<real,d>* x, 
		Mat<real, d>* e, 
		real* vol, 
		SPHKernel* sph_kernel,
		SPH_Kernel_Type kernel_type= SPH_Kernel_Type::SPIKY_SPH_KERNEL
	) {
		real phi_real = 0;
		for (int k = 0; k < nbs_num; k++) {
			int j = nbs[k];
			Vec<real,d> r_ij = x[i] - x[j];
			Vec<real, d-1> coord = projectPlane<d>(r_ij, e[i]);
			phi_real += vol[j] * sph_kernel->Weight<d - 1>(coord.Length(), kernel_type);
		}
		return phi_real;
	}

	template<int d>
#ifdef  RESEARCHP_ENABLE_CUDA
	__host__ __device__
#endif
	real codim1AreaSPH(
		int i,
		int nbs_num,
		int* nbs,
		Vec<real, d>* x,
		Mat<real, d>* e,
		SPHKernel* sph_kernel,
		SPH_Kernel_Type kernel_type = SPH_Kernel_Type::SPIKY_SPH_KERNEL
	) {
		real phi_real = 0;
		for (int k = 0; k < nbs_num; k++) {
			int j = nbs[k];
			Vec<real,d> r_ij = x[i] - x[j];
			Vec<real,d-1> coord = projectPlane<d>(r_ij, e[i]);
			phi_real += sph_kernel->Weight<d - 1>(coord.Length(), kernel_type);
		}
		return 1 / phi_real;
	}

	template<typename Fint, int d>
#ifdef  RESEARCHP_ENABLE_CUDA
	__host__ __device__
#endif
	Vec<real, d> codim1GradSymmSPH(
		int i,
		int nbs_num,
		int* nbs,
		Vec<real, d>* x,
		Mat<real, d>* e,
		Mat<real, d-1>* g,
		real* vol,
		real* h,
		const Fint& phi,
		SPHKernel* sph_kernel,
		SPH_Kernel_Type kernel_type = SPH_Kernel_Type::SPIKY_SPH_KERNEL
	) {
		Vec<real,d> grad_phi;
		for (int k = 0; k < nbs_num; k++) {
			int j = nbs[k];
			if (i == j) continue;
			Vec<real,d> r_ij = x[i] - x[j];
			Vec<real,d-1> coord = projectPlane<d>(r_ij, e[i]);
			Vec<real, d> grad=planeVector<d>(g[i].Inverse() * sph_kernel->Grad<d - 1>(coord, kernel_type), e[i]);
			grad_phi += vol[j] * h[i] * grad * (phi(i) / h[i] / h[i] + phi(j) / h[j] / h[j]);
		}
		return grad_phi;
	}

	template<typename Fint, int d>
#ifdef  RESEARCHP_ENABLE_CUDA
	__host__ __device__
#endif
	Vec<real, d> codim1GradDiffSPH(
		int i,
		int nbs_num,
		int* nbs,
		Vec<real, d>* x,
		Mat<real, d>* e,
		Mat<real, d - 1>* g,
		real* vol,
		real* h,
		const Fint& phi,
		SPHKernel* sph_kernel,
		SPH_Kernel_Type kernel_type = SPH_Kernel_Type::SPIKY_SPH_KERNEL
	) {
		Vec<real, d> grad_phi;
		for (int k = 0; k < nbs_num; k++) {
			int j = nbs[k];
			if (i == j) continue;
			Vec<real, d> r_ij = x[i] - x[j];
			Vec<real, d - 1> coord = projectPlane<d>(r_ij, e[i]);
			Vec<real, d> grad = planeVector<d>(g[i].Inverse() * sph_kernel->Grad<d - 1>(coord, kernel_type), e[i]);
			grad_phi += vol[j] / h[j] * grad * (phi(j) - phi(i));
		}
		return grad_phi;
	}

	template<typename Fint, int d>
#ifdef  RESEARCHP_ENABLE_CUDA
	__host__ __device__
#endif
	real codim1DivSPH(
		int i,
		int nbs_num,
		int* nbs,
		Vec<real, d>* x,
		Mat<real, d>* e,
		Mat<real, d - 1>* g,
		real* vol,
		real* h,
		const Fint& phi,
		SPHKernel* sph_kernel,
		SPH_Kernel_Type kernel_type = SPH_Kernel_Type::SPIKY_SPH_KERNEL
	) {
		real div_phi = 0;
		for (int k = 0; k < nbs_num; k++) {
			int j = nbs[k];
			if (i == j) continue;
			Vec<real, d> r_ij = x[i] - x[j];
			Vec<real, d - 1> coord = projectPlane<d>(r_ij, e[i]);
			Vec<real, d> grad = planeVector<d>(g[i].Inverse() * sph_kernel->Grad<d - 1>(coord, kernel_type), e[i]);
			div_phi += vol[j] / h[j] * dot(grad,phi(i, j));
		}
		return div_phi;
	}

	template<typename Fint, int d>
#ifdef  RESEARCHP_ENABLE_CUDA
	__host__ __device__
#endif
	decltype(auto) codim1LapSPH(
		int i,
		int nbs_num,
		int* nbs,
		Vec<real, d>* x,
		Mat<real, d>* e,
		Mat<real, d - 1>* g,
		real* vol,
		real* h,
		const Fint& phi,
		SPHKernel* sph_kernel,
		SPH_Kernel_Type kernel_type = SPH_Kernel_Type::SPIKY_SPH_KERNEL
	) {
		decltype(phi(0, 0)) lap_phi;
		if constexpr (std::is_same<decltype(phi(0, 0)), real>::value) {
			lap_phi = 0;
		}
		for (int k = 0; k < nbs_num; k++) {
			int j = nbs[k];
			if (i == j) continue;
			Vec<real, d> r_ij = x[i] - x[j];
			Vec<real, d - 1> coord = projectPlane<d>(r_ij, e[i]);
			Vec<real, d> grad_kerenel = planeVector<d>((g[i].Inverse()) * sph_kernel->Grad<d - 1>(coord, kernel_type), e[i]);
			real kernel_lap = 2 * grad_kerenel.Length() / coord.Length();
			lap_phi += vol[j] / h[j] * phi(i, j)* kernel_lap;
		}
		return lap_phi;
	}

	template<int d>
#ifdef  RESEARCHP_ENABLE_CUDA
	__host__ __device__
#endif
	real calculateCurvature(
		int idx, 
		Vec<real,d>* x, 
		Mat<real,d>* e, 
		Mat<real,d-1>* g, 
		real* vol,
		real* h, 
		int** nbs_ptr, 
		int* nbs_num, 
		SPHKernel* sph_kernel,
		SPH_Kernel_Type kernel_type = SPH_Kernel_Type::SPIKY_SPH_KERNEL
	) {
		auto phi = [h] __device__ __host__(const int i, const int j)->real {
			return h[j] - h[i];
		};
		return codim1LapSPH<decltype(phi),d>(
			idx,
			nbs_num[idx],
			nbs_ptr[idx],
			x,
			e,
			g,
			vol,
			h,
			phi,
			sph_kernel,
			kernel_type
		);
	}

	template<int d>
#ifdef  RESEARCHP_ENABLE_CUDA
	__host__ __device__
#endif
	real calculateDiv(
		int idx, 
		Vec<real,d>* x, 
		Vec<real, d>* v,
		Mat<real, d>* e, 
		Mat<real, d-1>* g, 
		real* vol,
		real* h, 
		int** nbs_ptr, 
		int* nbs_num, 
		SPHKernel* sph_kernel,
		SPH_Kernel_Type kernel_type = SPH_Kernel_Type::SPIKY_SPH_KERNEL
	) {
		//TOFIX:
		auto phi = [e, v, x]
#ifdef  RESEARCHP_ENABLE_CUDA
			__host__ __device__
#endif
			(const int i, const int j)->Vec<real, d> {
			if constexpr (d == 2) {
				Vec<real, d> r_ij = x[i] - x[j];
				Vec<real, d> ei[2], ej[2];//
				ej[0] = ei[0] = r_ij.Normalize();
				ei[1] = cross(e[i].col(d - 1), ei[0]);
				ej[1] = cross(e[j].col(d - 1), ej[0]);
				Vec<real, d> t;
				for (int k = 0; k < d - 1; k++) {
					t += (dot(ej[k], v[j]) - dot(ei[k], v[j]))* ei[k];
				}
				return t;
			}
			else if constexpr (d == 3) {
				return v[j] - v[i];
			}

		};
		return codim1DivSPH<decltype(phi),d>(
			idx,
			nbs_num[idx],
			nbs_ptr[idx],
			x,
			e,
			g,
			vol,
			h,
			phi,
			sph_kernel,
			kernel_type
		);
	}

#ifdef  RESEARCHP_ENABLE_CUDA
	__host__ __device__
#endif
	real WPCA(const real r, const real v_r)
	{
		if (r < v_r)return (real)1 - pow(r / v_r, 3);
		else return (real)0;
	}

	template<int d>
#ifdef  RESEARCHP_ENABLE_CUDA
	__host__ __device__
#endif
	Mat<real,d> framePCA(
		int idx, 
		int nbs_num, 
		int* nbs, 
		Vec<real,d>* x, 
		Mat<real,d>* e, 
		real v_r
	) {
		Vec<real, d> x_tidle;
		Mat<real, d> Cp;
		real w_total = 0;		
		for (int i = 0; i < nbs_num; i++) {
			real d = (x[nbs[i]] - x[idx]).Length();
			real wij = WPCA(d, v_r);
			x_tidle += WPCA(d, v_r) * x[nbs[i]];
			w_total += wij;
		}
		if (w_total == 0) { printf("Init normal error! invalid w_total"); exit(1); }
		x_tidle /= w_total;
		for (int i = 0; i < nbs_num; i++) {
			real dist = (x[nbs[i]] - x[idx]).Length();
			real wij = WPCA(dist, v_r);
			Cp += wij * mul<real,real,d>(x[nbs[i]] - x_tidle,x[nbs[i]] - x_tidle);
		}Cp /= w_total;
		/// Note that: the normal cannot change dramastically!
		/// the direction for new normal must be close to the old one
		Vec<real,d> eigen_val;
		Mat<real, d> eigen_vec;
		Cp.EigenSymm(eigen_vec,eigen_val);
		Vec<real, d> normal = eigen_vec.col(d - 1);
		if (dot(e[idx].col(d - 1),normal) < 0) normal = -normal;
		return frameFromNormal<d>(normal);
	}

	template<typename F,int d>
#ifdef  RESEARCHP_ENABLE_CUDA
	__host__ __device__
#endif
	void fitLocalMLS(
		int idx, 
		F& phi, 
		int nbs_num, 
		int* nbs, 
		Vec<real,d>* x, 
		Mat<real, d>* e,
		real* coef
	) {
		real* data = new real[d * nbs_num];
		for (int i = 0; i < nbs_num; i++) {
			int nb = nbs[i];
			Vec<real,d-1> tang = projectPlane<d>(x[nb] - x[idx], e[idx]);
			if constexpr (d == 2) {
				data[2* i] = tang[0];
				data[2* i+1] = phi(nb);
			}
			else if constexpr (d == 3) {
				data[3 * i] = tang[0];
				data[3 * i + 1] = tang[1];
				data[3 * i + 2] = phi(nb);
			}
		}
		if constexpr (d == 2) {
			Vec<real, d> zero;
			zq::fitMovingLeastSquare2D<real>(coef, data, nbs_num,&(zero[0]));
		}
		else if constexpr (d == 3) {
			Vec<real, d> zero;
			zq::fitMovingLeastSquare3D<real>(coef, data, nbs_num,&(zero[0]));
		}
		delete[] data;
	}
	template<int d>
#ifdef  RESEARCHP_ENABLE_CUDA
	__host__ __device__
#endif
	void fitLocalShape(
		int idx, 
		int nbs_num, 
		int* nbs, 
		Vec<real,d>* x, 
		Mat<real, d>* e,
		real* coef
	) {
		auto phi = [e, x, idx]
#ifdef  RESEARCHP_ENABLE_CUDA
			__host__ __device__
#endif
		(int j) {
			return dot(e[idx].col(d - 1),x[j] - x[idx]);
		};
		return fitLocalMLS<decltype(phi),d>(
			idx,
			phi,
			nbs_num,
			nbs,
			x,
			e,
			coef
		);
	}
	template<int d>
#ifdef  RESEARCHP_ENABLE_CUDA
	__host__ __device__
#endif
	Mat<real,d-1> calculateTensor(
		int idx, 
		int nbs_num, 
		int* nbs, 
		Vec<real,d>* x, 
		Mat<real, d>* e
	) {
		real coef[3 * (d - 1)];
		fitLocalShape<d>(
			idx,
			nbs_num,
			nbs,
			x,
			e,
			coef
		);
		Mat<real,d-1> metrix_tensor;
		if constexpr (d == 2) {
			metrix_tensor = Mat<real, d - 1>(1 + coef[1] * coef[1]);
		}
		else if constexpr (d == 3) {
			metrix_tensor= Mat < real, d - 1>(1 + coef[1] * coef[1], coef[2] * coef[1], coef[2] * coef[1], 1 + coef[2] * coef[2]);
		}
		return metrix_tensor;
	}
}}	
#endif	