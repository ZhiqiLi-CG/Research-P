/***********************************************************/
/**	\file       Physics Process
	\brief		Physics Process
	\author		Zhiqi Li
	\date		1/9/2022
*/
/***********************************************************/
#ifndef __CODIM1_SPH_OBJECT_H__
#define __PHYSICS_SPH_OBJECT_H__

#include <iostream>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <cstdlib>
#include <string.h>
#include<zqBasicMath/math_matrix.h>
#include<zqBasicMath/math_vector.h>
#include<zqPhysics/physics_boundary.h>
#include<zqPhysics/physics_neighbor_search.h>
#include<zqPhysics/physics_sph_kernel.h>
#include<ResearchP_config.h>
#ifdef  RESEARCHP_ENABLE_CUDA
#include <zqBasicUtils/utils_cuda.h>
#endif //  RESEARCHP_ENABLE_CUDA
#include <zqBasicUtils/utils_array_func.h>
#include<zqPhysics/physics_object.h>
#include<zqPhysics/physics_codim1_sph.h>
namespace zq{ namespace physics{
	template<int d, int side = HOST>
	class SimpleFilm :public ProcessObjects<d, side> {
		Typedef_VectorTTDDi(d);
		Typedef_MatrixTTDDi(d);
	#define TransferToLocal() \
		auto t_dot=this->t_dot;auto t_r = this->t_r;auto v_r = this->v_r;auto rho = this->rho;\
		auto gravity= this->gravity;auto gravity_vec= this->gravity_vec;auto gamma_0 = this->gamma_0;auto gamma_a = this->gamma_a;\
		auto alpha_h = this->alpha_h;auto alpha_k = this->alpha_k;auto alpha_d = this->alpha_d;auto alpha_c = this->alpha_c;\
		auto h_0 = this->h_0;auto p_0 =this->p_0;auto mu = this->mu;auto vol_0 = this->vol_0;\
		auto vol_b = this->vol_b;auto t_kernel=this->t_kernel;\
		auto x_ptr = this->x_ptr;auto v_ptr = this->v_ptr;auto f_ptr = this->f_ptr;auto m_ptr = this->m_ptr;\
		auto h_ptr = this->h_ptr;auto vol_ptr = this->vol_ptr;auto s_ptr = this->s_ptr;auto e_ptr = this->e_ptr;\
		auto g_ptr = this->g_ptr;auto ah_ptr = this->ah_ptr;auto gm_ptr = this->gm_ptr;auto GM_ptr = this->GM_ptr;\
		auto p_ptr = this->p_ptr;auto vo_ptr = this->vo_ptr;auto nbs_ptr_ptr = this->nbs_ptr_ptr;auto nbs_num_ptr = this->nbs_num_ptr;

	public:
		real t_dot = (real).2;										////threshold for the dot product between two normals
		real t_r = (real)0.1;										////local tangential radius, the supporting radius to search for the tangential neighbors, initialized based on _dx
		real v_r = (real)0.1;										////volumetric radius, the supporting radius to search for the volumetric neighbors, initialized based on _dx
		real rho = 1e3;
		real gravity = 9.81;										// the gravitational acceleration
		VecD gravity_vec;										// the gravitational acceleration
		real gamma_0 = 7.275e-2;
		real gamma_a = 8.3144598 * 298.15;							// used in update surface tension coeef
		real alpha_h = 1e1;											// unknown, used in updating pressure
		real alpha_k = 1e1;											// unknown, used in updating pressure
		real alpha_d = 1e1;											// unknown, used in updating pressure
		real alpha_c = 1e1;											// unknown, used in updating concentration
		real h_0 = 6e-7;											// rest thickness of the film 
		real p_0 = 10132.5;											// the standard atmospheric pressure
		real mu = 8.9e-4;											// used when calculate viscosity force
		real vol_0 = 1e1;											// the volume at the beginning
		real vol_b = 1e1;
		SPHKernel* t_kernel;												////tangential SPH kernel
	public:
		/// Attribute
		Array<VecD, side> x, v, f;										///x and v are required to init, f is internal result 
		VecD* x_ptr, * v_ptr, * f_ptr;
		Array<real, side> m, h, vol, s;									///m, h and vol are required to init
		real* m_ptr, * h_ptr, * vol_ptr, * s_ptr;
		Array<MatD, side> e;
		Array<VecD, side> init_n;										///g and init_n are required to init, e is set by init_n
		Array<MatT, side> g;
		MatD* e_ptr;
		VecD* init_n_ptr;
		MatT* g_ptr;
		Array<real, side> ah; // advected height					///init by h
		real* ah_ptr;
		Array<real, side> gm; // surface tension coefficient		///internal result
		real* gm_ptr;
		Array<real, side> GM;  // surfactant concentration			///Need init
		real* GM_ptr;
		Array<real, side> p; // pressure							///internal result
		real* p_ptr;
		Array<real, side> vo; // vorticity							///Need init
		real* vo_ptr;
		Array<Array<int, side>> nbs;
		Array<int*, side> nbs_ptr;
		Array<int, side> nbs_num;
		int** nbs_ptr_ptr;
		int* nbs_num_ptr;
		SpatialHashing<d> nbs_searcher;
	public:
		SimpleFilm() {
			/// Init the var map
			key_var = std::unordered_map<std::string, real*>{
				{"t_dot", &t_dot},
				{"t_r",&t_r},
				{"v_r",&v_r},
				{"gravity",&gravity},
				{"gamma_0",&gamma_0},
				{"gamma_a",&gamma_a},
				{"alpha_h",&alpha_h },
				{"alpha_k", &alpha_k},
				{"alpha_d", &alpha_d},
				{"alpha_c", &alpha_c},
				{"h_0", &h_0},
				{"p_0", &p_0},
				{"mu", &mu},
				{"rho", &rho},
			};
			key_array = std::unordered_map<std::string, void*>{
				{"x",(void*)&x},
				{"v",(void*)&v},
				{"h",(void*)&h},
				{"GM",(void*)&GM},
				{"vo",(void*)&vo},
				{"init_n",(void*)&init_n},
				{"g",(void*)&g},
			};
		}
		void init(
			VecD MaxPosition,
			VecD  MinPosition
		) {
#ifdef  RESEARCHP_ENABLE_CUDA
			SPHKernel* tem_kernel = new SPHKernel(t_r);
			cpyh2d(t_kernel, tem_kernel, sizeof(SPHKernel));
#else
			t_kernel = new SPHKernel(t_r);
#endif
			if constexpr (d == 2) {
				gravity_vec = VecD(0, gravity);
			}
			else if constexpr (d == 3) {
				gravity_vec = VecD(0, 0, gravity);
			}
			this->ah = h;
			Array<MatD> local_e(init_n.size());
			Array<VecD> local_n=this->init_n;
			for (int i = 0; i < local_n.size(); i++) {
				local_e[i] = frameFromNormal<d>(local_n[i]);
			}
			this->e = local_e;
			nbs.resize(x.size());
			updateNb(MaxPosition, MinPosition);
			getPtr();
			updateVol(0);  printRealArray<side>(vol, "Init Vol", false);
			updateM(0);  printRealArray<side>(m, "Init Mass", false);
			updateInitS(0);
			vol_0 = vol_b = calculateVolume<d,side>(x, e, s);
			updateG(0);
		}
		void update(real dt) {
			getPtr();
#ifdef  RESEARCHP_ENABLE_CUDA
			thrust::fill(f.begin(), f.end(), VecD());
#else
			std::fill(f.begin(), f.end(), VecD());
#endif
			//checkUpdateReal(updateH(), h, "test ori_h  original",true);
			updateLowercaseGamma(dt); printRealArray<side>(gm, "test gamma", false);
			updatePressure(dt); printRealArray<side>(p, "test pressure", false);
			updateAdvectedHeight(dt); printRealArray<side>(ah, "test ah", false);
			updateVorticity(dt); printRealArray<side>(vo, "test Vo", false);
			updateConcentration(dt); printRealArray<side>(GM, "test Con", false);
			/*
			updateExternalForce(dt); printVecArray<side,d>(f, "test externa force", false);
			updateVorticityForce(dt); printVecArray<side, d>(f, "test VorticityForce", false);
			updatePressureForce(dt); printVecArray<side, d>(f, "test PressureForce", false);
			updateMarangoniForce(dt); printVecArray<side, d>(f, "test Maran Force", false);
			updateCapillaryForces(dt); printVecArray<side, d>(f, "test CapillaryForces", false);
			updateViscosityForces(dt); printVecArray<side, d>(f, "test ViscosityForces", false);
			*/
			updateVelocity(dt); printVecArray<side, d>(v, "test velocity", true);
			updatePosition(dt); printVecArray<side, d>(x, "test Position", true);
			updateNb(); printArrayArray<side>(nbs, "test nbs", false);

			updateFrame(dt); printMatArray<side, d>(e, "test frame", false);
			updateG(dt); printMatArray<side, d-1>(g, "test g", false);
			updateH(dt); printRealArray<side>(h, "test ori_h", false);
			updateS(dt); printRealArray<side>(s, "test s", false);
			vol_b = calculateVolume<d, side>(x, e, s);
		}
	public:
		void getPtr() {
			/// first, resize,
			v.resize(x.size());
			f.resize(x.size());
			m.resize(x.size());
			h.resize(x.size());
			vol.resize(x.size());
			s.resize(x.size());
			e.resize(x.size());
			g.resize(x.size());
			ah.resize(x.size());
			gm.resize(x.size());
			GM.resize(x.size());
			p.resize(x.size());
			vo.resize(x.size());
			nbs.resize(x.size());
			/// This is because resize will move the ptr
#ifdef  RESEARCHP_ENABLE_CUDA
			x_ptr = thrust::raw_pointer_cast(&x[0]);
			v_ptr = thrust::raw_pointer_cast(&v[0]);
			f_ptr = thrust::raw_pointer_cast(&f[0]);
			m_ptr = thrust::raw_pointer_cast(&m[0]);
			h_ptr = thrust::raw_pointer_cast(&h[0]);
			vol_ptr = thrust::raw_pointer_cast(&vol[0]);
			s_ptr = thrust::raw_pointer_cast(&s[0]);
			e_ptr = thrust::raw_pointer_cast(&e[0]);
			g_ptr = thrust::raw_pointer_cast(&g[0]);
			ah_ptr = thrust::raw_pointer_cast(&ah[0]);
			gm_ptr = thrust::raw_pointer_cast(&gm[0]);
			GM_ptr = thrust::raw_pointer_cast(&GM[0]);
			p_ptr = thrust::raw_pointer_cast(&p[0]);
			vo_ptr = thrust::raw_pointer_cast(&vo[0]);
			nbs_ptr_ptr = thrust::raw_pointer_cast(&nbs_ptr[0]);
			nbs_num_ptr = thrust::raw_pointer_cast(&nbs_num[0]);
#else
			x_ptr = &x[0];
			v_ptr = &v[0];
			f_ptr = &f[0];
			m_ptr = &m[0];
			h_ptr = &h[0];
			vol_ptr = &vol[0];
			s_ptr = &s[0];
			e_ptr = &e[0];
			g_ptr = &g[0];
			ah_ptr = &ah[0];
			gm_ptr = &gm[0];
			GM_ptr = &GM[0];
			p_ptr = &p[0];
			vo_ptr = &vo[0];
			nbs_ptr_ptr = &nbs_ptr[0];
			nbs_num_ptr = &nbs_num[0];
#endif
		}
		void updateNb(VecD MaxPosition, VecD  MinPosition) {
			nbs_searcher = SpatialHashing<d>(t_r, MaxPosition, MinPosition);
			updateNb();
		}
		void updateNb() {
			Array<MatD> e_host = e;
			Array<VecD> x_host = x;
			nbs_searcher.updatePoints(x_host);
			Array<int*> nbs_ptr_local(x_host.size());
			Array<int> nbs_num_local(x_host.size());
			for (int i = 0; i < x_host.size(); i++) {
				Array<int> tem;
				findTangNeighbor(i, e_host, x_host, tem);
				nbs[i] = tem;
#ifdef  RESEARCHP_ENABLE_CUDA
				nbs_ptr_local[i] = thrust::raw_pointer_cast(&nbs[i][0]);
#else
				nbs_ptr_local[i] = &nbs[i][0];
#endif
				nbs_num_local[i] = nbs[i].size();
			}
			nbs_ptr = nbs_ptr_local;
			nbs_num = nbs_num_local;
		}
		void updateLowercaseGamma(real dt) {
			/// why + ?
			TransferToLocal();
			zq::utils::Calc_Each(
				[GM_ptr, gamma_0, gamma_a] __device__ __host__(const int idx)->real {
					return gamma_0 - gamma_a * GM_ptr[idx];
				},
			 gm
			);
		}
		void updatePressure(real dt) {
			TransferToLocal();
			zq::utils::Calc_Each(
				[x_ptr, e_ptr, g_ptr, vol_ptr, h_ptr, nbs_ptr_ptr, nbs_num_ptr, t_kernel, alpha_h, h_0, alpha_k, gm_ptr, alpha_d, v_ptr] __device__ __host__(const int idx)->real {
				real k = calculateCurvature<d>(
					idx,
					x_ptr,
					e_ptr,
					g_ptr,
					vol_ptr,
					h_ptr,
					nbs_ptr_ptr,
					nbs_num_ptr,
					t_kernel
				);
				return alpha_h * (h_ptr[idx] / h_0 - 1) + alpha_k * gm_ptr[idx] * k + alpha_d *
					calculateDiv<d>(
						idx,
						x_ptr,
						v_ptr,
						e_ptr,
						g_ptr,
						vol_ptr,
						h_ptr,
						nbs_ptr_ptr,
						nbs_num_ptr,
						t_kernel
					);
			},
				p
				);
		}
		void updateAdvectedHeight(real dt) {
			TransferToLocal();
			zq::utils::Calc_Each(
				[dt, x_ptr, e_ptr, g_ptr, vol_ptr, v_ptr, h_ptr, ah_ptr, t_kernel, nbs_ptr_ptr, nbs_num_ptr, this] __device__ __host__(const int idx)->real {
				return ah_ptr[idx] - ah_ptr[idx] *
					calculateDiv<d>(
						idx,
						x_ptr,
						v_ptr,
						e_ptr,
						g_ptr,
						vol_ptr,
						h_ptr,
						nbs_ptr_ptr,
						nbs_num_ptr,
						t_kernel
					)* dt;
				},
				ah
				);
		}
		void updateVorticity(real dt) {
			TransferToLocal();
			auto vo_ij = [vo_ptr] __device__ __host__(const int i, const int j)->real {
				return vo_ptr[j] - vo_ptr[i];
			};
			zq::utils::Calc_Each(
				[dt, vo_ij, t_kernel, nbs_num_ptr, nbs_ptr_ptr, x_ptr, e_ptr, g_ptr, vol_ptr, h_ptr, alpha_c, vo_ptr]__device__ __host__(const int idx)->real {
				real new_vol = codim1LapSPH<decltype(vo_ij),d>(
					idx,
					nbs_num_ptr[idx],
					nbs_ptr_ptr[idx],
					x_ptr,
					e_ptr,
					g_ptr,
					vol_ptr,
					h_ptr,
					vo_ij,
					t_kernel
				);
				return alpha_c * dt * new_vol + vo_ptr[idx];
			}
			, vo
				);
		}
		void updateConcentration(real dt) {
			TransferToLocal();
			auto GM_ij = [GM_ptr]
				__device__ __host__
				(const int i, const int j)->real {
				return GM_ptr[j] - GM_ptr[i];
			};
			zq::utils::Calc_Each(
				[GM_ptr, alpha_c, dt, nbs_num_ptr, nbs_ptr_ptr, x_ptr, e_ptr, g_ptr, vol_ptr, h_ptr, GM_ij, t_kernel]
			__device__ __host__
			(const int idx)->real {
				return GM_ptr[idx] + alpha_c * dt * codim1LapSPH<decltype(GM_ij), d>(
					idx,
					nbs_num_ptr[idx],
					nbs_ptr_ptr[idx],
					x_ptr,
					e_ptr,
					g_ptr,
					vol_ptr,
					h_ptr,
					GM_ij,
					t_kernel
				);
			},
			GM
				);
		}
		void updateExternalForce(real dt) {
			TransferToLocal();
			zq::utils::Calc_Each(
				[vol_0, vol_b, p_0, f_ptr, m_ptr, gravity_vec, s_ptr, e_ptr]__device__ __host__(const int idx)->VecD {
				real p_b = (vol_0 / vol_b) * p_0;
				VecD newForce = f_ptr[idx] + m_ptr[idx] * gravity_vec + (p_b - p_0) * s_ptr[idx] * (e_ptr[idx]).col(d - 1) / 1000;
				return newForce;
			}
			, f
				);
		}
		void updateVorticityForce(real dt) {
			TransferToLocal();
			zq::utils::Calc_Each(
				[nbs_num_ptr, nbs_ptr_ptr, x_ptr, e_ptr, f_ptr, vo_ptr] __device__ __host__(const int idx) ->VecD {
				for (int k = 0; k < nbs_num_ptr[idx]; k++) {
					int j = nbs_ptr_ptr[idx][k];
					VecD rij = x_ptr[idx] - x_ptr[j];
					VecD rt_ij = planeVector<d>(projectPlane<d>(rij, e_ptr[idx]), e_ptr[idx]);
					VecD newForce = f_ptr[idx] - cross(rt_ij,vo_ptr[j] * e_ptr[j].col(d - 1));
					return newForce;
		}
	}
			, f
				);
}
		void updatePressureForce(real dt) {
			TransferToLocal();
			auto p_f = [p_ptr] __device__ __host__(const int i) -> real {
				return p_ptr[i];
			};
			zq::utils::Calc_Each(
				[p_f, x_ptr, f_ptr, vol_ptr, h_ptr, nbs_num_ptr, nbs_ptr_ptr, e_ptr, g_ptr, t_kernel] __device__ __host__(const int idx)->VecD {
				VecD newForce = f_ptr[idx] + 2 * vol_ptr[idx] *
					codim1GradSymmSPH<decltype(p_f),d>(
						idx,
						nbs_num_ptr[idx],
						nbs_ptr_ptr[idx],
						x_ptr,
						e_ptr,
						g_ptr,
						vol_ptr,
						h_ptr,
						p_f,
						t_kernel
					);
				return newForce;
			}
			, f
				);
		}
		void updateMarangoniForce(real dt) {
			TransferToLocal();
			auto gm_f = [gm_ptr] __device__ __host__(const int i)->real {
				return gm_ptr[i];
			};
			zq::utils::Calc_Each(
				[f_ptr, vol_ptr, h_ptr, gm_f, t_kernel, nbs_num_ptr, nbs_ptr_ptr, x_ptr, e_ptr, g_ptr] __device__ __host__(const int idx)->VecD {
				return  f_ptr[idx] + (vol_ptr[idx] / h_ptr[idx]) *
					codim1GradDiffSPH<decltype(gm_f),d>(
						idx,
						nbs_num_ptr[idx],
						nbs_ptr_ptr[idx],
						x_ptr,
						e_ptr,
						g_ptr,
						vol_ptr,
						h_ptr,
						gm_f,
						t_kernel
					);
			}
			, f
				);
		}
		void updateCapillaryForces(real dt) {
			TransferToLocal();
			auto xt_f = [x_ptr, e_ptr] __device__ __host__(const int i, const int j)->real {
				real value = -dot(x_ptr[i] - x_ptr[j], e_ptr[i].col(d - 1));
				return value;
			};
			zq::utils::Calc_Each(
				[nbs_num_ptr, nbs_ptr_ptr, x_ptr, e_ptr, g_ptr, vol_ptr, h_ptr, xt_f, t_kernel, f_ptr, gm_ptr] __device__ __host__(const int idx)->VecD {
				real lap = codim1LapSPH<decltype(xt_f),d>(
					idx,
					nbs_num_ptr[idx],
					nbs_ptr_ptr[idx],
					x_ptr,
					e_ptr,
					g_ptr,
					vol_ptr,
					h_ptr,
					xt_f,
					t_kernel
				);
				return  f_ptr[idx] + (vol_ptr[idx] * gm_ptr[idx] / h_ptr[idx]) * e_ptr[idx].col(d - 1) * lap;
			}
			, f
				);
		}
		void updateViscosityForces(real dt) {
			TransferToLocal();
			auto uij_f = [v_ptr, e_ptr] __device__ __host__(const int i, const int j)->VecD {
				VecD u_ij = v_ptr[j] - v_ptr[i];
				VecD returnV = u_ij - dot(u_ij, e_ptr[i].col(d - 1)) * e_ptr[i].col(d - 1);
				return returnV;
			};
			zq::utils::Calc_Each(
				[f_ptr, vol_ptr, mu, nbs_num_ptr, nbs_ptr_ptr, x_ptr, e_ptr, g_ptr,  h_ptr, uij_f, t_kernel] __device__ __host__(const int idx)->VecD {
				return  f_ptr[idx] + (vol_ptr[idx] * mu) *
					codim1LapSPH<decltype(uij_f), d>(
						idx,
						nbs_num_ptr[idx],
						nbs_ptr_ptr[idx],
						x_ptr,
						e_ptr,
						g_ptr,
						vol_ptr,
						h_ptr,
						uij_f,
						t_kernel
					);
			}
			, f
				);
		}
		void updateVelocity(real dt) {
			TransferToLocal();
			zq::utils::Calc_Each(
				[f_ptr, v_ptr, m_ptr, dt] __device__ __host__(const int idx)->VecD {
				return v_ptr[idx] + (f_ptr[idx] / m_ptr[idx]) * dt;
			}
			, v
				);
		}
		void updatePosition(real dt) {
			TransferToLocal();
			zq::utils::Calc_Each(
				[x_ptr, v_ptr, dt] __device__ __host__(const int idx)->VecD {
				return x_ptr[idx] + v_ptr[idx] * dt;
			}
			, x
				);
		}
		void updateFrame(real dt) {
			TransferToLocal();
			zq::utils::Calc_Each(
				[nbs_num_ptr, nbs_ptr_ptr, x_ptr, e_ptr, v_r] __device__ __host__(const int idx)->MatD {
				return framePCA<d>(
					idx,
					nbs_num_ptr[idx],
					nbs_ptr_ptr[idx],
					x_ptr,
					e_ptr,
					v_r
				);
			}
			, e
				);
		}
		void updateG(real dt) {
			TransferToLocal();
			zq::utils::Calc_Each(
				[nbs_num_ptr, nbs_ptr_ptr, x_ptr, e_ptr] __device__ __host__(const int idx)->MatT {
				return calculateTensor<d>(
					idx,
					nbs_num_ptr[idx],
					nbs_ptr_ptr[idx],
					x_ptr,
					e_ptr
				);
			}
			, g
				);
		}
		void updateH(real dt) {
			TransferToLocal();
			zq::utils::Calc_Each(
				[nbs_num_ptr, nbs_ptr_ptr, x_ptr, e_ptr, vol_ptr, t_kernel] __device__ __host__(const int idx)->real {
				return codim1HeightSPH<d>(
					idx,
					nbs_num_ptr[idx],
					nbs_ptr_ptr[idx],
					x_ptr,
					e_ptr,
					vol_ptr,
					t_kernel
				);
			}
			, h
				);
		}
		void updateVol(real dt) {
			TransferToLocal();
			zq::utils::Calc_Each(
				[h_ptr, nbs_num_ptr, nbs_ptr_ptr, x_ptr, e_ptr, t_kernel] __device__ __host__(const int idx)->real {
				return h_ptr[idx] * codim1AreaSPH<d>(
					idx,
					nbs_num_ptr[idx],
					nbs_ptr_ptr[idx],
					x_ptr,
					e_ptr,
					t_kernel
				);
			}
			, vol
				);
		}
		void updateInitS(real dt) {
			TransferToLocal();
			zq::utils::Calc_Each(
				[nbs_num_ptr, nbs_ptr_ptr, x_ptr, e_ptr, t_kernel] __device__ __host__(const int idx)->real {
				return codim1AreaSPH<d>(
					idx,
					nbs_num_ptr[idx],
					nbs_ptr_ptr[idx],
					x_ptr,
					e_ptr,
					t_kernel
				);
			}
			, s
				);
		}
		void updateS(real dt) {
			TransferToLocal();
			zq::utils::Calc_Each(
				[vol_ptr,h_ptr] __device__ __host__(const int idx)->real {
				return vol_ptr[idx] / h_ptr[idx];
			}
			, s
		);
		}
		void updateM(real dt) {
			TransferToLocal();
			zq::utils::Calc_Each(
				[vol_ptr,rho] __device__ __host__(const int idx)->real {
				return vol_ptr[idx] * rho;
			}
			, m
				);
		}
	public:
		void findTangNeighbor(
			int i,
			const Array<MatD>& e_host,
			const Array<VecD>& x_host,
			Array<int>& nbs
		) {
			auto cond = [&](int j) {
				return true;
				return isTangNeighbor(x_host[i], e_host[i], x_host[j], e_host[j], t_r, t_dot);
			};
			nbs_searcher.findNbs(x_host[i], nbs, cond);
		}
		bool isTangNeighbor(
			const VecD& pos0,
			const MatD& e0,
			const VecD& pos,
			const MatD& e,
			real t_r,
			real t_dot
		) {
			////check angle
			VecD n = e0.col(d - 1);
			VecD n_p = e.col(d - 1);
			real angle = dot(n,n_p);
			if (angle < t_dot) return false;	////skip the points with large angles
			////check distance
			VecD u = pos - pos0;
			VecT t = projectPlane<d>(u, e0);
			return t.Length() < t_r;
		}
};
}}	
#endif	