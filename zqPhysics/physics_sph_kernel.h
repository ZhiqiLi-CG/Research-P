/***********************************************************/
/**	\file
	\brief		sph kernel
	\author		Zhiqi Li
	\date		1/12/2023
*/
/***********************************************************/
#ifndef __PHYSICS_SPH_KERNEL_H__
#define __PHYSICS_SPH_KERNEL_H__

#include <iostream>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <cstdlib>
#include <string.h>
#include<ResearchP_config.h>
#ifdef  RESEARCHP_ENABLE_CUDA
#include <zqBasicUtils/utils_cuda.h>
#endif //  RESEARCHP_ENABLE_CUDA
#include<zqBasicUtils/utils_hash.h>
#include<zqBasicUtils/utils_aux_func.h>
#include<zqBasicMath/math_vector.h>
#include<zqBasicMath/math_const.h>
#include<zqBasicMath/math_utils.h>

namespace zq{ namespace physics{
	class UnitSPIKY {
		Typedef_VectorDDi(3);
	public:
		VecD alpha;//0,1,2, -> d=1,2,3

#ifdef  RESEARCHP_ENABLE_CUDA
		__host__ __device__
#endif
		UnitSPIKY() { alpha = VecD(2.0, 10.0 / ZQ_PI, 15.0 / ZQ_PI); }

#ifdef  RESEARCHP_ENABLE_CUDA
		__host__ __device__
#endif
		real Weight(const int d, const real r)const { return r < 1 ? alpha[d - 1] * QuickPower3(1 - r) : 0; }

#ifdef  RESEARCHP_ENABLE_CUDA
		__host__ __device__
#endif
		real Grad(const int d, const real r)const { return r < 1 ? alpha[d - 1] * (-3) * QuickPower2(1 - r) : 0; }
	};
	class UnitCUBIC {
		Typedef_VectorDDi(3);
	public:
		VecD alpha;//0,1,2, -> d=1,2,3

#ifdef  RESEARCHP_ENABLE_CUDA
		__host__ __device__
#endif
		UnitCUBIC() { alpha = VecD(4.0 / 3.0, 40.0 / (7.0 * ZQ_PI), 8.0 / ZQ_PI); }
		
#ifdef  RESEARCHP_ENABLE_CUDA
		__host__ __device__
#endif
		real Weight(const int d, const real r)const {
			if (0 <= r && r < 0.5) return alpha[d - 1] * ((real)6 * QuickPower3(r) - (real)6 * QuickPower2(r) + 1);
			else if (0.5 <= r && r < 1) return alpha[d - 1] * 2 * QuickPower3(1 - r);
			else return 0;
		}

#ifdef  RESEARCHP_ENABLE_CUDA
		__host__ __device__
#endif
		real Grad(const int d, const real r)const {
			if (0 <= r && r < 0.5) return alpha[d - 1] * 6.0 * r * (3.0 * r - 2.0);
			else if (0.5 <= r && r < 1) return alpha[d - 1] * (-6.0) * QuickPower2(1.0 - r);
			else return 0;
		}
	};

	enum SPH_Kernel_Type
	{
		SPIKY_SPH_KERNEL,
		CUBIC_SPH_KERNEL
	};

	class SPHKernel {
		real h;
		real h_pows_inv[5];//3d in maximum, so we must have h_pows[5]
	public:
		
		UnitSPIKY unitSPIKY;
		UnitCUBIC unitCUBIC;
		
#ifdef  RESEARCHP_ENABLE_CUDA
		__host__ __device__
#endif
		SPHKernel(const real _h = 1.0) :h(_h) {
			h_pows_inv[0] = 1;
			for (int i = 1; i < 5; i++) { h_pows_inv[i] = h_pows_inv[i - 1] / h; }
		}

		template<int d> 
#ifdef  RESEARCHP_ENABLE_CUDA
		__host__ __device__
#endif
		real Weight(real r, SPH_Kernel_Type type)const {
			if (type == SPH_Kernel_Type::SPIKY_SPH_KERNEL)
				return unitSPIKY.Weight(d, fabs(r / h)) * h_pows_inv[d];
			else if (type == SPH_Kernel_Type::CUBIC_SPH_KERNEL)
				return unitCUBIC.Weight(d, fabs(r / h)) * h_pows_inv[d];
		}

		template<int d> 
#ifdef  RESEARCHP_ENABLE_CUDA
		__host__ __device__
#endif		
		real Grad_Norm(const real r, SPH_Kernel_Type type) const {
			if (type == SPH_Kernel_Type::SPIKY_SPH_KERNEL)
				return unitSPIKY.Grad(d, fabs(r / h)) * h_pows_inv[d + 1];
			else if (type == SPH_Kernel_Type::CUBIC_SPH_KERNEL)
				return unitCUBIC.Grad(d, fabs(r / h)) * h_pows_inv[d + 1];
		}

		template<int d>
#ifdef  RESEARCHP_ENABLE_CUDA
		__host__ __device__
#endif
		Vec<real, d> Grad(const Vec<real,d>& r, SPH_Kernel_Type type)const {
			real r_len = r.Length();
			if (r_len == 0) return r;
			real grad_coeff = Grad_Norm<d>(r_len,type) / r_len;
			return  grad_coeff * r;
		}

	};

}}	
#endif	