/***********************************************************/
/**	\file
	\brief		Physics Grid
	\author		Zhiqi Li
	\date		1/9/2023
*/
/***********************************************************/
#ifndef __PHYSICS_GRID_H__
#define __PHYSICS_GRID_H__

#include <iostream>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <cstdlib>
#include <string.h>
#include<ResearchP_config.h>
#ifdef  RESEARCHP_ENABLE_CUDA
#include <zqBasicUtils/utils_cuda.h>
#endif //  RESEARCHM_ENABLE_CUDA
#include<zqBasicUtils/utils_hash.h>
#include<zqBasicUtils/utils_aux_func.h>
#include<zqBasicMath/math_vector.h>

namespace zq{ namespace physics{

	/// Grid
	template<int d>
	class Grid {
		Typedef_VectorDDi(d);
	public:
		VecD MinPosition;
		VecD MaxPosition;
		real dx;
		int cellNumber;
	public:
		/// Constructor
#ifdef  RESEARCHM_ENABLE_CUDA
		__host__ __device__
#endif
		Grid() {}
		
#ifdef  RESEARCHM_ENABLE_CUDA
		__host__ __device__
#endif
		Grid(real dx, const VecD& MaxPosition, const VecD& MinPosition) :
			MinPosition(MinPosition), MaxPosition(MaxPosition), dx(dx)
		{
			cellNumber = CellNumber();
		}
		
#ifdef  RESEARCHM_ENABLE_CUDA
		__host__ __device__
#endif
		Grid(Grid& grid) {
			this->MinPosition = grid->MinPosition;
			this->MaxPosition = grid->MaxPosition;
			this->dx = grid->dx;
			this->cellNumber = grid->cellNumber;
		}
		
#ifdef  RESEARCHM_ENABLE_CUDA
		__host__ __device__
#endif
		VecDi getCell(VecD pos) {
			return VecDi(((pos - MinPosition) / dx));
		}
		
#ifdef  RESEARCHM_ENABLE_CUDA
		__host__ __device__
#endif
		int getCellId(VecDi cell) {
			VecDi tem = CellNumberDim();
			if constexpr (d == 1) return cell[0];
			else if constexpr (d == 2) return tem[0] * cell[1] + cell[0];
			else if constexpr (d == 3) return tem[1] * tem[0] * cell[2] + tem[0] * cell[1] + cell[0];
		}
		
#ifdef  RESEARCHM_ENABLE_CUDA
		__host__ __device__
#endif
		VecDi getCell(int cellId) {
			VecDi tem = CellNumberDim();
			if constexpr (d == 1) {
				tem[0] = cellId % tem[0];
			}
			else if constexpr (d == 2) {
				tem[1] = cellId / tem[0];
				tem[0] = cellId % tem[0];
			}
			else if constexpr (d == 3) {
				tem[2] = cellId / (tem[1] * tem[0]);
				cellId = cellId % (tem[1] * tem[0]);
				tem[1] = cellId / tem[0];
				tem[0] = cellId % tem[0];
			}
			return tem;
		}

#ifdef  RESEARCHM_ENABLE_CUDA
		__host__ __device__
#endif
		bool getNextCell(VecDi cell, int index, VecDi& ans) {
			if constexpr (d == 1) {
				Assert(index <= 2, "getNextCell out of index:{}", index);
				VecDi add_term_1D[3]{
						VecDi(-1),VecDi(0),VecDi(1),
				};
				ans = add_term_1D[index] + cell;
			}
			if constexpr (d == 2) {
				Assert(index <= 8, "getNextCell out of index:{}", index);
				VecDi add_term_2D[9]{
					VecDi(-1,-1),VecDi(0,-1),VecDi(1,-1),
					VecDi(-1,0),VecDi(0,0),VecDi(1,0),
					VecDi(-1,1),VecDi(0,1),VecDi(1,1)
				};
				ans = add_term_2D[index] + cell;
			}
			else if constexpr (d == 3) {
				Assert(index <= 26, "getNextCell out of index:{}", index);
				VecDi add_term_3D[27]{
					VecDi(-1,-1,-1),VecDi(0,-1,-1),VecDi(1,-1,-1),
					VecDi(-1,0,-1),VecDi(0,0,-1), VecDi(1,0,-1),
					VecDi(-1,1,-1),VecDi(0,1,-1),VecDi(1,1,-1),
					VecDi(-1,-1,0),VecDi(0,-1,0),VecDi(1,-1,0),
					VecDi(-1,0,0),VecDi(0,0,0),VecDi(1,0,0),
					VecDi(-1,1,0),VecDi(0,1,0),VecDi(1,1,0),
					VecDi(-1,-1,1),VecDi(0,-1,1),VecDi(1,-1,1),
					VecDi(-1,0,1),VecDi(0,0,1), VecDi(1,0,1),
					VecDi(-1,1,1),VecDi(0,1,1),VecDi(1,1,1)
				};
				ans = add_term_3D[index] + cell;
			}
			return ValidCell(ans);
		}
		
#ifdef  RESEARCHM_ENABLE_CUDA
		__host__ __device__
#endif
		int CellNumber() {
			VecDi tem = CellNumberDim();
			if constexpr (d == 1) return tem[0];
			else if constexpr (d == 2) return tem[0] * tem[1];
			else if constexpr (d == 3) return tem[0] * tem[1] * tem[2];
		}
		
#ifdef  RESEARCHM_ENABLE_CUDA
		__host__ __device__
#endif
		VecDi CellNumberDim() {
			/// Note: all the size will increase one!
			return VecDi((MaxPosition - MinPosition) / dx) + vec_one<real,d>();
		}
		
#ifdef  RESEARCHM_ENABLE_CUDA
		__host__ __device__
#endif
		bool ValidCell(VecDi cell) {
			VecDi tem = CellNumberDim();
			if constexpr (d == 1) {
				if (cell[0] < 0) return false;
				if (cell[0] >= tem[0]) return false;
			}
			else if constexpr (d == 2) {
				if (cell[0] < 0 || cell[1] < 0) return false;
				if (cell[0] >= tem[0] || cell[1] >= tem[1]) return false;
			}
			else if constexpr (d == 3) {
				if (cell[0] < 0 || cell[1] < 0 || cell[2] < 0) return false;
				if (cell[0] >= tem[0] || cell[1] >= tem[1] || cell[2] >= tem[2]) return false;
			}
			return true;
		}
		
	};

}}	
#endif	