/***********************************************************/
/**	\file
	\brief		Array Check Functions
	\author		Zhiqi Li, based on yzLib of Dr.Yizhong Zhang
	\date		9/28/2012
*/
/***********************************************************/
#ifndef __PHYSICS_NEIGHBOR_SEARCH_H__
#define __PHYSICS_NEIGHBOR_SEARCH_H__

#include <iostream>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <cstdlib>
#include <string.h>
#include<zqBasicUtils/utils_hash.h>
#include<zqBasicUtils/utils_common.h>
#include<ResearchP_config.h>
/*
#ifdef  RESEARCHP_ENABLE_CUDA
#include <zqBasicUtils/utils_cuda.h>
#endif //  RESEARCHM_ENABLE_CUDA
*/
#include<zqBasicUtils/utils_hash.h>
#include<zqBasicUtils/utils_aux_func.h>
#include<zqBasicMath/math_vector.h>
#include<zqPhysics/physics_grid.h>

namespace zq{ namespace physics{
	
	/// This version of NbsSearcher can not be used in kernel function,
	/// because Array could not be used in kernel function
	template<int d>
	class NbsSearcher {
		Typedef_VectorDDi(d);
	public:
		NbsSearcher(){}

		void update_points(Array<VecD>& points) {
			Assert(false,"update_points of NbsSearcher must be override");
			exit(0);
		}

		void find_nbs(const VecD& pos, Array<int>& ans) {
			Assert(false, "update_points of NbsSearcher must be override");
			exit(0);
		}

		template<class Cond>
		void find_nbs(const VecD& pos, Array<int>& ans, const Cond& phi) {
			Assert(false, "update_points of NbsSearcher must be override");
			exit(0);
		}

		int find_nearest(const VecD& pos, Array<int>& ans) {
			Assert(false, "find_nearest of NbsSearcher must be override");
			exit(0);
		}

		template<class Cond>
		int find_nearest(const VecD& pos, Array<int>& ans, const Cond& phi) {
			Assert(false, "update_points of NbsSearcher must be override");
			exit(0);
		}

	};

	/// Spatial Hasing
	template<int d>
	class SpatialHashing:NbsSearcher<d> {
		Typedef_VectorDDi(d);
	public:
		Grid<d> grid;
		std::unordered_map<int, Array<int>> hash_table;
		Array<VecD> points;
	public:
		SpatialHashing() { hash_table.clear(); }
		SpatialHashing(real dx, const Array<VecD>& points){
			VecD max_pos = getMaxPosition(points);
			VecD min_pos = getMinPosition(points);
			VecD d = (max_pos - min_pos);
			max_pos = d * 2 + min_pos;
			min_pos = min_pos - d;
			SpatialHashing(dx, max_pos, min_pos);
		}
		SpatialHashing(real dx, const VecD& max_pos, const VecD& min_pos) {
			hash_table.clear();
			grid = Grid<d>(dx, max_pos, min_pos);
		}

		void updatePoints(Array<VecD>& points) {
			this->points = points;
			hash_table.clear();
			for (int i = 0; i < points.size(); i++) {
				VecDi cell = grid.getCell(points[i]);
				int cellID = grid.getCellId(cell);
				if (hash_table.find(cellID) == hash_table.end()) {
					hash_table.insert(std::pair<int, Array<int>>(cellID, Array<int>()));
				}
				hash_table[cellID].push_back(i);
			}
		}

		void findNbs(const VecD& pos, Array<int>& ans) {
			VecDi cell = grid.getCell(pos);
			int nb_cell_num = d == 1 ? 3 : (d == 2 ? 9 : 27);
			for (int i = 0; i < nb_cell_num; i++) {
				VecDi adjCell;
				if (grid.getNextCell(cell, i, adjCell)) {
					int cellId = grid.getCellId(adjCell);
					if (hash_table.find(cellId) != hash_table.end()) {
						for (int j = 0; j < hash_table[cellId].size(); j++) {
							int point_index = hash_table[cellId][j];
							if ((points[point_index] - pos).Length() <= grid.dx && f(point_index)) ans.push_back(point_index);
						}
					}
				}
			}
		}

		template<class Cond>
		void findNbs(const VecD& pos, Array<int>& ans, const Cond& phi) {
			VecDi cell = grid.getCell(pos);
			int nb_cell_num = d == 1 ? 3 : (d == 2 ? 9 : 27);
			for (int i = 0; i < nb_cell_num; i++) {
				VecDi adjCell;
				if (grid.getNextCell(cell, i, adjCell)) {
					int cellId = grid.getCellId(adjCell);
					if (hash_table.find(cellId) != hash_table.end()) {
						for (int j = 0; j < hash_table[cellId].size(); j++) {
							int point_index = hash_table[cellId][j];
							if ((points[point_index] - pos).Length() <= grid.dx && phi(point_index)) ans.push_back(point_index);
						}
					}
				}
			}
		}

		int findNearest(const VecD& pos, Array<int>& ans) {
			findNbs(pos, ans);
			return findNearestIn(pos, ans);
		}

		template<class Cond>
		int findNearest(const VecD& pos, Array<int>& ans, const Cond& phi) {
			findNbs(pos, ans,phi);
			return findNearestIn(pos, ans);
		}
	protected:
		real findMinCoord(Array<VecD> ans, int idx=0) {
			Assert(ans.size() != 0, "findMinCoord cannot takte zero points");
			real min_coord = ans[0][idx];
			for (int i = 1; i < ans.size(); i++) {
				if (ans[i][idx] < min_coord) min_coord = ans[i][idx];
			}
			return min_coord;
		}
		real findMaxCoord(Array<VecD> ans, int idx = 0) {
			Assert(ans.size() != 0, "findMinCoord cannot takte zero points");
			real max_coord = ans[0][idx];
			for (int i = 1; i < ans.size(); i++) {
				if (ans[i][idx] > max_coord) max_coord = ans[i][idx];
			}
			return max_coord;
		}
		VecD getMinPosition(const Array<VecD>& points) {
			/// return the lower bound for bounding box
			if constexpr (d == 1) {
				real least1 = findMinCoord(points, 0);
				return VecD(least1);
			}
			else if constexpr (d == 2) {
				real least1 = findMinCoord(points, 0);
				real least2 = findMinCoord(points, 1);
				return VecD(least1, least2);
			}
			else if constexpr (d == 3) {
				real least1 = findMinCoord(points, 0);
				real least2 = findMinCoord(points, 1);
				real least3 = findMinCoord(points, 2);
				return VecD(least1, least2, least3);
			}
		}
		VecD getMaxPosition() {
			///return the upper bound for bounding box
			if constexpr (d == 1) {
				real largest1 = findMaxCoord(points, 0);
				return VecD(largest1);
			}
			else if constexpr (d == 2) {
				real largest1 = findMaxCoord(points, 0);
				real largest2 = findMaxCoord(points, 1);
				return VecD(largest1, largest2);
			}
			else if constexpr (d == 3) {
				real largest1 = findMaxCoord(points, 0);
				real largest2 = findMaxCoord(points, 1);
				real largest3 = findMaxCoord(points, 2);
				return VecD(largest1, largest2, largest3);
			}
		}
		int findNearestIn(const VecD& pos, Array<int>& point_index) {
			int nearst_idx = -1;
			for (int i = 0; i < point_index.size(); i++) {
				int idx = point_index[i];
				if (nearst_idx<0 || (points[nearst_idx] - pos).Length()>(points[idx] - pos).Length()) {
					nearest_idx = idx;
				}
			}
			return idx;
		}
	};

	/// Kd tree
	/// Todo
	

}}	
#endif	