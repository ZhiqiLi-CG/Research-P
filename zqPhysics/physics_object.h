/***********************************************************/
/**	\file       Physics Object
	\brief		Physics Object
	\author		Zhiqi Li
	\date		1/9/2022
*/
/***********************************************************/
#ifndef __PHYSICS_OBJECT_H__
#define __PHYSICS_OBJECT_H__

#include <iostream>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <cstdlib>
#include <string.h>
#include<zqBasicMath/math_matrix.h>
#include<zqBasicMath/math_vector.h>
#include<zqPhysics/physics_boundary.h>
#include<zqPhysics/physics_neighbor_search.h>
#include<ResearchP_config.h>
#include<zqPhysics/physics_sph_kernel.h>
#ifdef  RESEARCHP_ENABLE_CUDA
#include <zqBasicUtils/utils_cuda.h>
#endif //  RESEARCHP_ENABLE_CUDA
#include <zqBasicUtils/utils_array_func.h>
#include <zqBasicUtils/utils_json.h>
namespace zq{ namespace physics{
	template<int side,int d>
	void printVecArray(const Array<Vec<real,d>,side>& var,const char* promp,bool force) {
		if (force) {
			Array<Vec<real, d>> new_var = var;
			printf("%s",promp); 
			printf("--------------------------------------------\n");
			for (int i = 0; i < new_var.size(); i++) {
				printf("(");
				for (int j = 0; j < d; j++) {
					printf("%f,", new_var[i][j]);
				}
				printf(")\t");
			}
			printf("--------------------------------------------\n");
		}
	}
	template<int side, int d>
	void printMatArray(const Array<Mat<real, d>,side>& var, const char* promp, bool force) {
		if (force) {
			Array<Mat<real, d>> new_var = var;
			printf("%s", promp);
			printf("--------------------------------------------\n");
			for (int i = 0; i < new_var.size(); i++) {
				printf("(");
				for (int j = 0; j < d; j++) {
					for (int k = 0; k < d; k++) {
						printf("%f,", new_var[i](j, k));
					}
				}
				printf(")\t");
			}
			printf("--------------------------------------------\n");
		}
	}
	template<int side>
	void printRealArray(const Array<real,side>& var, const char* promp, bool force) {
		if (force) {
			Array<real> new_var = var;
			printf("%s", promp);
			printf("--------------------------------------------\n");
			for (int i = 0; i < new_var.size(); i++) {
				printf("(");
				printf("%f,", new_var[i]);
				printf(")\t");
			}
			printf("--------------------------------------------\n");
		}
	}
	template<int side>
	void printArrayArray(const Array<Array<int, side>>& var, const char* promp, bool force) {
		if (force) {
			Array<Array<int>> new_var = var;
			printf("%s", promp);
			printf("--------------------------------------------\n");
			for (int i = 0; i < new_var.size(); i++) {
				printf("(");
				for (int j = 0; j < new_var[i].size(); j++) {
					printf("%f,", new_var[i][j]);
				}
				printf(")\t");
			}
			printf("--------------------------------------------\n");
		}
	}
	/// Process Objects comprise six parts:
	///		1. Objects, which are described by attribute array and const 
	///		2. Update Componenets
	///		3. Update Main Function
	///		4. Bounday Modification
	///		5. Input Methods(Including init methods)
	///		6. Output Methods
	/// Note that each process objects describe the objects that take the same calculation process
	/// For example, 
	///		Obj1: codim_sph, Obj2: euler,....
	template<int d, int side = HOST>
	class ProcessObjects {
		Typedef_VectorDDi(d);
		Typedef_MatrixDDi(d);
	public:
		/// 1. array and const
		std::unordered_map<std::string, real*> key_var;
		std::unordered_map<std::string, void*> key_array;
	public:
		/// 2. update components, which is not set by base
	public:
		/// 3. update main 
		void update(real dt){}
	public:
		/// 4. Bounadry Control
		BoundaryControl* bc_arr[3];
		void set_boundary_control(
			BoundaryControl* bc1,
			BoundaryControl* bc2,
			BoundaryControl* bc3
		) {
			bc_arr[0] = bc1;
			bc_arr[1] = bc2;
			bc_arr[2] = bc3;
		}
	public:
		/// 5. Input Methods(Including init methods)
		void Init_Constant(const std::unordered_map<std::string, real>& key_value) {
			std::unordered_set < std::string> name_set;
			for (auto& iter : key_value) {
				if (key_var.find(iter.first) != key_var.end()) {
					(*(key_var[iter.first])) = iter.second;
					name_set.insert(iter.first);
				}
				else {
					throw "the value you set for Basic Particles is not in the key_var map";
				}
			}
			/// Then check if all the constant have its own value
			if (!Check_Init_Finish(name_set, key_var)) {
				printf("Your init of constant is not complete\n");
				throw "Your init of constant is not complete";
			}
		}
		template<class T>
		bool Check_Init_Finish(const std::unordered_set < std::string>& name_set, const std::unordered_map<std::string, T>& name_map) {
			for (auto& iter : name_map) {
				if (name_set.find(iter.first) == name_set.end()) {
					return false;
				}
			}
			return true;
		}
		template<class T,int this_side>
		void Init_Array(std::string name, const Array<T, this_side>& value) {
			if (key_array.find(name) == key_array.end()) {
				throw "the array you want init is not in Basic particles";
			}
			Array<T,side>* ptr = (Array<T,side>*)(key_array[name]);
			(*((Array<T,side>*)ptr)) = value;
		}
		json settingExample() {
			json j;
			for (auto& iter : key_var) {
				j[(iter.first).c_str()] = 0;
			}
			return j;
		}
	public:
		///	6. Output Methods
	};

}}	
#endif	