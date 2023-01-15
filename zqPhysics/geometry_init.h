/***********************************************************/
/**	\file       geometry init
	\brief		geometry init
	\author		Zhiqi Li
	\date		1/9/2022
*/
/***********************************************************/
#ifndef __GEOMETRY_INIT_H__
#define __GEOMETRY_INIT_H__

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
#include<ResearchP_config.h>
#include<zqPhysics/physics_sph_kernel.h>
#include <zqBasicUtils/utils_array_func.h>
#include <queue>
namespace zq{ namespace physics{
	/// Only 3D is available, 2D:TODO
	template<int d>
	class InitSphere {
	public:
		Typedef_VectorDDi(3);
		Typedef_MatrixTTDDi(3);
		real radius;
		real dx;
		Array<VecD> points;
		Array<VecD> normals;
		Array<MatT> g;
		Array<Vec<int,3>> triangles;
		VecD maxPosition, minPosition;
		real V, S;
		InitSphere() {}
		InitSphere(real radius, real stop_dx) :radius(radius) {
			Assert(d == 3, "Init Sphere is only valid for d==3");
			dx = splitTriagle(radius, stop_dx, points, triangles);
			V = 4 * radius * radius * ZQ_PI;
			S = V / points.size();
			for (int i = 0; i < points.size(); i++) {
				normals.push_back(points[i].Normalize());
			}
			for (int i = 0; i < points.size(); i++) {
				g.push_back(MatT(1, 0, 0, 1));
			}
			maxPosition = getMaxPosition(points);
			minPosition = getMinPosition(points);
			VecD d = maxPosition - minPosition;
			maxPosition += 10 * d;
			minPosition -= 10 * d;
			printf("finish initSphere\n");
		}
		real ratio(VecD v1, VecD v2, VecD u1, VecD u2) {
			real r1 = (v1 - v2).Length() / (u1 - u2).Length();
			real r2 = (u1 - u2).Length() / (v1 - v2).Length();
			if (r1 >= 1) return r1;
			else return r2;
		}
		bool validTriangle(VecD v1, VecD v2, VecD v3, VecD u1, VecD u2, VecD u3) {
			real threshhold = 4;
			if (ratio(v1, v2, u1, u2) >= threshhold || ratio(v3, v2, u3, u2) >= threshhold || ratio(v1, v3, u1, u3) >= threshhold)
				return false;
			return true;
		}
		static VecD midArc(VecD a, VecD b)
		{
			VecD c(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
			return c / c.Length() * a.Length();
		}
		static int find_pair(std::unordered_map<long long, int>& pairs, int x, int y, int id) {
			long long pair;
			if (x >= y) std::swap(x, y);
			pair = ((((long long)1L) << 32) * ((long long)x)) + y;
			if (pairs.find(pair) != pairs.end()) return pairs[pair];
			else {
				pairs.insert(std::pair<long long, int>(pair, id));
				return id;
			}
		}
		static real splitTriagle(real radius, real stop_dx, Array<VecD>& points, Array<Vec<int,3>>& triangles) {
			std::queue<Vec<int, 3>> triangl_list; points.clear(); triangles.clear();
			points.push_back(VecD(radius, 0, 0));
			points.push_back(VecD(0, radius, 0));
			points.push_back(VecD(0, 0, radius));
			points.push_back(VecD(-radius, 0, 0));
			points.push_back(VecD(0, -radius, 0));
			points.push_back(VecD(0, 0, -radius));
			triangl_list.push(Vec<int,3>(0, 1, 2));
			triangl_list.push(Vec<int, 3>(0, 4, 2));
			triangl_list.push(Vec<int, 3>(3, 4, 2));
			triangl_list.push(Vec<int, 3>(3, 2, 1));

			triangl_list.push(Vec<int, 3>(0, 1, 5));
			triangl_list.push(Vec<int, 3>(0, 4, 5));
			triangl_list.push(Vec<int, 3>(3, 4, 5));
			triangl_list.push(Vec<int, 3>(3, 5, 1));

			real real_dx;
			std::unordered_map<long long, int> pairs;
			while (1) {
				Vec<int, 3> triangle = triangl_list.front();
				VecD a = points[triangle[0]];
				VecD b = points[triangle[1]];
				VecD c = points[triangle[2]];
				if ((a - b).Length() < stop_dx) {
					real_dx = (a - b).Length();
					break;
				}
				VecD ab = midArc(a, b);
				int ab_index = find_pair(pairs, triangle[0], triangle[1], points.size());
				if (ab_index == points.size()) { points.push_back(ab); }
				VecD bc = midArc(b, c);
				int bc_index = find_pair(pairs, triangle[1], triangle[2], points.size());
				if (bc_index == points.size()) { points.push_back(bc); }
				VecD ca = midArc(c, a);
				int ca_index = find_pair(pairs, triangle[2], triangle[0], points.size());
				if (ca_index == points.size()) { points.push_back(ca); }
				triangl_list.push(Vec<int, 3>(triangle[0], ca_index, ab_index));
				triangl_list.push(Vec<int, 3>(triangle[1], ab_index, bc_index));
				triangl_list.push(Vec<int, 3>(triangle[2], bc_index, ca_index));
				triangl_list.push(Vec<int, 3>(ab_index, bc_index, ca_index));
				triangl_list.pop();
			}
			while (!triangl_list.empty()) {
				triangles.push_back(triangl_list.front());
				triangl_list.pop();
			}
			return real_dx;
		}
		real findMinCoord(Array<VecD> ans, int idx = 0) {
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
		VecD getMaxPosition(const Array<VecD>& points) {
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
	};


	/*
	* 		void augmentH(real aug, int side = 0) { // 0:low, 1:high
			for (int i = 0; i < points.size(); i++) {
				if (side == 0 && points[i][2] < 0)
					h[i] += aug;// *abs(points[i][2]);
				if (side == 1 && points[i][2] > 0)
					h[i] += aug;// *abs(points[i][2]);
			}
		}
		void augmentH2(real aug, int side = 0) { // 0:low, 1:high
			for (int i = 0; i < points.size(); i++) {
				if (side == 0 && points[i][2] < 0)
					h[i] = aug;// *abs(points[i][2]);
				if (side == 1 && points[i][2] > 0)
					h[i] = aug;// *abs(points[i][2]);
			}
		}
		void sinH(real aug) { // 0:low, 1:high
			for (int i = 0; i < points.size(); i++) {
				h[i] += aug * sin(points[i][0] + points[i][1] + points[i][2]);
			}
		}
		void sinH2(real aug, real range) { // 0:low, 1:high
			for (int i = 0; i < points.size(); i++) {
				if (points[i][2] > -range && points[i][2] < range)
					h[i] += aug * sin(10 * (points[i][0] + points[i][1]));
			}
		}
		void augmentGM(real aug, int side = 0) { // 0:low, 1:high
			for (int i = 0; i < points.size(); i++) {
				if (side == 0 && points[i][2] < 0)
					GM[i] += aug;
				if (side == 1 && points[i][2] > 0)
					GM[i] += aug;
			}
		}
		void augmentVorticity(real aug, real range) { // 0:low, 1:high
			for (int i = 0; i < points.size(); i++) {
				if (points[i][2]<range && points[i][2] > -range) {
					vo[i] += aug;
				}
			}
		}
		void augmentVorticity(real aug, real range1, real range2) { // 0:low, 1:high
			for (int i = 0; i < points.size(); i++) {
				if (points[i][2]<range1 && points[i][2] > range2) {
					vo[i] += aug * abs(points[i][2]);
				}
			}
		}
		void rotateV(real aug) {
			for (int i = 0; i < points.size(); i++) {
				v[i] = VectorD(-points[i][1], points[i][0], 0).normalized() * (points[i][1] * points[i][1] + points[i][0] * points[i][0]) * aug;
			}
		}
		void augmentV(real aug, int side = 0) { // 0:low, 1:high
			for (int i = 0; i < points.size(); i++) {
				if (side == 0 && points[i][2] < 0)
					v[i] += VectorD(0, 0, aug * abs(points[i][2]) * 1e-3);
				if (side == 1 && points[i][2] > 0)
					v[i] += VectorD(0, 0, aug * abs(points[i][2]) * 1e-3);
				if (side == 2 && points[i][1] < 0)
					v[i] += VectorD(0, aug * abs(points[i][2]) * 1e-3, 0);
				if (side == 3 && points[i][1] > 0)
					v[i] += VectorD(0, aug * abs(points[i][2]) * 1e-3, 0);
			}
		}
		void draw_points() {
			for (int i = 0; i < points.size(); i++) {
				glColor3f(1, 0, 0);
				glBegin(GL_POINTS);
				glVertex3f(points[i][0], points[i][1], points[i][2]);
				glEnd();
			}
		}
	* 		void draw_lines() {
			for (int i = 0; i < triangles.size(); i++) {
				glColor3f(1, 0, 0);
				glBegin(GL_LINES);
				glVertex3f(points[triangles[i][0]][0], points[triangles[i][0]][1], points[triangles[i][0]][2]);
				glVertex3f(points[triangles[i][1]][0], points[triangles[i][1]][1], points[triangles[i][1]][2]);
				glEnd();
				glBegin(GL_LINES);
				glVertex3f(points[triangles[i][2]][0], points[triangles[i][2]][1], points[triangles[i][2]][2]);
				glVertex3f(points[triangles[i][1]][0], points[triangles[i][1]][1], points[triangles[i][1]][2]);
				glEnd();
				glBegin(GL_LINES);
				glVertex3f(points[triangles[i][0]][0], points[triangles[i][0]][1], points[triangles[i][0]][2]);
				glVertex3f(points[triangles[i][2]][0], points[triangles[i][2]][1], points[triangles[i][2]][2]);
				glEnd();
			}
		}
	* 		Array<VecD> getSphereColor(Camera& camera) {
			Array<VecD> colors;
			VectorD position = VectorD(camera.camera_x, camera.camera_y, camera.camera_z);
			real gamma = 1;
			for (int i = 0; i < triangles.size(); i++) {
				VectorD normal;
				real color[3];
				getArtifactColor(position, points[triangles[i][0]], h[triangles[i][0]], gamma, color);
				colors.push_back(VectorD(color[0], color[1], color[2]));
				getArtifactColor(position, points[triangles[i][1]], h[triangles[i][1]], gamma, color);
				colors.push_back(VectorD(color[0], color[1], color[2]));
				getArtifactColor(position, points[triangles[i][2]], h[triangles[i][2]], gamma, color);
				colors.push_back(VectorD(color[0], color[1], color[2]));
			}
			return colors;
			//printf("%d %d %d\n", color[0], color[1], color[2]);

		}
void draw_triangles(Array<VectorD> color, Array<VectorD> initPoints) {
	//glColor3f(1.0, 0, 0);
	for (int i = 0; i < triangles.size(); i++) {
		if (validTriangle(points[triangles[i][0]], points[triangles[i][1]], points[triangles[i][2]], initPoints[triangles[i][0]], initPoints[triangles[i][1]], initPoints[triangles[i][2]])) {
			glBegin(GL_TRIANGLES);
			glColor3f(color[i * 3][0], color[i * 3][1], color[i * 3][2]);
			glVertex3f(points[triangles[i][0]][0], points[triangles[i][0]][1], points[triangles[i][0]][2]);
			glColor3f(color[i * 3 + 1][0], color[i * 3 + 1][1], color[i * 3 + 1][2]);
			glVertex3f(points[triangles[i][1]][0], points[triangles[i][1]][1], points[triangles[i][1]][2]);
			glColor3f(color[i * 3 + 2][0], color[i * 3 + 2][1], color[i * 3 + 2][2]);
			glVertex3f(points[triangles[i][2]][0], points[triangles[i][2]][1], points[triangles[i][2]][2]);
			glEnd();
		}
		else {
			glBegin(GL_POINTS);
			glColor3f(color[i * 3][0], color[i * 3][1], color[i * 3][2]);
			glVertex3f(points[triangles[i][0]][0], points[triangles[i][0]][1], points[triangles[i][0]][2]);
			glEnd();
			glBegin(GL_POINTS);
			glColor3f(color[i * 3 + 1][0], color[i * 3 + 1][1], color[i * 3 + 1][2]);
			glVertex3f(points[triangles[i][1]][0], points[triangles[i][1]][1], points[triangles[i][1]][2]);
			glEnd();
			glBegin(GL_POINTS);
			glColor3f(color[i * 3 + 2][0], color[i * 3 + 2][1], color[i * 3 + 2][2]);
			glVertex3f(points[triangles[i][2]][0], points[triangles[i][2]][1], points[triangles[i][2]][2]);
			glEnd();
		}
	}
}
static void getArtifactColor(VectorD camera_position, VectorD point, real h, real gamma, real color[3]) {
	VectorD normal = point; normal = normal.normalized();
	real cos_theta = (point - camera_position).normalized().dot(normal);
	//printf("-------%0.9lf %f\n", h,cos_theta);
	real thick = artifact_thick(h * 2, gamma, cos_theta);
	//printf("-------%0.9lf\n", thick);
	getColor(thick * 1e9, color);

}
*/
}}	
#endif	