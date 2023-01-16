#ifndef __SIMPLE_BUBBLE_SIMULATION_H__
#define __SIMPLE_BUBBLE_SIMULATION_H__
#include<zqPhysics/physics_simulation.h>
#include<zqPhysics/geometry_init.h>
#include<zqPhysics/codim1_sph_object.h>
#include<zqBasicMath/math_io.h>
#include<zqBasicUtils/utils_io.h>
#include<physics_data_path.h>
template<int d, int side>
class SimpleBubbleSimulation :public zq::physics::PhysicsSimulation {
	using Base = zq::physics::PhysicsSimulation;
public:
	int mode = 0;
	zq::physics::SimpleFilm<d, side> bubble;
	int N;
	zq::Array<zq::Vec<int, 3>> triangle;
	zq::real bubble_radius, stop_dx, transparent;
	virtual void settingExample() {
		zq::json j;
		j["test"] = 0;
		j["output_dir"] = "output";
		j["last_frame"] = 10000;
		j["frame_rate"] = 200;
		j["snapshot_stride"] = 1;
		j["const_value"] = bubble.settingExample();
		j["stop_dx"] = 0.03;
		j["transparent"] = 0.5;
		j["bubble_radius"] = 1;
		std::ofstream jfile(setting_file.c_str());
		jfile << j.dump(2);
	}
	virtual void advanceOneTimeStep(const zq::real dt, const zq::real time) {
		zq::Info("[Frame] dt, time:{} {}", dt, time);
		bubble.update(dt);
	}
	virtual void initialize()
	{
		if (mode == -1) {
			settingExample();
			return;
		}
		else {
			zq::json j;
			std::ifstream jfile(setting_file.c_str());
			jfile >> j;
			int test = j["test"];
			output_dir = std::string(physics_data_path) + std::string("/")+std::string(j["output_dir"]);
			last_frame = j["last_frame"];
			frame_rate = j["frame_rate"];// 200;
			verbose = true;
			cfl = 0.1;
			snapshot_stride = j["snapshot_stride"];
			const_value = j["const_value"];
			stop_dx = j["stop_dx"];
			bubble_radius = j["bubble_radius"];
			transparent = j["transparent"];
			zq::Info("enter test case {}", test);
			switch (test) {
				case 0: Case_0(); break; /// Bubble cap
				case 1: Case_1(); break; /// Bubble mag
				case 2: Case_2(); break; /// Bubble vor
				case 3: Case_3(); break; /// Bubble pressure
			}
		}
	}
	virtual void writeOutputFiles(const int frame)
	{
		zq::Info("frame {}", frame);
		Base::writeOutputFiles(frame);
		if (frame == 0) {
			zq::json j;
			j ={
					{"type","point_cloud"},
					{"dim",3},
					{"points",{
						{"name","points.bin"},
						{"type","base"},
						{"dim",d},
						{"data","real"}
					}},
					{"triangle",{
						{"name","triangle.bin"},
						{"type","triangle"},
						{"dim",3},
						{"data","int"}
					}},
					{"h",{
						{"name","h.bin"},
						{"type","color"},
						{"dim",1},
						{"data","real"},
						{"color_type","film"},
						{"transparent",transparent}
					}},
					{"n",{
						{"name","n.bin"},
						{"type","vector"},
						{"dim",d},
						{"data","real"}
					}},
					{"v",{
						{"name","v.bin"},
						{"type","vector"},
						{"dim",d},
						{"data","real"}
					}},
					{"f",{
						{"name","f.bin"},
						{"type","vector"},
						{"dim",d},
						{"data","real"}
					}}
			};
			std::ofstream fout(frame_dir + "/format.json");
			fout << j.dump(2);
		}
		zq::Array<zq::Vec<zq::real,d>, HOST> host_x = bubble.x;
		zq::Array<zq::Vec<zq::real, d>, HOST> host_v = bubble.v;
		zq::Array<zq::Vec<zq::real, d>, HOST> host_f = bubble.f;
		zq::Array<zq::Mat<zq::real, d>, HOST> host_e = bubble.e;
		zq::Array<zq::real, HOST> host_h = bubble.ah;
		zq::Array<zq::Vec<zq::real, d>, HOST> host_n(host_e.size());
		for (int i = 0; i < host_e.size(); i++) {
			host_n[i] = host_e[i].col(d-1);
		}
		zq::utils::Write_Vector_Array_3D< zq::real,d > (frame_dir + "/points.bin", host_x);
		zq::utils::Write_Vector_Array_3D< int, d >(frame_dir + "/triangle.bin", triangle);
		zq::utils::Write_Vector_Array_3D< zq::real, d >(frame_dir + "/v.bin", host_v);
		zq::utils::Write_Vector_Array_3D< zq::real, d >(frame_dir + "/f.bin", host_f);
		zq::utils::Write_Vector_Array_3D< zq::real, d >(frame_dir + "/n.bin", host_n);
		zq::utils::Write_Array < zq::real > (frame_dir + "/h.bin", host_h);
	}
	void Case_0() {
		zq::physics::InitSphere<d> sphere(bubble_radius, stop_dx);
		triangle = sphere.triangles;
		zq::Array<zq::Vec<zq::real,d>> v;
		zq::Array<zq::real> h;
		zq::Array<zq::real> GM;
		zq::Array<zq::real> vo;
		bubble.template  Init_Array<zq::Vec<zq::real,d>,HOST>("x", sphere.points);
		bubble.template Init_Array<zq::Vec<zq::real, d>, HOST>("init_n", sphere.normals);
		bubble.template Init_Array<zq::Mat<zq::real, d-1>, HOST>("g", sphere.g);
		v.resize(sphere.points.size());
		h.resize(sphere.points.size());
		GM.resize(sphere.points.size());
		vo.resize(sphere.points.size());
		for (int i = 0; i < sphere.points.size(); i++) {
			v[i] = zq::Vec<zq::real,d>();
			h[i] = const_value["h_0"];
			GM[i] = 1e-7;
			//if (sphere.points[i].z < bubble_radius * 0.9);
			//	vo[i] = 0.0001;

		}
		/*
		zq::real h_aug= 5e-8,range=bubble_radius*0.1;
		for (int i = 0; i < h.size(); i++) {
			zq::real theta = atan(sphere.points[i][1] / sphere.points[i][0]);
			zq::real this_range = range * abs(cos(theta*4));
			if (abs(sphere.points[i][2]) < this_range)
				h[i] = const_value["h_0"] + h_aug;// (1 - abs(sphere.points[i][2]) / this_range)* h_aug;
			if (sphere.points[i][2] > 0)
				h[i] = const_value["h_0"]+h_aug;// *sin(10 * (sphere.points[i][0] + sphere.points[i][1]));
//			h[i] += h_aug * sin(sphere.points[i][0] + sphere.points[i][1] + sphere.points[i][2]);
		}
		*/
		zq::real gm_aug = 1e-7, range = bubble_radius * 0.3;
		for (int i = 0; i < GM.size(); i++) {
			zq::real theta = atan(sphere.points[i][1] / sphere.points[i][0]);
			zq::real this_range = range * abs(cos(theta * 4))*zq::randNumber(0.7,1.2);
			if (abs(sphere.points[i][2]) < this_range)
				GM[i] = 1e-7 + gm_aug * zq::randNumber(0.7, 1.2);// (1 - abs(sphere.points[i][2]) / this_range)* h_aug;
			if (sphere.points[i][2] > 0)
				GM[i] = 1e-7 + gm_aug;// *sin(10 * (sphere.points[i][0] + sphere.points[i][1]));
//			h[i] += h_aug * sin(sphere.points[i][0] + sphere.points[i][1] + sphere.points[i][2]);
		}
		/*
		zq::real gm_aug = 5e-8;
		int gm_side = 2;
		for (int i = 0; i < GM.size(); i++) {
			if (gm_side == 0 && sphere.points[i][2] < 0)
				GM[i] += gm_aug;
			if (gm_side == 1 && sphere.points[i][2] > 0)
				GM[i] += gm_aug;
		}*/
		bubble.template Init_Array<zq::Vec<zq::real, d>, HOST>("v", v);
		bubble.template Init_Array<zq::real, HOST>("h", h);
		bubble.template Init_Array<zq::real, HOST>("GM", GM);
		bubble.template Init_Array<zq::real, HOST>("vo", vo);

		

		bubble.Init_Constant(
			std::unordered_map<std::string, zq::real>{
				{"t_dot", const_value["t_dot"]},
				{ "t_r",const_value["t_r"] },
				{ "v_r",const_value["v_r"] },
				{ "gravity",const_value["gravity"] },
				{ "gamma_0",const_value["gamma_0"] },
				{ "gamma_a",const_value["gamma_a"] },
				{ "alpha_h",const_value["alpha_h"] },
				{ "alpha_k", const_value["alpha_k"] },
				{ "alpha_d", const_value["alpha_d"] },
				{ "alpha_c", const_value["alpha_c"] },
				{ "h_0", const_value["h_0"] },
				{ "p_0", const_value["p_0"] },
				{ "mu", const_value["mu"] },
				{ "rho", const_value["rho"] },
				{"cap_coef",const_value["cap_coef"]}
			});
		bubble.init(sphere.maxPosition, sphere.minPosition);
	}
	void Case_1() {
		zq::physics::InitSphere<d> sphere(bubble_radius, stop_dx);
		triangle = sphere.triangles;
		zq::Array<zq::Vec<zq::real, d>> v;
		zq::Array<zq::real> h;
		zq::Array<zq::real> GM;
		zq::Array<zq::real> vo;
		bubble.template  Init_Array<zq::Vec<zq::real, d>, HOST>("x", sphere.points);
		bubble.template Init_Array<zq::Vec<zq::real, d>, HOST>("init_n", sphere.normals);
		bubble.template Init_Array<zq::Mat<zq::real, d - 1>, HOST>("g", sphere.g);
		v.resize(sphere.points.size());
		h.resize(sphere.points.size());
		GM.resize(sphere.points.size());
		vo.resize(sphere.points.size());
		for (int i = 0; i < sphere.points.size(); i++) {
			v[i] = zq::Vec<zq::real, d>();
			h[i] = const_value["h_0"];
			GM[i] = 1e-7;
		}
		zq::real gm_aug = 1e-7, range = bubble_radius * 0.3;
		for (int i = 0; i < GM.size(); i++) {
			zq::real theta = atan(sphere.points[i][1] / sphere.points[i][0]);
			zq::real this_range = range * abs(cos(theta * 4)) * zq::randNumber(0.7, 1.2);
			if (abs(sphere.points[i][2]) < this_range)
				GM[i] = 1e-7 + gm_aug * zq::randNumber(0.7, 1.2);// (1 - abs(sphere.points[i][2]) / this_range)* h_aug;
			if (sphere.points[i][2] > 0)
				GM[i] = 1e-7 + gm_aug;// *sin(10 * (sphere.points[i][0] + sphere.points[i][1]));
		}
		bubble.template Init_Array<zq::Vec<zq::real, d>, HOST>("v", v);
		bubble.template Init_Array<zq::real, HOST>("h", h);
		bubble.template Init_Array<zq::real, HOST>("GM", GM);
		bubble.template Init_Array<zq::real, HOST>("vo", vo);



		bubble.Init_Constant(
			std::unordered_map<std::string, zq::real>{
				{"t_dot", const_value["t_dot"]},
				{ "t_r",const_value["t_r"] },
				{ "v_r",const_value["v_r"] },
				{ "gravity",const_value["gravity"] },
				{ "gamma_0",const_value["gamma_0"] },
				{ "gamma_a",const_value["gamma_a"] },
				{ "alpha_h",const_value["alpha_h"] },
				{ "alpha_k", const_value["alpha_k"] },
				{ "alpha_d", const_value["alpha_d"] },
				{ "alpha_c", const_value["alpha_c"] },
				{ "h_0", const_value["h_0"] },
				{ "p_0", const_value["p_0"] },
				{ "mu", const_value["mu"] },
				{ "rho", const_value["rho"] },
				{ "cap_coef",const_value["cap_coef"] }
		});
		bubble.init(sphere.maxPosition, sphere.minPosition);
	}
	void Case_2() {
		zq::physics::InitSphere<d> sphere(bubble_radius, stop_dx);
		triangle = sphere.triangles;
		zq::Array<zq::Vec<zq::real, d>> v;
		zq::Array<zq::real> h;
		zq::Array<zq::real> GM;
		zq::Array<zq::real> vo;
		bubble.template  Init_Array<zq::Vec<zq::real, d>, HOST>("x", sphere.points);
		bubble.template Init_Array<zq::Vec<zq::real, d>, HOST>("init_n", sphere.normals);
		bubble.template Init_Array<zq::Mat<zq::real, d - 1>, HOST>("g", sphere.g);
		v.resize(sphere.points.size());
		h.resize(sphere.points.size());
		GM.resize(sphere.points.size());
		vo.resize(sphere.points.size());
		for (int i = 0; i < sphere.points.size(); i++) {
			v[i] = zq::Vec<zq::real, d>();
			h[i] = const_value["h_0"];
			GM[i] = 1e-7;
			if(sphere.points[i][2]<-bubble_radius*0.8)
				vo[i] = 3e-5;
		}
		bubble.template Init_Array<zq::Vec<zq::real, d>, HOST>("v", v);
		bubble.template Init_Array<zq::real, HOST>("h", h);
		bubble.template Init_Array<zq::real, HOST>("GM", GM);
		bubble.template Init_Array<zq::real, HOST>("vo", vo);

		bubble.Init_Constant(
			std::unordered_map<std::string, zq::real>{
				{"t_dot", const_value["t_dot"]},
				{ "t_r",const_value["t_r"] },
				{ "v_r",const_value["v_r"] },
				{ "gravity",const_value["gravity"] },
				{ "gamma_0",const_value["gamma_0"] },
				{ "gamma_a",const_value["gamma_a"] },
				{ "alpha_h",const_value["alpha_h"] },
				{ "alpha_k", const_value["alpha_k"] },
				{ "alpha_d", const_value["alpha_d"] },
				{ "alpha_c", const_value["alpha_c"] },
				{ "h_0", const_value["h_0"] },
				{ "p_0", const_value["p_0"] },
				{ "mu", const_value["mu"] },
				{ "rho", const_value["rho"] },
				{ "cap_coef",const_value["cap_coef"] }
		});
		bubble.init(sphere.maxPosition, sphere.minPosition);
	}
	void Case_3() {
		zq::physics::InitSphere<d> sphere(bubble_radius, stop_dx);
		triangle = sphere.triangles;
		zq::Array<zq::Vec<zq::real, d>> v;
		zq::Array<zq::real> h;
		zq::Array<zq::real> GM;
		zq::Array<zq::real> vo;
		bubble.template  Init_Array<zq::Vec<zq::real, d>, HOST>("x", sphere.points);
		bubble.template Init_Array<zq::Vec<zq::real, d>, HOST>("init_n", sphere.normals);
		bubble.template Init_Array<zq::Mat<zq::real, d - 1>, HOST>("g", sphere.g);
		v.resize(sphere.points.size());
		h.resize(sphere.points.size());
		GM.resize(sphere.points.size());
		vo.resize(sphere.points.size());
		for (int i = 0; i < sphere.points.size(); i++) {
			v[i] = zq::Vec<zq::real, d>();
			h[i] = const_value["h_0"];
			GM[i] = 1e-7;
			if (sphere.points[i][2] < -bubble_radius * 0.5) {
				zq::real ratio = (abs(sphere.points[i][2]) - bubble_radius * 0.5) / bubble_radius * 0.5;
				if constexpr (d == 3) {
					v[i] = zq::Vec<zq::real, d>(0, 0, ratio * 3e1);
				}
				else if constexpr (d == 2) {
					v[i] = zq::Vec<zq::real, d>(0, ratio * 3e1);
				}
			}
		}
		bubble.template Init_Array<zq::Vec<zq::real, d>, HOST>("v", v);
		bubble.template Init_Array<zq::real, HOST>("h", h);
		bubble.template Init_Array<zq::real, HOST>("GM", GM);
		bubble.template Init_Array<zq::real, HOST>("vo", vo);

		bubble.Init_Constant(
			std::unordered_map<std::string, zq::real>{
				{"t_dot", const_value["t_dot"]},
				{ "t_r",const_value["t_r"] },
				{ "v_r",const_value["v_r"] },
				{ "gravity",const_value["gravity"] },
				{ "gamma_0",const_value["gamma_0"] },
				{ "gamma_a",const_value["gamma_a"] },
				{ "alpha_h",const_value["alpha_h"] },
				{ "alpha_k", const_value["alpha_k"] },
				{ "alpha_d", const_value["alpha_d"] },
				{ "alpha_c", const_value["alpha_c"] },
				{ "h_0", const_value["h_0"] },
				{ "p_0", const_value["p_0"] },
				{ "mu", const_value["mu"] },
				{ "rho", const_value["rho"] },
				{ "cap_coef",const_value["cap_coef"] }
		});
		bubble.init(sphere.maxPosition, sphere.minPosition);
	}
};

#endif