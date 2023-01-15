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
	zq::real bubble_radius, stop_dx;
	virtual void settingExample() {
		zq::json j;
		j["test"] = 0;
		j["output_dir"] = "output";
		j["last_frame"] = 10000;
		j["frame_rate"] = 200;
		j["snapshot_stride"] = 1;
		j["const_value"] = bubble.settingExample();
		j["stop_dx"] = 0.03;
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
			zq::Info("enter test case {}", test);
			switch (test) {
				case 0: Case_0(); break; /// Bubble
			}
		}
	}
	virtual void writeOutputFiles(const int frame)
	{
		zq::Info("frame {}", frame);
		Base::writeOutputFiles(frame);
		if (frame == 0) {
			zq::json j;
			j["object"] = {
				{"type","triangle_cloud"},
				{"points","points.bin"},
				{"triangle","triangle.bin"},
				{"h","h.bin"},
				{"n","n.bin"},
				{"v","v.bin"},
				{"f","f.bin"},
			};
			std::ofstream fout(frame_dir + "/format.json");
			fout << j;
		}
		zq::Array<zq::Vec<zq::real,d>, HOST> host_x = bubble.x;
		zq::Array<zq::Vec<zq::real, d>, HOST> host_v = bubble.v;
		zq::Array<zq::Vec<zq::real, d>, HOST> host_f = bubble.f;
		zq::Array<zq::Mat<zq::real, d>, HOST> host_e = bubble.e;
		zq::Array<zq::real, HOST> host_h = bubble.h;
		zq::Array<zq::Vec<zq::real, d>, HOST> host_n(host_e.size());
		for (int i = 0; i < host_e.size(); i++) {
			host_n[i] = host_e[i].col(0);
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
			vo[i] = 0;

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
			});
		bubble.init(sphere.maxPosition, sphere.minPosition);
	}
	
	
};

#endif