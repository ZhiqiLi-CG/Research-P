/***********************************************************/
/**	\file
	\brief		Array Check Functions
	\author		Zhiqi Li, based on yzLib of Dr.Yizhong Zhang
	\date		9/28/2012
*/
/***********************************************************/
#ifndef __UTILS_ARRAY_CHECK_H__
#define __UTILS_ARRAY_CHECK_H__

#include <iostream>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <cstdlib>
#include <string.h>
#include<zqBasicUtils/utils_hash.h>
#include<zqBasicUtils/utils_timer.h>
#include<zqBasicUtils/utils_io.h>
#include<zqBasicUtils/utils_json.h>
#include<omp.h>

namespace zq{ namespace physics{

	class PhysicsSimulation {
	public:
		int test = 1;
		std::string output_dir = "output";
		std::string frame_dir;
		int first_frame = 0, last_frame = 250, current_frame = 0;
		int current_step = 0;
		real frame_rate = 25;
		real time = (real)0, current_time = (real)0;
		real cfl = (real)1;
		int max_iter_per_frame = 2000;//set to -1 to disable max iteration limit
		bool verbose = false;
		int snapshot_stride = 1;
		zq::json const_value;
		std::string setting_file;

		public:
		void setThreads(int number_threads) {
			omp_set_num_threads(number_threads);
			int max_threads = omp_get_max_threads();
			Info("Set {} threads, run with {} cores", number_threads, max_threads);
		}
		real timeAtFrame(const int frame) { 
			return (real)frame / frame_rate; 
		}
		int frameAtTime(const real time) { 
			return (int)((real)time * frame_rate); 
		}
		virtual real CFL() const { return cfl; }

		virtual void initialize() {}

		virtual void run(){
			zq::utils::Timer timer;
			timer.status = 0x01;// ticking
			if (current_frame == 0) writeOutputFiles(current_frame);
			while (current_frame < last_frame) {
				printf("%d %d begin\n", current_frame, last_frame);
				timer.Restart();
				current_frame++;
				advanceToTargetTime(timeAtFrame(current_frame));
				double frame_time = timer.Elapsed();
				Info("Frame time: {:.2f}s", frame_time);
				if(current_frame% snapshot_stride ==0)
					writeOutputFiles(current_frame);
				printf("%d %d end\n", current_frame, last_frame);
			}
		}

		virtual void load(const int frame)
		{
			current_frame = frame;
			std::cout << "Load Output file from frame: " << frame << std::endl;
			loadOutputFiles(current_frame);
		}

		virtual void advanceToTargetTime(const real target_time)
		{
			bool done = false;
			for (int substep = 1; !done; substep++) {
				real dt = this->CFL();
				if (time + dt >= target_time) { dt = target_time - time; done = true; }
				else if (time + 2 * dt >= target_time) { dt = (real).5 * (target_time - time); }
				advanceOneTimeStep(dt, time);
				time += dt;
			}
		}

		virtual void advanceOneTimeStep(const real dt, const real time) {}

		virtual void writeOutputFiles(const int frame)
		{
			if (frame == 0) {
				if (!zq::utils::Directory_Exists(output_dir.c_str()))
					zq::utils::Create_Directory(output_dir);
			}

			frame_dir = output_dir + "/" + std::to_string(frame);
			if (!zq::utils::Directory_Exists(frame_dir.c_str())) zq::utils::Create_Directory(frame_dir);
			{
				std::string file_name = output_dir + "/0/last_frame.txt";
				zq::utils::Write_Text_To_File(file_name, std::to_string(frame));
			}

			{
				std::string file_name = frame_dir + "/time";
				zq::utils::Write_Binary_To_File<real>(file_name, time); 
			}

			std::cout << "#     Write Frame " << frame << " to: " << frame_dir << std::endl;
		}

		virtual void loadOutputFiles(const int frame) {
			frame_dir = output_dir + "/" + std::to_string(frame);
			
			{
				std::string file_name = frame_dir + "/time";
				zq::utils::Read_Binary_From_File<real>(file_name, time); 
			}
			
			if (verbose)std::cout << "Load output files for frame " << frame << std::endl;
		}

		//////////////////////////////////////////////////////////////////////////
		////streaming the file IO
		virtual void runStream()
		{
			while (current_frame < last_frame) {
#pragma omp parallel sections
				{
#pragma omp section
					{writeOutputFiles(current_frame); }

#pragma omp section
					{advanceToTargetTime(timeAtFrame(current_frame + 1)); }
				}
				current_frame++;
			}
			writeOutputFiles(current_frame);
		}

		virtual void updateIOBufferStream(const int frame) {}

		virtual void advanceOneTimeStepStream(const real dt, const real time) {}

		virtual void writeOutputFilesStream(const int frame) { writeOutputFiles(frame); }

		virtual void settingExample() {
			json j;
			j["test"] = 0;
			j["output_dir"] = "output";
			j["last_frame"] = 10000;
			j["frame_rate"] = 200;
			j["snapshot_stride"] = 1;
			j["const_value"] = json();
			std::ofstream jfile(setting_file.c_str());
			jfile << j.dump(2);
		}

		virtual void setSettingFile(std::string setting_file) {
			this->setting_file = setting_file;
		}
	};

}}	
#endif	