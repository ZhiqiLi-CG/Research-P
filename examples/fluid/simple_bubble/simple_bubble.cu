#include<ResearchP_config.h>
#ifdef  RESEARCHP_ENABLE_CUDA
#include <zqBasicUtils/utils_cuda.h>
#include<simple_bubble_physics.h>
#include<physics_data_path.h>
#include<string>
#include<vector>
std::vector<std::string> setting_file_list{
	std::string(physics_data_path) + "/simple_bubble_capillary.json",
	std::string(physics_data_path) + "/simple_bubble_marangoni.json",
	std::string(physics_data_path) + "/simple_bubble_vorticity.json",
	std::string(physics_data_path) + "/simple_bubble_pressure.json",
	std::string(physics_data_path) + "/simple_bubble_cuda.json"
};

int main() {
	int mode = 0;
	int test = 0;
	printf("mode(-1:generate config example,0:reand):");
	scanf("%d", &mode);
	printf("\n test:");
	scanf("%d", &test);
	if (test != -1) {
		SimpleBubbleSimulation<3, DEVICE> bubble_simulation;
		bubble_simulation.setSettingFile(setting_file_list[test]);
		bubble_simulation.mode = mode;
		bubble_simulation.initialize();
		bubble_simulation.setThreads(10);
		bubble_simulation.run();
	}
}
#endif //  RESEARCHP_ENABLE_CUDA