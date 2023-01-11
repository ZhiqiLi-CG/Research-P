/***********************************************************/
/**	\file
	\brief		Example of persistent homology
	\details	Example of persistent homology
	\author		Zhiqi Li
	\date		12/24/2022
*/
/***********************************************************/

//	This example shows how to use GLUTWindow and FBO

#include <iostream>
#include <GL/glew.h>
#include <GL/glut.h>
#include<zqVisualization/visualization_opengl.h>
#include<zqVisualization/visualization_persistent_homology.h>
#include<zqBasicUtils/utils_timer.h>
#include<zqBasicMath/math_numerical_utils.h>
#include<zqBasicMath/math_io.h>
#include<zqBasicMath/math_homology.h>
#include<math_data_path.h>

using namespace std;

//	Two windows and a manager, ID of each window must be unique
zq::opengl::DemoWindowManager		manager;
zq::opengl::GLUTWindow3D<0>			win0(0, 0, 800, 600);
zq::opengl::GLUTWindow2D<1>			win1(820, 0, 400, 400);
zq::opengl::GLUTWindow2D<2>			win2(1240, 0, 400, 400);
// vector
std::vector<std::string> file_list{
	std::string(data_path) + std::string("/persistent_homology_one_cycle.dat"),
	std::string(data_path) + std::string("/persistent_homology_one_cycle_random.dat"),
	std::string(data_path) + std::string("/persistent_homology_two_cycle_random.dat"),
	std::string(data_path) + std::string("/persistent_homology_one_sphere_random.dat"),
};
std::vector<int> simplex_dimension{
	3,3,3,4
};

std::vector<int> witness_size{
	20,20,20,20
};

std::vector<bool> use_witness{
	false,false,false,true
};

zq::DenseVector<zq::Vec3f> points;

// homology
int example = 0;
float epsilon;
float max_epsilon = 0.6;
int reso_epsilon = 100;
std::vector<int> witness_index;

std::vector<std::pair<float, float>> epsilon_interval;
std::vector<int> feture_type;
std::vector<float> epsilon_list;
std::vector<zq::Simplical_Complex<zq::Vec3f>> complex_list;

// hci panal
float cur_epsilon = 1.0;

//	other variables
float				fps;
float* tex_ptr;

//	over-ridden functions
void calculateFps() {
	static zq::utils::FPSCalculator fps_calc;
	fps = fps_calc.GetFPS();
}

void win0_showFps() {
	glColor3f(1, 0, 1);
	zq::opengl::printInfo(0, 0, "fps: %f\nuse key a and d to adjust the epsilon, cur epsilon:%f/max epsilon:%f", fps, cur_epsilon, max_epsilon);
}

void win0_draw() {
	int complex_index = cur_epsilon / (max_epsilon / reso_epsilon);
	if (complex_index < 0) complex_index = 0;
	else if (complex_index >= complex_list.size()) complex_index = complex_list.size() - 1;
	//printf("complex_index:%d %d\n", complex_index, complex_list[complex_index].SimplexNumber());
	zq::homology::DrawSimplicalComplex(complex_list[complex_index]);
	if (use_witness[example]) {
		glColor4f(0.5, 0, 0, 1);
		for (int i = 0; i < points.Dim(); i++)
			zq::opengl::drawPointAsSphere(points.value[i], 0.01);
	}
}




void win3d_keyboard(unsigned char key, int x, int y) {
	switch (key) {
	case 27:
		exit(0);
	case 'a':
		cur_epsilon -= max_epsilon / reso_epsilon / 2;
		break;
	case 'd':
		cur_epsilon += max_epsilon / reso_epsilon / 2;
		break;
	}
}
void win1_showInfo() {
	glColor3f(0, 0, 1);
	zq::opengl::printInfo(0, 30, "persistent barcode");
}

void win2_showInfo() {
	glColor3f(0, 0, 1);
	zq::opengl::printInfo(0, 30, "persistent diagram");
}

void win1_draw() {
	//glEnable(GL_TEXTURE_2D);
	float display_max_epsilon = max_epsilon * 1.1;
	float win_x = 2;
	float win_y = 2;
	float win_ox = -1;
	float win_oy = -1;
	zq::homology::DrawPersistentBarCode(
		epsilon_interval,
		feture_type,
		max_epsilon,
		display_max_epsilon,
		win_x,
		win_y,
		win_ox,
		win_oy
	);
	// Then draw the line
	glColor3f(1, 0, 0);
	zq::opengl::drawLineSegment(
		zq::Vec2f(cur_epsilon / display_max_epsilon * win_x + win_ox, win_oy),
		zq::Vec2f(cur_epsilon / display_max_epsilon * win_x + win_ox, win_oy + win_y)
	);
}
void win2_draw() {
	//glEnable(GL_TEXTURE_2D);
	float display_max_epsilon = max_epsilon * 1.1;
	float win_x = 2;
	float win_y = 2;
	float win_ox = -1;
	float win_oy = -1;
	zq::homology::DrawPersistentDiagram(
		epsilon_interval,
		feture_type,
		max_epsilon,
		display_max_epsilon,
		max_epsilon * 0.1,
		win_x,
		win_y,
		win_ox,
		win_oy
	);
}
int main() {
	// 1. read points
	cur_epsilon = 0;

	printf("the example:");
	scanf("%d", &example);
	zq::readDenseVectorPointsFromFile(file_list[example].c_str(), points);
	// 2. set the epsilon list
	for (int i = 0; i < reso_epsilon; i++) {
		epsilon_list.push_back(max_epsilon / reso_epsilon * i);
	}
	if (use_witness[example]) {
		// 3. set the complex list
		// 3.1 find witness
		zq::FindWitness(
			points.Dim(),
			witness_size[example],
			witness_index
		);
		// 3.2 construct
		for (int i = 0; i < epsilon_list.size(); i++) {
			std::vector<std::vector<int>> results;
			zq::VRWitnessComplexConstruct(
				epsilon_list[i],
				points.value,
				witness_index,
				simplex_dimension[example],
				results
			);
			complex_list.push_back(zq::Simplical_Complex<zq::Vec3f>(&(points.value[0]), points.Dim(), results));
		}
	}
	else {
		// 3. set the complex list
		for (int i = 0; i < epsilon_list.size(); i++) {
			std::vector<std::vector<int>> results;
			zq::VRComplexConstruct(
				epsilon_list[i],
				points.value,
				simplex_dimension[example],
				results
			);
			complex_list.push_back(zq::Simplical_Complex<zq::Vec3f>(&(points.value[0]), points.Dim(), results));
		}

	}
	// 4. calculate the persistent diagram

	zq::CalculatePersistentDataSparse(
		epsilon_list,
		max_epsilon,
		complex_list,
		epsilon_interval,
		feture_type
	);
	auto tem_epsilon_interval = epsilon_interval;
	auto tem_feture_type = feture_type;
	feture_type.clear();
	epsilon_interval.clear();
	int maxv = -1;
	for (int i = 0; i < tem_feture_type.size(); i++) {
		maxv = maxv > tem_feture_type[i] ? maxv : tem_feture_type[i];
		if (tem_feture_type[i] != simplex_dimension[example]) {
			epsilon_interval.push_back(tem_epsilon_interval[i]);
			feture_type.push_back(tem_feture_type[i]);
		}
	}
	std::cout << maxv << std::endl;
	//	create window 0
	win0.SetDrawAppend(win0_showFps);				//	set draw append callback
	win0.SetDraw(win0_draw);						//	set draw callback
	win0.SetBackgroundColor(0, 0.4, 0.5, 0);		//	set background color
	win0.keyboardFunc = win3d_keyboard;
	win0.CreateGLUTWindow();						//	create window 0

	//	create window 1, then setup barcode
	win1.SetDraw(win1_draw);							//	set draw callback
	win1.SetDrawAppend(win1_showInfo);					//	set draw append callback
	win1.SetBackgroundColor(0, 1, 1, 1);				//	set background color
	win1.CreateGLUTWindow();							//	create window 1

	//	create window 1, then setup persistent diagram
	win2.SetDraw(win2_draw);							//	set draw callback
	win2.SetDrawAppend(win2_showInfo);					//	set draw append callback
	win2.SetBackgroundColor(0, 1, 1, 1);				//	set background color
	win2.CreateGLUTWindow();							//	create window 1

	//!	[Window Manager Setup]
	manager.AddIdleFunc(calculateFps);	//	add each idle function in sequence
	manager.AddIdleFunc(win0.idleFunc);
	manager.AddIdleFunc(win1.idleFunc);
	manager.AddIdleFunc(win2.idleFunc);

	manager.EnterMainLoop();			//	then enter main loop
	//!	[Window Manager Setup]
}