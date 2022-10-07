#include <cstdio>
#include <cmath>
#include <unistd.h>
#include <cstdlib>
#include <GL/glew.h>
#include <GLFW/glfw3.h>

//Measurements in kilograms and meters
//Angles in radians, rotation in radians per second
#define G_CONSTANT 6.674e-11
#define PI 3.1415
#define TICK 6000
#define TEST_MASS_1 10e11
#define TEST_MASS_2 10e11
#define TEST_MASS_3 10e9
//#define SCALE 2e12
#define RAD_ADJ PI / 180.0
#define MASS_EARTH 5.972e24
#define MASS_SUN 1.989e30
#define MASS_MOON 7.342e22
#define MASS_JUPITER 1.899e27
#define MASS_MARS 6.39e23
#define DISTANCE_EM 3.844e8
#define DISTANCE_ES 1.5e11
#define DISTANCE_JS 7.41e11
#define DISTANCE_MS 2.16e11
#define I_VEL_E 2.98e4
#define I_VEL_M 1.022e3
#define I_VEL_MARS 2.4e4
#define I_VEL_J 1.307e4
#define T_SLEEP 66000
#define SIGNIFICANT_RADIANS 0.1
#define R_EARTH 6.371e6
#define R_SUN 6.963e8
#define R_MOON 1.737e6
#define R_MARS 3.389e6
#define R_JUPITER 6.991e7
#define NUM_VERTICES_PER 18
#define WIDTH 1000
#define HEIGHT 1000
#define SCALE (DISTANCE_JS*1.05)
#define SCALE_UP 200
#define ROT_EARTH 7.272e-5 
#define ROT_SUN 2.972e-6 
#define ROT_MOON (ROT_EARTH * 3.385e-2) 
#define ROT_MARS (ROT_EARTH * 0.9715)
#define ROT_JUPITER (ROT_EARTH * 2.39)
#define MAX_ROTATION (PI * 1e16)

//Include after define constants
#include "mass_vectors.h"

const char *vertexShaderSource = "#version 330 core\n"
"layout (location = 0) in vec3 aPos;\n"
"void main()\n"
"{\n"
" gl_Position = vec4(aPos.x, aPos.y, aPos.z, 1.0);\n"
"}\0";

const char *fragmentShaderSource = "#version 330 core\n"
"out vec4 FragColor;\n"
"void main()\n"
"{\n"
"FragColor = vec4(1.0f, 0.5f, 0.2f, 1.0f);\n"
"}\0";

void framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
    // make sure the viewport matches the new window dimensions; note that width and 
    // height will be significantly larger than specified on retina displays.
    glViewport(0, 0, width, height);
}

void processInput(GLFWwindow *window)
{
    if(glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
        glfwSetWindowShouldClose(window, true);
}



/*
 * Formulas:
 * Newtonian gravity F = G (m1*m2)/(r^2)
 * Formulas for vector components
 * x = r (cos(al)*sin(be))/(sqrt(cos^2(be)*sin^2(al)*sin^2(be)))
 * y = r (cos(be)*sin(al)..etc
 * z = r (sin(al)*sin(be))/etc
 * where r is the magnitude of the vector and alpha is the angle measured from x along the z axis
 * and beta is measured from y along the z axis
 * z/x = tan alpha
 * z/y = tan beta
 * Sin theta = opp/hyp
 * Cos theta = adj/hyp
 * Tan theta = opp/adj
 * Useful math functions: sin, cos, tan, acos, asin, atan
 * fabs, abs, 
 */

double time_passed_total = 0;
/*
 * Find the vector created by one tick time of acceleration due to gravity for mo1
 * break it into components x, y, z
 * determine the sign of each component based on coordinates of mo1 and mo2
 */
struct vector_3 vector_acc_per_tick(struct massive_object mo1, struct massive_object mo2)
{
	struct vector_3 v_acc;
	v_acc.magnitude = acceleration_due_to_gravity(mo1, mo2) * TICK;
	find_angle_vector(&v_acc, mo1, mo2);
	component_3d_alt(&v_acc);
	
	if (mo1.loc_x > mo2.loc_x)
		v_acc.x_comp *= -1.0;
	if (mo1.loc_y > mo2.loc_y)
		v_acc.y_comp *= -1.0;
	if (mo1.loc_z > mo2.loc_z)
		v_acc.z_comp *= -1.0;
	//printf("Magnitude total %e\n", v_acc.magnitude);
	//printf("Vector x_comp: %f\n", v_acc.x_comp);
	//printf("Vector y_comp: %f\n", v_acc.y_comp);
	//printf("Vector z_comp: %f\n", v_acc.z_comp);

	return v_acc;
}

//Combine vectors by component values
struct vector_3 add_vectors(struct vector_3 *initial, struct vector_3 *v_acc)
{
	struct vector_3 v_final;
	v_final.x_comp = initial->x_comp + v_acc->x_comp;
	v_final.y_comp = initial->y_comp + v_acc->y_comp;
	v_final.z_comp = initial->z_comp + v_acc->z_comp;

	//printf("V_INIT: [%f, %f, %f]\n", initial->x_comp, initial->y_comp, initial->z_comp);
	//printf("V_FINAL: [%f, %f, %f]\n", v_final.x_comp, v_final.y_comp, v_final.z_comp);
	
	return v_final;
}

void update_location(struct mass_list *np)
{
	np->m_obj.loc_x += np->m_obj.vect.x_comp * TICK;
	np->m_obj.loc_y += np->m_obj.vect.y_comp * TICK;
	np->m_obj.loc_z += np->m_obj.vect.z_comp * TICK;
}


//Takes an index to decide which object will be calculating its new vector
void gravity_list_iterate(struct mass_list *root, int index)
{
	struct mass_list *np = root;
	for (int i = 0; i < index; i++){
		if (np->next != nullptr)
			np = np->next;
	}
	struct mass_list *tp = root;
	struct vector_3 v_change;
	for (int i = 0; tp !=nullptr; i++){
		if (tp != np){
			v_change = vector_acc_per_tick((np->m_obj), (tp->m_obj));
			np->m_obj.vect = add_vectors(&np->m_obj.vect, &v_change);
		}
		tp = tp->next; 
	}
	//Update rotation
	if ((np->m_obj.rotation_cur + np->m_obj.rotation * TICK) > MAX_ROTATION)
		np->m_obj.rotation_cur = np->m_obj.rotation;
	np->m_obj.rotation_cur += np->m_obj.rotation * TICK;
		
}

float *genVertices(struct massive_object mo)
{
	//24 float values 3 for each 8 vertices
	float *vertices = (float *) malloc(sizeof(float)*NUM_VERTICES_PER);
	float *vp = vertices;
	int scale_up = SCALE_UP;
	double cos_alongz = cos(mo.rotation_cur);
	double sin_alongz = sin(mo.rotation_cur);
	if (mo.radius == R_SUN)
		scale_up = 25;
	double scale_xprad = (mo.loc_x + (mo.radius * scale_up)) / SCALE;
	double scale_xmrad = (mo.loc_x - (mo.radius * scale_up)) / SCALE;
	double scale_yprad = (mo.loc_y + (mo.radius * scale_up)) / SCALE;
	double scale_ymrad = (mo.loc_y - (mo.radius * scale_up)) / SCALE;
	double scale_zprad = (mo.loc_z + (mo.radius * scale_up)) / SCALE;
	double scale_zmrad = (mo.loc_z - (mo.radius * scale_up)) / SCALE;
	double scale_x = mo.loc_x / SCALE;
	double scale_y = mo.loc_y / SCALE;
	double scale_z = mo.loc_z / SCALE;
	double local_x = scale_xprad - scale_x;
	double local_xm = scale_xmrad - scale_x;
	double local_y = scale_yprad - scale_y;
	double local_ym = scale_ymrad - scale_y;
	double local_z = 0.0;
	//least x, neutral y, neutral z
	*vertices = (float) scale_x/*mrad*/ + (local_xm * cos_alongz) - (sin_alongz * 0);
	vp++;
	*(vp++) = (float) scale_y + (0 * cos_alongz) + (sin_alongz * local_xm);
	*(vp++) = (float) scale_z; 
	//least y
	*(vp++) = (float) scale_x + (0 * cos_alongz) - (sin_alongz * local_ym);
	*(vp++) = (float) scale_y/*mrad*/ + (local_ym * cos_alongz) + (sin_alongz * 0);
	*(vp++) = (float) scale_z;
	//most y
	*(vp++) = (float) scale_x + (0 * cos_alongz) - (sin_alongz * local_y);
	*(vp++) = (float) scale_y/*prad*/ + (local_y * cos_alongz) + (sin_alongz * 0);
	*(vp++) = (float) scale_z;
	//most x
	*(vp++) = (float) scale_x/*prad*/ + (local_x * cos_alongz) - (sin_alongz * 0);
	*(vp++) = (float) scale_y + (0 * cos_alongz) + (sin_alongz * local_x);
	*(vp++) = (float) scale_z;
	//most z
	*(vp++) = (float) scale_x;
	*(vp++) = (float) scale_y;
	*(vp++) = (float) scale_zprad; 
	//least z
	*(vp++) = (float) scale_x;
	*(vp++) = (float) scale_y;
	*(vp++) = (float) scale_zmrad; 
	
	return vertices;
}

//Simulation setup
int main()
{
	//OpenGL initialization
	if (!glfwInit())
	{
		fprintf(stderr, "Failed GLFW\n");
		return 1;
	}
	glfwWindowHint(GLFW_SAMPLES, 4);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

	//Create a window
	GLFWwindow *window;
	window = glfwCreateWindow(WIDTH,HEIGHT, "Gravitas", NULL, NULL);

	if (window == nullptr){
		fprintf(stderr, "Failed window\n");
		glfwTerminate();
		return 1;
	}
	glfwMakeContextCurrent(window);

	glewExperimental = GL_TRUE;
	if (glewInit() != GLEW_OK){
		fprintf(stderr, "GLEW\n");
	return 1;
	}

	glViewport(0,0,WIDTH,HEIGHT);

	
	glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);

	//SHADERS
	unsigned int vertexShader;
	vertexShader = glCreateShader(GL_VERTEX_SHADER);

	glShaderSource(vertexShader, 1, &vertexShaderSource, NULL);
	glCompileShader(vertexShader);

	int success;
	char infoLog[512];
	glGetShaderiv(vertexShader, GL_COMPILE_STATUS, &success);

	unsigned int fragmentShader;
	fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);

	glShaderSource(fragmentShader, 1, &fragmentShaderSource, NULL);
	glCompileShader(fragmentShader);

	if (!success)
	{
		return 1;
	}

	unsigned int shaderProgram;
	shaderProgram = glCreateProgram();

	glAttachShader(shaderProgram, vertexShader);
	glAttachShader(shaderProgram, fragmentShader);
	glLinkProgram(shaderProgram);


	glDeleteShader(vertexShader);
	glDeleteShader(fragmentShader);

	//MASS INIT
	struct mass_list *root = new_node();
	struct mass_list *cp = new_massive_object(root, MASS_SUN);
	cp->m_obj.radius = R_SUN;
	cp->m_obj.rotation = ROT_SUN;
	cp->m_obj.rotation_cur = 0.0;

	cp->m_obj.loc_x = 0.0;
	cp->m_obj.loc_y = 0.0;
	cp->m_obj.loc_z = 0.0;

	struct vector_3 empty;
	empty.magnitude = 0.0;
	empty.x_comp = 0.0;
	empty.y_comp = 0.0;
	empty.z_comp = 0.0;

	cp->m_obj.vect = empty;

	struct mass_list *ep;
	cp = new_massive_object(root, MASS_EARTH);
	cp->m_obj.radius = R_EARTH;
	cp->m_obj.rotation = ROT_EARTH;
	cp->m_obj.rotation_cur = 0.0;
	
	cp->m_obj.loc_x = DISTANCE_ES;
	cp->m_obj.loc_y = 0.0;
	cp->m_obj.loc_z = 0.0;

	struct vector_3 earth_init;
	earth_init.magnitude = I_VEL_E;
	earth_init.x_comp = 0.0;
	earth_init.y_comp = I_VEL_E;
	earth_init.z_comp = 0.0;
	cp->m_obj.vect = earth_init;
	ep = cp;
	
	struct mass_list *mp;
	cp = new_massive_object(root, MASS_MOON);
	cp->m_obj.radius = R_MOON;
	cp->m_obj.rotation = ROT_MOON;
	cp->m_obj.rotation_cur = 0.0;
	cp->m_obj.loc_x = DISTANCE_ES;
	cp->m_obj.loc_y = DISTANCE_EM;
	cp->m_obj.loc_z = 0.0;

	struct vector_3 moon_init;
	moon_init.magnitude = I_VEL_M;
	moon_init.x_comp = I_VEL_M;
	moon_init.y_comp = 0.0;
	moon_init.z_comp = 0.0;
	moon_init = add_vectors(&earth_init, &moon_init);
	cp->m_obj.vect = moon_init;
	mp = cp;

	struct mass_list *jp;
	cp = new_massive_object(root, MASS_JUPITER);
	cp->m_obj.radius = R_JUPITER;
	cp->m_obj.rotation = ROT_JUPITER;
	cp->m_obj.rotation_cur = 0.0;
	cp->m_obj.loc_x = DISTANCE_JS;
	cp->m_obj.loc_y = 0.0;
	cp->m_obj.loc_z = 0.0;

	struct vector_3 jup_vel;
	jup_vel.magnitude = I_VEL_J;
	jup_vel.y_comp = I_VEL_J;
	jup_vel.x_comp = 0.0;
	jup_vel.z_comp = 0.0;
	cp->m_obj.vect = jup_vel;
	jp = cp;

	cp = new_massive_object(root, MASS_MARS);
	cp->m_obj.radius = R_MARS;
	cp->m_obj.rotation = ROT_MARS;
	cp->m_obj.rotation_cur = 0.0;
	cp->m_obj.loc_x = 0.0;
	cp->m_obj.loc_y = DISTANCE_MS * -1.0;
	cp->m_obj.loc_z = 0.0;
	cp->m_obj.vect.magnitude = I_VEL_MARS;
	cp->m_obj.vect.x_comp = I_VEL_MARS;
	cp->m_obj.vect.y_comp = 0.0;
	cp->m_obj.vect.z_comp = 0.0;
	struct mass_list *marsp = cp;
	
	float *drawMass;
	
	//A reference for size of drawMass and for layout of indices
	float testCoords[] = {
		-0.6f, 0.0f, 0.0f,
		0.0f, -0.6f, 0.0f,
		0.0f, 0.6f, 0.0f,
		0.6f, 0.0f, 0.0f,
		0.0f, 0.0f, 0.6f,
		0.0f, 0.0f, -0.6f

	};

	//To create a double pyramid shape
	unsigned int indices[] = {
		//0 left
		//1 bottom
		//2 top
		//3 right
		//4 front
		//5 back
		0, 2, 5, // left top back
		0, 4, 2, // top left front
		3, 2, 4, // top right front
		3, 5, 2, // top right back
		3, 1, 5, // bottom right back
		3, 1, 4, // bottom right front
		0, 1, 5, // bottom left back
		0, 1, 4 // bottom left front
	};

	//VERTEX ARRAYS, BUFFERS
	unsigned int VBO[2], VAO[2];
	unsigned int EBO[2];
	glGenBuffers(2, EBO);

	glGenVertexArrays(2, VAO);

	glGenBuffers(2, VBO);

	glBindVertexArray(VAO[0]);

	glBindBuffer(GL_ARRAY_BUFFER, VBO[0]);
//	glBufferData(GL_ARRAY_BUFFER, sizeof(float)*18, drawMass, GL_DYNAMIC_DRAW);

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO[0]);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices), indices, GL_STATIC_DRAW);

	//glGetIntegerv(GL_MAX_VERTEX_ATTRIBS, &addressOfAnInt);
	
	//0 layout vertex shader, 3 size of attrib vec3, type, normalized, space between vertex attrib, offset
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
	glEnableVertexAttribArray(0);

	//GL_FILL
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

	fprintf(stderr, "All init done\n");

	while (!glfwWindowShouldClose(window)){
		
		int i = 0;
		time_passed_total += TICK;
		for (cp = root; cp != nullptr; cp = cp->next){
			//Iterate through each object and get new vertex coordinates
			gravity_list_iterate(root, i);
			drawMass = genVertices(cp->m_obj);
			i++;
			
			//Pass new vertex values into buffer
			glBufferData(GL_ARRAY_BUFFER, sizeof(testCoords), drawMass, GL_DYNAMIC_DRAW);
			//can retrieve time
			//glfwGetTime();
			processInput(window);

			glClearColor(0.2f, 0.3f, 0.3f, 1.0f);

			glUseProgram(shaderProgram);
			glBindVertexArray(VAO[0]);
			//glDrawArrays(GL_TRIANGLES, 0, 3);
			glDrawElements(GL_TRIANGLES, 24, GL_UNSIGNED_INT, 0);

			glfwSwapBuffers(window);
			glfwPollEvents();

			//Don't do this until drawing is done
			//Removing this causes huge memory leak
			free(drawMass);
		}
		//Change the location of each object based on data from their vectors
		for (cp = root; cp != nullptr; cp = cp->next)
			update_location(cp);

		glClear(GL_COLOR_BUFFER_BIT);
	}

	glfwTerminate();

	//printf("Radius: %E\n", radius);

	return 0;
}
