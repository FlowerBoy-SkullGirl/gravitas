#include <stdio.h>
#include <math.h>
#include <ncurses.h>
#include <unistd.h>
#include <stdlib.h>

//Measurements in kilograms and meters
#define G_CONSTANT 6.674e-11
#define PI 3.1415
#define TICK 100000
#define TEST_MASS_1 10e11
#define TEST_MASS_2 10e11
#define TEST_MASS_3 10e9
#define SCALE 2e10
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

struct vector_3{
	double magnitude;
	double x_comp;
	double y_comp;
	double z_comp;
	double alpha;
	double beta;
	double gamma;
};

struct massive_object{
	double mass;
	struct vector_3 vect;
	double loc_x;
	double loc_y;
	double loc_z;
};

struct mass_list{
	struct massive_object m_obj;
	struct mass_list *next;
};

//Find the length of a vector given the components of each axis
double hypot_3d(struct massive_object mo1, struct massive_object mo2)
{
	double x_diff = fabs(mo2.loc_x - mo1.loc_x);
	double y_diff = fabs(mo2.loc_y - mo1.loc_y);
	double z_diff = fabs(mo2.loc_z - mo1.loc_z);
	double h_xy = sqrt((x_diff * x_diff) + (y_diff * y_diff));
	return sqrt((z_diff * z_diff) + (h_xy * h_xy));
}

void component_3d(struct vector_3 *vect)
{
	double sin_alpha = sin(vect->alpha);
	double sin_beta = sin(vect->beta);
	double cos_alpha = cos(vect->alpha);
	double cos_beta = cos(vect->beta);
	vect->x_comp = vect->magnitude * ((cos_alpha*sin_beta)/(sqrt((cos_beta*cos_beta)*(sin_alpha*sin_alpha)+(sin_beta*sin_beta))));
	vect->z_comp = tan(vect->alpha) * vect->x_comp;
	vect->y_comp = vect->z_comp / tan(vect->beta);
	
	//printf("Sina, sinb, cosa, cosb: %f, %f, %f, %f\n", sin_alpha, sin_beta, cos_alpha, cos_beta);
	//printf("Vector components: [%e, %e, %e]\n", vect->x_comp, vect->z_comp, vect->y_comp);
}

void component_3d_alt(struct vector_3 *vect)
{
	if (vect->gamma && vect->gamma > SIGNIFICANT_RADIANS){
		vect->y_comp = vect->magnitude * cos(vect->gamma);
		double hypot_xz = vect->magnitude * sin(vect->gamma);
		vect->x_comp = hypot_xz * cos(vect->alpha);
		vect->z_comp = hypot_xz * sin(vect->alpha);
	}else if (vect->alpha && vect->alpha > SIGNIFICANT_RADIANS){
		vect->x_comp = vect->magnitude * cos(vect->alpha);
		vect->z_comp = vect->magnitude * sin(vect->alpha);
	}else{
		vect->x_comp = vect->magnitude;
		vect->y_comp = 0.0;
		vect->z_comp = 0.0;
	}
}

//Find alpha and beta angles of a vector such that they corrsepsond to x-z angle and y-z angle
void find_angle_vector(struct vector_3 *dest, struct massive_object mo1, struct massive_object mo2)
{
	double x_diff = fabs(mo2.loc_x - mo1.loc_x);
	double y_diff = fabs(mo2.loc_y - mo1.loc_y);
	double z_diff = fabs(mo2.loc_z - mo1.loc_z);
	if (x_diff == 0.0)
		dest->alpha = 0.0;
	else
		dest->alpha = atan(z_diff/x_diff);
	if (y_diff == 0.0){
		dest->beta = 0.0;
		dest->gamma = 0.0;
	}else{
		dest->beta = atan(z_diff/y_diff);
		dest->gamma = atan(x_diff/y_diff);
	}
	//printf("Angles a, b, g: %f, %f, %f\n", dest->alpha, dest->beta, dest->gamma);
}

//find the force of gravity between two masses
double force_gravity(struct massive_object mo1, struct massive_object mo2)
{
	double radius = hypot_3d(mo1, mo2);
	if (mo1.mass == MASS_EARTH)
		mvprintw(5, 12, "Radius ES: %e", radius);
	//printf("Object radius: %f\n", radius);
	return G_CONSTANT * ((mo1.mass * mo2.mass)/(radius * radius));
}

//Find the acceleration due to gravity on mass 1 
double acceleration_due_to_gravity(struct massive_object mo1, struct massive_object mo2)
{
	//printf("Force of gravity between MOs: %EN\n", force_gravity(mo1, mo2));
	if (mo1.mass == MASS_EARTH && mo2.mass == MASS_SUN)
		mvprintw(7, 30, "FG: %e", force_gravity(mo1,mo2));
	return (force_gravity(mo1, mo2)/(mo1.mass));
	
}

/*
 * Find the vector created by one tick time of acceleration due to gravity for mo1
 * break it into components x, y, z
 * determine the sign of each component based on coordinates of mo1 and mo2
 */
struct vector_3 vector_acc_per_tick(struct massive_object mo1, struct massive_object mo2)
{
	struct vector_3 v_acc;
	v_acc.magnitude = acceleration_due_to_gravity(mo1, mo2) * TICK;
	if (mo1.mass == MASS_EARTH)
		mvprintw(9, 5, "E V_ACC: %e", v_acc.magnitude);
	find_angle_vector(&v_acc, mo1, mo2);
	component_3d_alt(&v_acc);
	
	if (mo1.loc_x > mo2.loc_x)
		v_acc.x_comp *= -1.0;
	if (mo1.loc_y > mo2.loc_y)
		v_acc.y_comp *= -1.0;
	if (mo1.loc_z > mo2.loc_z)
		v_acc.z_comp *= -1.0;
	if (mo1.mass == MASS_EARTH){
		mvprintw(15,5, "E V_ACC X, Y: %e, %e", v_acc.x_comp, v_acc.y_comp);
		mvprintw(13,5, "E V_ACC y, a, b: %e, %e, %e", v_acc.gamma, v_acc.alpha, v_acc.beta);
	}
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

struct mass_list *new_node()
{
	struct mass_list *np = (struct mass_list *) malloc(sizeof(struct mass_list));
	np->m_obj.mass = 0.0;
	return np;
}

//RETURNS pointer to new obj
struct mass_list *new_massive_object(struct mass_list *root, double d_mass)
{
	struct mass_list *np;
	if (root->m_obj.mass == 0.0){
		root->m_obj.mass = d_mass;
		return root;
	}
	if (root->next == NULL)
		root->next = new_node();

	np = new_massive_object(root->next, d_mass);
	return np;
	
}

//Takes an index to remove selected massive object
struct mass_list *remove_massive_object(struct mass_list *root, int index)
{
	struct mass_list *np = root;
	for (int i = 0; i < index - 1; i++){
		if (np->next != NULL)
			np = np->next;
	}
	if (np != NULL && np->next != NULL){
		np->next = np->next->next;
		free(np->next);
	}

	return root;
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
		if (np->next != NULL)
			np = np->next;
	}
	struct mass_list *tp = root;
	struct vector_3 v_change;
	for (int i = 0; tp !=NULL; i++){
		if (tp != np){
			v_change = vector_acc_per_tick((np->m_obj), (tp->m_obj));
			np->m_obj.vect = add_vectors(&np->m_obj.vect, &v_change);
		}
		tp = tp->next; 
	}
		
}

//Simulation setup
int main()
{

	//ncurses init
	initscr();
	int ncur_x = 0;
	int ncur_y = 0;
	getmaxyx(stdscr, ncur_y, ncur_x);
	ncur_x /= 2;
	ncur_y /= 2;

/*
	struct massive_object mo1, mo2, mo3;
	mo1.mass = TEST_MASS_1;
	mo2.mass = TEST_MASS_2;
	mo3.mass = TEST_MASS_3;
	mo1.loc_x = 10.0;
	mo1.loc_y = 0.0;
	mo1.loc_z = 0.0;
	
	mo2.loc_x = -10.0;
	mo2.loc_y = 0.0;
	mo2.loc_z = 0.0;
*/
	struct mass_list *root = new_node();
	struct mass_list *cp = new_massive_object(root, MASS_SUN);

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
	cp->m_obj.loc_x = 0.0;
	cp->m_obj.loc_y = DISTANCE_MS * -1.0;
	cp->m_obj.loc_z = 0.0;
	cp->m_obj.vect.magnitude = I_VEL_MARS;
	cp->m_obj.vect.x_comp = I_VEL_MARS;
	cp->m_obj.vect.y_comp = 0.0;
	cp->m_obj.vect.z_comp = 0.0;
	struct mass_list *marsp = cp;
	
	/*
	cp = new_massive_object(root, MASS_SUN / 25.0);

	cp->m_obj.loc_x = DISTANCE_ES * -6.0;
	cp->m_obj.loc_y = 0.0;
	cp->m_obj.loc_z = 0.0;

	struct vector_3 wowie;
	wowie.magnitude = 0.0;
	wowie.x_comp = I_VEL_MARS * 0.3;
	wowie.y_comp = I_VEL_J * -0.3; 
	wowie.z_comp = 0.0;

	cp->m_obj.vect = wowie;
	*/
/*
	struct vector_3 v_init_mo1, v_init_mo2, v_init_mo3;
	v_init_mo1.magnitude = 2.0;
	v_init_mo2.magnitude = 2.0;
	
	v_init_mo1.alpha = 0.0;
	v_init_mo1.beta = 0.0;
	v_init_mo2.alpha = 0.0;
	v_init_mo2.beta = 0.0;

	v_init_mo3.alpha = 0.0;
	v_init_mo3.beta = 0.0;
	puts("Comp vi1");
	component_3d(&v_init_mo1);
	puts("Comp vi2");
	component_3d(&v_init_mo2);
	v_init_mo1.x_comp = 1.0;
	v_init_mo1.y_comp = 1.0;
	v_init_mo1.z_comp = 0.0;
	v_init_mo2.x_comp = 1.0;
	v_init_mo2.y_comp = 1.0;
	v_init_mo2.z_comp = 0.0;
	v_init_mo2.y_comp *= -1.0;
	v_init_mo2.x_comp *= -1.0;

	v_init_mo3.magnitude = 0.0;
	v_init_mo3.x_comp = 0.0;
	v_init_mo3.y_comp = 0.0;
	v_init_mo3.z_comp = 0.0;
	
	mo1.vect = v_init_mo1;
	mo2.vect = v_init_mo2;
	mo3.vect = v_init_mo3;
*/	
	long coord_x, coord_y, coord_z = 0;

	erase();

	char SYMBOLS[] = {'S', 'E', 'm', 'J', 'M', '$'};
	while (true){
		
		int i = 0;
		time_passed_total += TICK;
		for (cp = root; cp != NULL; cp = cp->next){
			coord_x = cp->m_obj.loc_x / SCALE;
			coord_y = cp->m_obj.loc_y / SCALE;
			coord_z = cp->m_obj.loc_z / SCALE;
			
			gravity_list_iterate(root, i);
			mvprintw(coord_y + ncur_y, coord_x + ncur_x, "%c", SYMBOLS[i]);

			if (cp->m_obj.mass == MASS_EARTH){
				mvprintw(4, ncur_x + 20, "Current v: [%e,%e,%e]", root->next->m_obj.vect.x_comp, root->next->m_obj.vect.y_comp, root->next->m_obj.vect.z_comp);
				mvprintw(6, ncur_x + 20, "Mass: %e", root->next->m_obj.mass);
				mvprintw(8, ncur_x + 20, "Pos Earth: [%e, %e, %e]", root->next->m_obj.loc_x, root->next->m_obj.loc_y, root->next->m_obj.loc_z); 
				mvprintw(9, ncur_x + 20, "Time s: %e", time_passed_total);

			}
		i++;
		}
		for (cp = root; cp != NULL; cp = cp->next)
			update_location(cp);

		refresh();
		usleep(T_SLEEP);
		erase();
	}

	endwin();

	//printf("Radius: %E\n", radius);

	return 0;
}
