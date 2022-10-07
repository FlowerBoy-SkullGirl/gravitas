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
	double radius;
	double rotation;
	double rotation_cur;
};

struct mass_list{
	struct massive_object m_obj;
	struct mass_list *next;
};


struct mass_list *new_node()
{
	struct mass_list *np = (struct mass_list *) malloc(sizeof(struct mass_list));
	np->m_obj.mass = 0.0;
	np->next = nullptr;
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
	if (root->next == nullptr)
		root->next = new_node();

	np = new_massive_object(root->next, d_mass);
	return np;
	
}

//Takes an index to remove selected massive object
struct mass_list *remove_massive_object(struct mass_list *root, int index)
{
	struct mass_list *np = root;
	for (int i = 0; i < index - 1; i++){
		if (np->next != nullptr)
			np = np->next;
	}
	if (np != nullptr && np->next != nullptr){
		np->next = np->next->next;
		free(np->next);
	}

	return root;
}

//Find the length of a vector given the components of each axis
double hypot_3d(struct massive_object mo1, struct massive_object mo2)
{
	double x_diff = fabs(mo2.loc_x - mo1.loc_x);
	double y_diff = fabs(mo2.loc_y - mo1.loc_y);
	double z_diff = fabs(mo2.loc_z - mo1.loc_z);
	double h_xy = sqrt((x_diff * x_diff) + (y_diff * y_diff));
	return sqrt((z_diff * z_diff) + (h_xy * h_xy));
}

//NO LONGER USED: KEPT IN CASE OF USEFULNESS
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


//Breaks a vector into 3 components, x, y, and z based on given angle information
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
	//printf("Object radius: %f\n", radius);
	return G_CONSTANT * ((mo1.mass * mo2.mass)/(radius * radius));
}

//Find the acceleration due to gravity on mass 1 
double acceleration_due_to_gravity(struct massive_object mo1, struct massive_object mo2)
{
	//printf("Force of gravity between MOs: %EN\n", force_gravity(mo1, mo2));
	return (force_gravity(mo1, mo2)/(mo1.mass));
	
}


