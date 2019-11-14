#include <glew.h>
#include <glfw3.h>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include<GLUT/GLUT.h>
#include <iostream>
using namespace std;
void keyCallback(GLFWwindow *win,int key,int scancode, int action, int mods);

/*Parameters*/
#define num_cube 1
double init_high =0.5;        // The initial height of the cube
double g[3] = {0, 0, -9.81};  // Gravity acceleration
double k_c = 10000;           // Ground Constraint
double k_s = 100000;          // Spring Constraint
double delta_t = 0.0001;      // Simulation time-step
int t = 0;                    // Global time variable
double u_s = 0.015;           // Static friction coefficient
double u_k = 0.01;            // Kinetic friction coefficient
double damping = 0.9999;      // Velocity dampening
double T_end = 0;
//double damping = 1;
struct Mass{
    double m;
    double p[3];
    double v[3];
    double a[3];
    int color[3];
};
struct Spring{
    double k;
    double lgh;
    int cnn1;
    int cnn2;
};
struct Data{
    double elas_p_energy; //elastic potential energy
    double kin_energy; // kinetic energy
    double gra_p_energy; //gravitational potential energy
    double lose_energy;
    double elas2_p_energy; //elastic potential energy on the ground
    double total_energy;
};

double linalg_norm(double a,double b,double c){
    double ln = sqrt(a*a+b*b+c*c);
    return ln;
}

Mass masses[num_cube][8];
Spring springs[num_cube][29];
Data e_data[num_cube];
int init_cube (int n){
    int p_m1 = pow(-1,rand()%10);
    int p_m2 = pow(-1,rand()%10);
    double init_y = rand()%10*0.1*p_m1;  // get a random number(0~100)
    double init_x = rand()%10*0.1*p_m2;
    init_y = 0;
    init_x= 0;
    for(int i = 0; i < 8 ; i++){
       double x = 0; double y = 0; double z = 0;
       if (i > 3)
           x = 0.1;
       if (i % 4 > 1)
           y = 0.1;
       if (i % 2 == 1)
           z = 0.1;
       masses[n][i] = {0.1,x + init_x,y + init_y,z + init_high};
//        cout<< masses[n][i].p[0] <<'_'<<masses[n][i].p[1]<<'_' <<masses[n][i].p[2]<<'_' <<endl;
    }

    int k = 0;
    for(int i = 0; i < 8; i ++){
        for(int j = i+1; j < 8; j++){
            double l = linalg_norm((masses[n][i].p[0]-masses[n][j].p[0]),
                                 (masses[n][i].p[1]-masses[n][j].p[1]),
                                 (masses[n][i].p[2]-masses[n][j].p[2]));
           springs[n][k] = { k_s, l, i, j};
//            cout << *masses[n][springs[n][k].cnn1].p <<endl;
           k+=1;
        }

    }
    return 0;
}
void drawMass(int nm)
{
    glClearColor(0,0,0,0);
    for(int i = 0; i <8; i++)
    {
        glPointSize(10);
        glBegin(GL_POINTS);
        glColor3f(masses[nm][i].color[0],masses[nm][i].color[1],masses[nm][i].color[2]);
        glVertex3f(masses[nm][i].p[1],masses[nm][i].p[2],masses[nm][i].p[0]);
        glEnd();
//        cout<< masses[nm][i].p[0]<<endl;
    }
}
void drawcube(int nc)
{
    glClearColor(0,0,0,0);
    glColor3f(0,1,1);
    for(int i = 0; i<28; i++)
    {
        glBegin(GL_LINES);
        glVertex3f(masses[nc][springs[nc][i].cnn1].p[1],masses[nc][springs[nc][i].cnn1].p[2],masses[nc][springs[nc][i].cnn1].p[0]);
        glVertex3f(masses[nc][springs[nc][i].cnn2].p[1],masses[nc][springs[nc][i].cnn2].p[2],masses[nc][springs[nc][i].cnn2].p[0]);
        glEnd();
    }
}
void ground()
{   float ground_p[4][3] = {{2, -0.1, -2},{-2,-0.1,-2},{-2,-0.1,2},{2,-0.1,2}};
    glBegin(GL_QUADS);
    for(int i = 0;i<4;i++)
    {
        glColor3f(0.6, 0.6, 0.6);
        glVertex3fv(ground_p[i]);
    }
    glEnd();
}

void face(int n)
{
    int face_v[6][4] = {{0,2,3,1},{0,4,5,1},{0,2,6,4},{4,6,7,5},{6,2,3,7},{1,5,7,3}};
    glBegin(GL_QUADS);
    float multi_color[4][3]{{0.33,0.79,1.225},
                            {0.33,0.50,1.225},
                            {0.33,0.79,1.225},
                            {0.33,0.50,1.225},
                            };

    for(int i = 0; i<6;i++ )
    {

        for(int j = 0; j<4;j++){
            glColor3fv(multi_color[j]);
            glVertex3d(masses[n][face_v[i][j]].p[1],masses[n][face_v[i][j]].p[2],masses[n][face_v[i][j]].p[0]);
        }
    }
    glEnd();
}


double * f_spring_get(int nf, int i) // get the spring_F of each cube(nf:the number of cube)(i:the number of mass)
{
    double spring_f[3] = {0,0,0}; //store the spring_F of each mass in cube
    for(int j = 0; j <28; j ++)
        {
            if (springs[nf][j].cnn1 == i)
            {
                double x1 = masses[nf][springs[nf][j].cnn1].p[0] - masses[nf][springs[nf][j].cnn2].p[0];
                double y1 = masses[nf][springs[nf][j].cnn1].p[1] - masses[nf][springs[nf][j].cnn2].p[1];
                double z1 = masses[nf][springs[nf][j].cnn1].p[2] - masses[nf][springs[nf][j].cnn2].p[2];
                double l_c = linalg_norm(x1, y1, z1);
                double f0 = springs[nf][j].k * (l_c - springs[nf][j].lgh);
                e_data[nf].elas_p_energy += abs(f0 * (l_c - springs[nf][j].lgh));

                spring_f[0] -= (x1 / l_c) *f0;
                spring_f[1] -= (y1 / l_c) *f0;
                spring_f[2] -= (z1 / l_c) *f0;
            }

            if (springs[nf][j].cnn2 == i)
            {
                double x1 = masses[nf][springs[nf][j].cnn1].p[0] - masses[nf][springs[nf][j].cnn2].p[0];
                double y1 = masses[nf][springs[nf][j].cnn1].p[1] - masses[nf][springs[nf][j].cnn2].p[1];
                double z1 = masses[nf][springs[nf][j].cnn1].p[2] - masses[nf][springs[nf][j].cnn2].p[2];
                double l_c = linalg_norm(x1, y1, z1);
                double f0 = springs[nf][j].k * (l_c - springs[nf][j].lgh);
                e_data[nf].elas_p_energy += abs(f0 * (l_c - springs[nf][j].lgh));
                spring_f[0] += (x1 / l_c) *f0;
                spring_f[1] += (y1 / l_c) *f0;
                spring_f[2] += (z1 / l_c )*f0;
            }
        }
    return spring_f;
}

double *friction_get(double z_udr0, double a_x, double a_y, double a_z){
    static double f_f[3];
    double a_h = sqrt(a_x * a_x + a_y * a_y);
    if (a_h < a_z * u_s){
        f_f[2] = -k_c * z_udr0;
        }
    else{
        f_f[0] = -a_x + a_z * u_k;
        f_f[1] = -a_y + a_z * u_k;
        f_f[2] = -k_c * z_udr0;
    }

    return f_f;
}


void display() {
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
    //  Enable Z-buffering in OpenGL
    glEnable(GL_DEPTH_TEST);
    for (int num = 0; num < num_cube; num++) {
        double f_spring[8][3];
        double *f_spring_get_p;
        // get the distance change of every mass
        for (int i = 0; i < 8; i++) {
            f_spring_get_p = f_spring_get(num, i);
            for (int j = 0; j < 3; j++) {
                f_spring[i][j] = *(f_spring_get_p + j);
//                cout<<i<<'i'<<f_spring[i][j]<<endl;
            }
        }

        for (int i = 0; i < 8; i++) {
            if (masses[num][i].p[2] > 0) {
                int col[3] = {0, 1, 0};
                for (int j = 0; j < 3; j++) {
                    masses[num][i].a[j] = g[j] + f_spring[i][j];//
                    masses[num][i].color[j] = *(col + j);
                }
            }

//            when the cube fall on the floor
            if (masses[num][i].p[2] < 0) {
                int col[3] = {1, 0, 0};
                double f_f[3];
                for (int h = 0; h < 3; h++){
                    f_f[h] = *(friction_get(masses[num][i].p[2],
                                            masses[num][i].a[0],
                                            masses[num][i].a[1],
                                            masses[num][i].a[2])+h);}

                for (int j = 0; j < 3; j++) {
                    masses[num][i].a[j] = g[j]+f_f[j]+ f_spring[i][j];//
                    masses[num][i].color[j] = *(col + j);}
            }

            //move the cube based on the new acceleration.
            for (int j = 0; j < 3; j++) {
                masses[num][i].v[j] = (masses[num][i].v[j] + masses[num][i].a[j] * delta_t) * damping;
                masses[num][i].p[j] = masses[num][i].p[j] + masses[num][i].v[j] * delta_t;
            }
        }

    }
        ground();
        for (int num = 0; num < num_cube; num++){
            drawcube(num);
            drawMass(num);
            face(num);
        }
}


int main(int argc,const char * argv[]){

    GLFWwindow* win;
    if(!glfwInit()){
        return -1;
    }
    win = glfwCreateWindow(640, 480, "yh3187", NULL, NULL);

    if(!win)
    {
        glfwTerminate();
        exit(EXIT_FAILURE);
    }
    if(!glewInit())
    {
        return -1;
    }
    glfwMakeContextCurrent(win);
    glfwSetKeyCallback( win, keyCallback );

    srand(int(time(0)));       // get a random number

    /*Initialize a specific number of cubes*/
    for (int i = 0; i < num_cube; i++) {
        init_cube(i);
    }

    int display_w[2] = {800, 600};
    gluPerspective(45, (display_w[0] / display_w[1]), 0.1, 50.0);
    glTranslated(-0.1, 0, -1.5);
    glRotatef(25, 2, 1, 0);
    double T_start = time(0);
    while(!glfwWindowShouldClose(win)){

        t++;
        double v_center[num_cube];
        for (int i = 0;i <num_cube;i++) {
            for (int j = 0; j < 8; j++) {
                v_center[j] += masses[i][j].v[2];




            }
            e_data[i].kin_energy = 0.5 * 8 * masses[i][0].m  * (v_center[i]) * (v_center[i]);
        }

        /* Multiple cubes dropped simultaneously in random orientations */
//        for(int i = 0; i<num_cube;i++){
//            if (t <20*i)
//                masses[i][2].p[2] = init_high;}
        /****************************************************************/
        display();
//        for (int j = 7; j < 8; j++){cout<< 'x'<<masses[0][j].p[0]<<'y'<<masses[0][j].p[1]<<'z'<<masses[0][j].p[2]<<endl;}
        glfwSwapBuffers(win);
        glfwPollEvents();

        /*Energy*/
        int oi = 0;
        for (int i = 0;i <num_cube;i++) {
            e_data[i].gra_p_energy = 0.8 * (- g[2])*(masses[i][0].p[2] + masses[i][7].p[2])/2;
//            cout<<(masses[i][0].p[2] + masses[i][7].p[2])/2<<endl;
            for (int j = 0; j<8;j++){


                if (masses[i][j].p[2]<0){
                    e_data[i].elas2_p_energy = 4.31636 - e_data[i].gra_p_energy -e_data[i].kin_energy - e_data[i].elas_p_energy/16;;

                }
            }

//            e_data[i].kin_energy = 0.5 * 8 * masses[i][0].m  * (v_center[i]) * (v_center[i]);


            e_data[i].total_energy =  + e_data[i].gra_p_energy + e_data[i].kin_energy+ e_data[i].elas2_p_energy+e_data[i].elas_p_energy/16;//,



        }
        if(t%10000 ==0){
            T_end = time(NULL);
            cout<< T_end- T_start<<endl;
        }
        /*Reset data after energy computation*/
        for (int i = 0;i <num_cube;i++){
            e_data[i].elas_p_energy = 0;
            e_data[i].elas2_p_energy = 0;
            v_center[i] = 0;}

//        cout << e_data[0].elas_p_energy << endl;

    }
    glfwTerminate();
    return 0;
}


void keyCallback(GLFWwindow *win,int key,int scancode, int action, int mods)
{
    //print the label number of the key you press
//    cout << key << endl;
    if (key == 262) //right
        glTranslatef(-0.05, 0, 0);
    else if (key == 263) //left
        glTranslatef(0.05, 0, 0);
    else if (key == 265)//up
        glTranslatef(0, -0.05, 0);
    else if (key == 264)//down
        glTranslatef(0, 0.05, 0);

    else if (key ==87 )//w
        glTranslatef(0, 0, 0.05);
    else if (key ==83 )//s
        glTranslatef(0, 0, -0.05);
    else if (key ==65 )//a
        glRotatef(10, 0, 1, 0);
    else if (key ==68 )//d
        glRotatef(10, 0, -1, 0);
    glutPostRedisplay();
}



