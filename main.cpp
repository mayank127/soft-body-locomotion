#include "corotationalLinearFEM.h"
#include "volumetricMeshLoader.h"
#include "tetMesh.h"
#include "vec3d.h"
#include "sparseMatrix.h"
#include "muscle.h"
#include "generateMassMatrix.h"
#include "renderVolumetricMesh.h"

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <string>
#include <GL/glut.h>
#include <vector>
using namespace std;

VolumetricMesh * volumetricMesh = NULL;
TetMesh * tetMesh = NULL;
vector<muscle> muscles;
vector< vector<double> > elongations;
CorotationalLinearFEM * corotationalLinearFEM;
RenderVolumetricMesh renderVolumetricMesh;

SparseMatrix * stiffnessMatrix;
SparseMatrix * massMatrix_sparse;
SparseMatrix * Mtilda_sparse;

double * velocity;
double * external_force;
double * massMatrix;
double * internal_force;
double * displacement;
double * muscle_force;
double * Mtilda;
double * force_calculation;
int numVertices;


int win_width = 512;
int win_height = 512;

void drawTriangle(Vec3d a, Vec3d b, Vec3d c){
	Vec3d normal = norm(cross(b-a, c-b));
	glBegin(GL_POLYGON);
		glNormal3f(normal[0], normal[1], normal[2]);
		glVertex3f(a[0],a[1],a[2]);
		glVertex3f(b[0],b[1],b[2]);
		glVertex3f(c[0],c[1],c[2]);
	glEnd();
}

void draw() {
	int numEle = tetMesh->getNumElements();
	/*for(int i=0;i<numEle; i++){
		Vec3d & a = *(tetMesh->getVertex(i, 0));
		Vec3d & b = *(tetMesh->getVertex(i, 1));
		Vec3d & c = *(tetMesh->getVertex(i, 2));
		Vec3d & d = *(tetMesh->getVertex(i, 3));
		glColor3f(0.0, 1.0, 0.0);
		drawTriangle(a,b,c);
		//glColor3f(1.0, 1.0, 0.0);
		drawTriangle(a,c,d);
		//glColor3f(0.0, 1.0, 1.0);
		drawTriangle(a,d,b);
		//glColor3f(1.0, 0.0, 1.0);
		drawTriangle(b,d,c);
	}*/
	for (int i=0;i<muscles.size();i++) {
		for (int j=0;j<muscles[i].muscle_segments.size();j++) {
			Vec3d & p = *(tetMesh->getVertex(muscles[i].muscle_segments[j].centerIndex));
			Vec3d dir = muscles[i].muscle_segments[j].direction;
			dir = dir*muscles[i].muscle_segments[j].cur_length / 2;
			Vec3d p1 = p + dir;
			Vec3d p2 = p - dir;
			glBegin(GL_LINE);
				glColor3f(1.0, 0.0, 0.0);
				glVertex3f(p1[0], p1[1], p1[2]);
				glVertex3f(p2[0], p2[1], p2[2]);
			glEnd();

		}
	}
	glBegin(GL_POLYGON);
		// glPointSize(8);
		glColor3f(1.0, 0.0, 0.0);
		glVertex3f(-100, -0.25, 0);
		glVertex3f(-100, -0.5, 0);
		glVertex3f(100, -0.5, 0);
		glVertex3f(100, -0.25, 0);
		glVertex3f(-100, -0.25, 0);
	glEnd();
}

void display( void )
{
	glClear(GL_COLOR_BUFFER_BIT);
	glPushMatrix();
		glScalef(3,3,3);
		draw();
		renderVolumetricMesh.Render(tetMesh);
	glPopMatrix();
	glutSwapBuffers();
}

void reshape( int w, int h )
{
	if  ( h == 0 ) h = 1;
	glMatrixMode( GL_PROJECTION );
	glLoadIdentity();
	glOrtho( -4.0, 4.0, -4.0, 4.0, -1., 1. );
	glViewport( 0, 0, w, h );
	win_width = w;
	win_height = h;
	glutPostRedisplay();
}

Vec3d multiply_vec(Vec3d vec, double * mat){
	Vec3d result(0.0);
	for(int i=0; i < 3; i++){
		for(int j=0;j<3;j++){
			result[j] += vec[i] * mat[3*i + j];
		}
	}
	return result;
}


void get_muscle_force(double * displacement, double * muscle_f){
	int numVertices = tetMesh->getNumVertices();
	for(int i=0; i< tetMesh->getNumVertices() * 3; i++)
		muscle_f[i] = 0;

	for(int i = 0; i < muscles.size(); i++){
		muscles[i].force(displacement, elongations[i]);
	}
	for(int el = 0; el < tetMesh->getNumElements();el++){
		vector<double> radWeight, longWeight;
		for(int i = 0; i < muscles.size(); i++){
			if(muscles[i].type == 0){
				longWeight.push_back(muscles[i].get_distance(el));
			}
			else{
				radWeight.push_back(muscles[i].get_distance(el));
			}
		}
		double lSum = 0, rSum = 0;
		for(int i =0; i < longWeight.size(); i++){
			lSum += 1/longWeight[i];
		}
		for(int i =0; i < radWeight.size(); i++){
			rSum += 1/radWeight[i];
		}
		int j=0,k=0;
		double force[9];
		for(int l=0; l< 9; l++){
			force[l] = 0;
		}

		for(int i = 0; i < muscles.size(); i++){
			if(muscles[i].type == 0){
				double d = 1 / (longWeight[j] * lSum);
				for(int l=0; l< 9; l++){
					force[l] += d * muscles[i].Force[l];
				}
				j++;
			}
			else{
				double d = 1 / (longWeight[k] * lSum);
				for(int l=0; l< 9; l++){
					force[l] += d * muscles[i].Force[l];
				}
				k++;
			}
		}

		int vtxIndex[4];
		for (int vtx=0; vtx<4; vtx++)
			vtxIndex[vtx] = tetMesh->getVertexIndex(el, vtx);
		Vec3d & a = *(tetMesh->getVertex(vtxIndex[0]));
		Vec3d & b = *(tetMesh->getVertex(vtxIndex[1]));
		Vec3d & c = *(tetMesh->getVertex(vtxIndex[2]));
		Vec3d & d = *(tetMesh->getVertex(vtxIndex[3]));
		Vec3d normal1 = cross(b-a, c-b);
		Vec3d normal2 = cross(c-a, d-c);
		Vec3d normal3 = cross(d-a, b-d);
		Vec3d normal4 = cross(d-b, c-d);
		double sum = len(normal1) + len(normal2) + len(normal3) + len(normal4);

		Vec3d face_force = multiply_vec((normal1 + normal2 +  normal3), force);
		muscle_f[3 * vtxIndex[0] + 0] += 1/(3*sum) * face_force[0];
		muscle_f[3 * vtxIndex[0] + 1] += 1/(3*sum) * face_force[1];
		muscle_f[3 * vtxIndex[0] + 2] += 1/(3*sum) * face_force[2];

		face_force = multiply_vec((normal1 + normal4 +  normal3), force);
		muscle_f[3 * vtxIndex[1] + 0] += 1/(3*sum) * face_force[0];
		muscle_f[3 * vtxIndex[1] + 1] += 1/(3*sum) * face_force[1];
		muscle_f[3 * vtxIndex[1] + 2] += 1/(3*sum) * face_force[2];

		face_force = multiply_vec((normal1 + normal2 +  normal4), force);
		muscle_f[3 * vtxIndex[2] + 0] += 1/(3*sum) * face_force[0];
		muscle_f[3 * vtxIndex[2] + 1] += 1/(3*sum) * face_force[1];
		muscle_f[3 * vtxIndex[2] + 2] += 1/(3*sum) * face_force[2];

		face_force = multiply_vec((normal4 + normal2 +  normal3), force);
		muscle_f[3 * vtxIndex[3] + 0] += 1/(3*sum) * face_force[0];
		muscle_f[3 * vtxIndex[3] + 1] += 1/(3*sum) * face_force[1];
		muscle_f[3 * vtxIndex[3] + 2] += 1/(3*sum) * face_force[2];

	}
}

void invertMatrix(double* data, int actualsize, int maxsize)  {
    if (actualsize <= 0) return;  // sanity check
    if (actualsize == 1) return;  // must be of dimension >= 2
    for (int i=1; i < actualsize; i++) data[i] /= data[0]; // normalize row 0
    for (int i=1; i < actualsize; i++)  { 
      for (int j=i; j < actualsize; j++)  { // do a column of L
        double sum = 0.0;
        for (int k = 0; k < i; k++)  
            sum += data[j*maxsize+k] * data[k*maxsize+i];
        data[j*maxsize+i] -= sum;
        }
      if (i == actualsize-1) continue;
      for (int j=i+1; j < actualsize; j++)  {  // do a row of U
        double sum = 0.0;
        for (int k = 0; k < i; k++)
            sum += data[i*maxsize+k]*data[k*maxsize+j];
        data[i*maxsize+j] = 
           (data[i*maxsize+j]-sum) / data[i*maxsize+i];
        }
      }
    for ( int i = 0; i < actualsize; i++ )  // invert L
      for ( int j = i; j < actualsize; j++ )  {
        double x = 1.0;
        if ( i != j ) {
          x = 0.0;
          for ( int k = i; k < j; k++ ) 
              x -= data[j*maxsize+k]*data[k*maxsize+i];
          }
        data[j*maxsize+i] = x / data[j*maxsize+j];
        }
    for ( int i = 0; i < actualsize; i++ )   // invert U
      for ( int j = i; j < actualsize; j++ )  {
        if ( i == j ) continue;
        double sum = 0.0;
        for ( int k = i; k < j; k++ )
            sum += data[k*maxsize+j]*( (i==k) ? 1.0 : data[i*maxsize+k] );
        data[i*maxsize+j] = -sum;
        }
    for ( int i = 0; i < actualsize; i++ )   // final inversion
      for ( int j = 0; j < actualsize; j++ )  {
        double sum = 0.0;
        for ( int k = ((i>j)?i:j); k < actualsize; k++ )  
            sum += ((j==k)?1.0:data[j*maxsize+k])*data[k*maxsize+i];
        data[j*maxsize+i] = sum;
        }
    };

void multiply_matrix(double * a, double * b, double * c){
	for(int i=0; i< 3 * numVertices; i++){
		c[i] = 0;
		for(int j=0; j< 3 * numVertices; j++){
			c[i] += a[i*3*numVertices + j] * b[j];
		}
	}
}

void calculate_elongations() {
	for (int i=0;i<muscles.size();i++) {
		for (int j=0;j<muscles[i].muscle_segments.size();j++) {
			double f[3] = {external_force[muscles[i].muscle_segments[j].centerIndex*3], external_force[muscles[i].muscle_segments[j].centerIndex*3 + 1], external_force[muscles[i].muscle_segments[j].centerIndex*3 + 2]};
			Vec3d dir = muscles[i].muscle_segments[j].direction;
			double dot = dir[0]*f[0] + dir[1]*f[1] + dir[2]*f[2];
			double len = dot / muscles[i].muscle_segments[j].stiffness;
			if (abs(len) > 0.5*muscles[i].muscle_segments[j].init_length) {
				len = 0.5*muscles[i].muscle_segments[j].init_length * len/abs(len);
			}
			elongations[i][j] = len + muscles[i].muscle_segments[j].init_length;
			cout << elongations[i][j] << "  ";
		}
		cout << endl;
	}
	cout << endl << endl;
}

void genetrate_mass(){
	GenerateMassMatrix::computeMassMatrix(tetMesh, &massMatrix_sparse, true);
	massMatrix = new double[3 * tetMesh->getNumVertices() * 3 * tetMesh->getNumVertices()];
	massMatrix_sparse->MakeDenseMatrix(massMatrix);
}


void update_velocity(int t){
	double dt = 0.01;
	for(int i=0; i< numVertices * 3; i++){
		displacement[i] = dt * velocity[i];
	}
	tetMesh->applyDeformation(displacement);
	corotationalLinearFEM->ComputeForceAndStiffnessMatrix(displacement, internal_force, stiffnessMatrix, 1);
	calculate_elongations();
	get_muscle_force(displacement, muscle_force);

	stiffnessMatrix->MakeDenseMatrix(Mtilda);
	for(int i=0; i<3*numVertices;i++){
		for(int j=0; j<3*numVertices;j++){
			Mtilda[i*3*numVertices + j] *= dt * 0.2 + dt*dt;
			Mtilda[i*3*numVertices + j] += massMatrix[i*3*numVertices + j];
		}
	}
	invertMatrix(Mtilda, 3*numVertices, 3*numVertices);


	multiply_matrix(massMatrix, velocity, force_calculation);
	for(int i = 0; i< 3 *numVertices; i++){
		force_calculation[i] += dt * (external_force[i] - internal_force[i] + muscle_force[i]);
	}
	multiply_matrix(Mtilda, force_calculation, velocity);
	for(int i=0;i<numVertices; i++){
		Vec3d & p = *(tetMesh->getVertex(i));
		if (p[1]<-0.2) {
			for (int j=0;j<4;j++) {
				external_force[i*3 + 1] = -1 * force_calculation[i*3 + 1];
				velocity[i*3+1] = 0;
			}
		}
	}
	// for(int i = 0; i< 3 *numVertices; i++){
	// 	// cout<<velocity[i]<<" ";
	// }
	cout<<endl;
	cout<<"HERE"<<endl;
	glutPostRedisplay();
	glutTimerFunc(1, update_velocity, 1);
}


int main (int argc, char *argv[]) 
{

	glutInit( &argc, argv );
	glutInitDisplayMode( GLUT_RGBA | GLUT_DOUBLE );
	glutInitWindowSize( win_width, win_height );

	glutCreateWindow( "Soft Body Locomotion" );

	glClearColor(0.0f, 0.0f, 0.0f, 1.0f);


	char * volumetricMeshFilename = "examples/beam3_tet.veg";
	volumetricMesh = VolumetricMeshLoader::load(volumetricMeshFilename);
	tetMesh = dynamic_cast<TetMesh*>(volumetricMesh);
	if (tetMesh == NULL)
	{
		printf("Error: the input mesh is not a tet mesh (CLFEM deformable model).\n");
		exit(1);
	}

	numVertices = tetMesh->getNumVertices();
	corotationalLinearFEM = new CorotationalLinearFEM(tetMesh);
	genetrate_mass();
	velocity = new double[3 * numVertices];
	external_force = new double[3 * numVertices];
	internal_force = new double[3 * numVertices];
	muscle_force= new double[3 * numVertices];
	displacement = new double[3 * numVertices];
	Mtilda = new double[3 * numVertices * 3 * numVertices];
	force_calculation = new double[3 * numVertices];
	for(int i=0;i<3*numVertices; i++){
		velocity[i] = 0;
		external_force[i] = 0;
		internal_force[i] = 0;
		displacement[i] = 0;
		muscle_force[i] = 0;
		force_calculation[i] = 0;
	}
	corotationalLinearFEM->GetStiffnessMatrixTopology(&stiffnessMatrix);

	for(int i=0;i<3*numVertices;i+=3){
		external_force[i+1] = -0.1;
	}
	external_force[200*3+0] = 10;
	// external_force[207 * 3 + 1] = -100;
	// for(int i=0;i<3*numVertices/2; i+=3){
	// 	external_force[i] = -1;
	// }


	muscle_segment ms1(15, 0.5, Vec3d(0,1,0), 5);
	muscle m(tetMesh, LONGITUDINAL);
	m.muscle_segments.push_back(ms1);

	muscle_segment ms2(10, 0.5, Vec3d(0,1,0), 2);
	muscle m2(tetMesh, LONGITUDINAL);
	m2.muscle_segments.push_back(ms2);

	muscles.push_back(m);
	muscles.push_back(m2);
	elongations.push_back(vector<double>(1,0.5));
	elongations.push_back(vector<double>(1,0.5));

	renderVolumetricMesh = RenderVolumetricMesh();

	glutTimerFunc(1, update_velocity, 1);
	glutDisplayFunc( display );
	glutReshapeFunc( reshape );
	glutMainLoop();
}
