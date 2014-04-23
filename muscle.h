#ifndef MUSCLE
#define MUSCLE

#include <vector>
#include "vec3d.h"
#include "tetMesh.h"
using namespace std;

enum MUSCLE_TYPE
{
	LONGITUDINAL, RADIAL
};

class muscle_segment
{
public:
	int centerIndex;
	double init_length;
	double cur_length;
	Vec3d direction;
	double stiffness;
	muscle_segment(int centerIndex, double init_length, Vec3d direction, double stiffness);
	void force(double * f, double len);
	~muscle_segment();
};

class muscle
{
public:
	int numVertices, numElements;
	TetMesh * tetMesh;
	double * undeformedPositions;
	double ** MInverse;
	vector<muscle_segment> muscle_segments;
	double Force[9];
	MUSCLE_TYPE type;

	muscle(TetMesh * tetMesh_, MUSCLE_TYPE type);
	void inverse4x4(double * A, double * AInv);
	void force(double * u, vector<double>);
	double get_distance(int el);
	~muscle();
};


#endif
