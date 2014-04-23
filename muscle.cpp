#include "muscle.h"
#include "polarDecomposition.h"
#include <iostream>
using namespace std;


void multiply(double * A, double * B, double * AB){
	for(int i=0; i < 3; i++){
		for(int j = 0; j < 3; j++){
			AB[i*3 + j] = 0;
			for(int k = 0; k< 3; k++){
				AB[i*3 + j] += A[3*i + k] * B[3*k + j];
			}
		}
	}
}


muscle::muscle(TetMesh * tetMesh_ , MUSCLE_TYPE t)  : tetMesh(tetMesh_){
	type = t;
	numVertices = tetMesh->getNumVertices();

	numElements = tetMesh->getNumElements();
	undeformedPositions = (double*) malloc (sizeof(double) * 3 * numVertices);
	for(int i=0; i < numVertices; i++)
	{
		Vec3d * v = tetMesh->getVertex(i);
		for(int j=0; j<3; j++)
			undeformedPositions[3*i+j] = (*v)[j];
	}

	MInverse = (double**) malloc (sizeof(double*) * numElements);
	for(int el = 0; el < numElements; el++)
	{
		int vtxIndex[4];
		for(int vtx=0; vtx<4; vtx++)
			vtxIndex[vtx] = tetMesh->getVertexIndex(el, vtx);
		/*
		Form matrix: 
		M = [ v0   v1   v2   v3 ]
			[  1    1    1    1 ]
		*/
		double M[16]; // row-major
		for(int vtx=0; vtx<4; vtx++)
			for(int dim=0; dim<3; dim++)
				M[4 * dim + vtx] = undeformedPositions[3 * vtxIndex[vtx] + dim];
		M[12] = M[13] = M[14] = M[15] = 1.0;

		MInverse[el] = (double*) malloc (sizeof(double) * 16);
		inverse4x4(M, MInverse[el]);
	}
}

void muscle::inverse4x4(double * A, double * AInv)
{
  // converted to C from Mathematica output   
  AInv[0] = -A[11] * A[14] * A[5] + A[10] * A[15] * A[5] + A[11] * A[13] * A[6] - A[10] * A[13] * A[7] - A[15] * A[6] * A[9] + A[14] * A[7] * A[9];
  AInv[1] = A[1] * A[11] * A[14] - A[1] * A[10] * A[15] - A[11] * A[13] * A[2] + A[10] * A[13] * A[3] + A[15] * A[2] * A[9] - A[14] * A[3] * A[9];
  AInv[2] = -A[15] * A[2] * A[5] + A[14] * A[3] * A[5] + A[1] * A[15] * A[6] - A[13] * A[3] * A[6] - A[1] * A[14] * A[7] + A[13] * A[2] * A[7];
  AInv[3] = A[11] * A[2] * A[5] - A[10] * A[3] * A[5] - A[1] * A[11] * A[6] + A[1] * A[10] * A[7] + A[3] * A[6] * A[9] - A[2] * A[7] * A[9];
  AInv[4] = A[11] * A[14] * A[4] - A[10] * A[15] * A[4] - A[11] * A[12] * A[6] + A[10] * A[12] * A[7] + A[15] * A[6] * A[8] - A[14] * A[7] * A[8];
  AInv[5] = -A[0] * A[11] * A[14] + A[0] * A[10] * A[15] + A[11] * A[12] * A[2] - A[10] * A[12] * A[3] - A[15] * A[2] * A[8] + A[14] * A[3] * A[8];
  AInv[6] = A[15] * A[2] * A[4] - A[14] * A[3] * A[4] - A[0] * A[15] * A[6] + A[12] * A[3] * A[6] + A[0] * A[14] * A[7] - A[12] * A[2] * A[7];
  AInv[7] = -A[11] * A[2] * A[4] + A[10] * A[3] * A[4] + A[0] * A[11] * A[6] - A[0] * A[10] * A[7] - A[3] * A[6] * A[8] + A[2] * A[7] * A[8];
  AInv[8] = -A[11] * A[13] * A[4] + A[11] * A[12] * A[5] - A[15] * A[5] * A[8] + A[13] * A[7] * A[8] + A[15] * A[4] * A[9] - A[12] * A[7] * A[9];
  AInv[9] = -A[1] * A[11] * A[12] + A[0] * A[11] * A[13] + A[1] * A[15] * A[8] - A[13] * A[3] * A[8] - A[0] * A[15] * A[9] + A[12] * A[3] * A[9];
  AInv[10] = -A[1] * A[15] * A[4] + A[13] * A[3] * A[4] + A[0] * A[15] * A[5] - A[12] * A[3] * A[5] + A[1] * A[12] * A[7] - A[0] * A[13] * A[7];
  AInv[11] = A[1] * A[11] * A[4] - A[0] * A[11] * A[5] + A[3] * A[5] * A[8] - A[1] * A[7] * A[8] - A[3] * A[4] * A[9] + A[0] * A[7] * A[9]; 
  AInv[12] = A[10] * A[13] * A[4] - A[10] * A[12] * A[5] + A[14] * A[5] * A[8] - A[13] * A[6] * A[8] - A[14] * A[4] * A[9] + A[12] * A[6] * A[9];
  AInv[13] = A[1] * A[10] * A[12] - A[0] * A[10] * A[13] - A[1] * A[14] * A[8] + A[13] * A[2] * A[8] + A[0] * A[14] * A[9] - A[12] * A[2] * A[9]; 
  AInv[14] = A[1] * A[14] * A[4] - A[13] * A[2] * A[4] - A[0] * A[14] * A[5] + A[12] * A[2] * A[5] - A[1] * A[12] * A[6] + A[0] * A[13] * A[6];
  AInv[15] = -A[1] * A[10] * A[4] + A[0] * A[10] * A[5] - A[2] * A[5] * A[8] + A[1] * A[6] * A[8] + A[2] * A[4] *A[9] - A[0] * A[6] * A[9];

  double invDet = 1.0 / (A[0] * AInv[0] + A[1] * AInv[4] + A[2] * AInv[8] + A[3] * AInv[12]);

  for(int i=0; i<16; i++)
    AInv[i] *= invDet;
}

void muscle::force(double * u, vector<double> length){
	double FT[9];
	for(int k = 0; k < 9; k++)
		Force[k]=0;
	for(int i=0;i<muscle_segments.size();i++){
		muscle_segments[i].force(FT, length[i]);
		int vtxIndex[4];
		int el = muscle_segments[i].centerIndex;
		for (int vtx=0; vtx<4; vtx++)
			vtxIndex[vtx] = tetMesh->getVertexIndex(el, vtx);
		double P[16];
		for(int i=0; i<3; i++)
			for(int j=0; j<4; j++)
				P[4 * i + j] = undeformedPositions[3 * vtxIndex[j] + i] + u[3 * vtxIndex[j] + i];
		// row 4
		for(int j=0; j<4; j++)
			P[12 + j] = 1;

		// F = P * Inverse(M)
		double F[9]; // upper-left 3x3 block
		for(int i=0; i<3; i++) 
			for(int j=0; j<3; j++) {
				F[3 * i + j] = 0;
				for(int k=0; k<4; k++)
					F[3 * i + j] += P[4 * i + k] * MInverse[el][4 * k + j];
			}

		double R[9]; // rotation (row-major)
		double S[9]; // symmetric (row-major)
		double det = PolarDecomposition::Compute(F, R, S, 1E-6);
		if (det < 0)
		{
			// flip R so that it becomes orthogonal
			for(int i=0; i<9; i++)
				R[i] *= -1.0;
		}

		double RT[9] = {R[0], R[3], R[6], R[1], R[4], R[7], R[2], R[3], R[8]};
		double FTT[9];
		multiply(R, FT, FTT);
		multiply(FTT, RT, FT);

		for(int k = 0; k < 9; k++)
			Force[k] += FT[k];
	}
}

double muscle::get_distance(int el){
	Vec3d center = tetMesh->getElementCenter(el);
	double minD = 1111111111;
	for(int i = 0; i < muscle_segments.size(); i++){
		int ci = muscle_segments[i].centerIndex;
		if(el != ci){
			Vec3d p = tetMesh->getElementCenter(ci);
			double len = len2(p - center);
			if(len < minD)
				minD = len;
		}
	}
	return minD;
}

muscle::~muscle(){
}


//////////////////////////////////////////////////////////////////////////////////////////////////


muscle_segment::muscle_segment(int ci, double len, Vec3d dir, double stiff){
	centerIndex = ci;
	init_length = len;
	direction = dir;
	cur_length = len;
	stiffness = stiff;
}


void muscle_segment::force(double * F, double length){
	cur_length = length;
	double f = stiffness * (cur_length - init_length);

	double U[9] = {direction[0], 0, 0, direction[1], 0, 0, direction[2], 0, 0 };
	double UT[9] = {direction[0], direction[1], direction[2], 0, 0, 0, 0, 0, 0};
	double FF[9] = {f, 0, 0, 0, 0, 0, 0, 0, 0};
	double FT[9];
	multiply(U, FF, FT);
	multiply(FT, UT, F);
}

muscle_segment::~muscle_segment(){
}