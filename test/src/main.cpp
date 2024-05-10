#include <iostream>
#include <iomanip>

#include "math/matrix.h"

using namespace std;
using namespace FI;
using namespace MATH;

void printRM( const Matrix &m ) {
	int s = 13;
	const float *o = m.data();
	cout << right;
	cout << setw(s)<<o[0] << setw(s)<<o[1] << setw(s)<<o[2] << setw(s)<<o[3] << endl;
	cout << setw(s)<<o[4] << setw(s)<<o[5] << setw(s)<<o[6] << setw(s)<<o[7] << endl;
	cout << setw(s)<<o[8] << setw(s)<<o[9] << setw(s)<<o[10]<< setw(s)<<o[11]<< endl;
	cout << setw(s)<<o[12]<< setw(s)<<o[13]<< setw(s)<<o[14]<< setw(s)<<o[15]<< endl;
}

void printCM( const Matrix &m ) {
	int s=12;
	const float *o = m.data();
	cout << right << setprecision(2) << setiosflags(ios::fixed) << setfill(' ');
	cout << setw(s)<<o[0] << setw(s)<<o[4] << setw(s)<<o[8] << setw(s)<<o[12] << endl;
	cout << setw(s)<<o[1] << setw(s)<<o[5] << setw(s)<<o[9] << setw(s)<<o[13] << endl;
	cout << setw(s)<<o[2] << setw(s)<<o[6] << setw(s)<<o[10] << setw(s)<<o[14] << endl;
	cout << setw(s)<<o[3] << setw(s)<<o[7] << setw(s)<<o[11] << setw(s)<<o[15] << endl;
}

int main(int argc, char *argv[]) 
{

	cout << "Matrix Test." << endl;
	
	float m[16] = { 1.0f,  1.0f,  1.0f,  1.0f,
			-1.0f,  2.0f,  2.0f,  2.0f,
			-1.0f, -2.0f,  3.0f,  3.0f,
			-1.0f, -2.0f, -3.0f,  4.0f };

	Matrix mat1;
	mat1 = m;
	cout << "Original Matrix: " << endl;

	printCM(mat1);
	mat1.transpose();
	cout << "transposed matrix:" << endl;
	printCM(mat1);

	mat1 = m;
	Matrix mat2 = mat1;
	cout << "matrix 2 created from matrix 1: " << endl;
	printCM(mat2);

	cout << "determinant = " << mat1.determinant() << endl;

	cout << "inverse: " << endl;
	mat2.invert();
	printCM(mat2);
	mat2 = m;
	cout << "input matrix multiplied by its inverse is:" << endl;
	printCM(mat1*mat2.inverse());
}

