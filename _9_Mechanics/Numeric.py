#include "Numeric.h"


// �]���藝�̌����D
double Numeric::Law_of_cosines(double a, double b, double c) {
	double cos_th = (b*b + c*c - a*a) / (2*b*c);
	return cos_th;
}

// 2�悷��D
double Numeric::Square(double x) {
	return x * x;
}

// x��3��
double Numeric::Cube(double x) {
	return x * x * x;
}

// �������ʂ��l�̌ܓ����Cint�^�ɂ��ďo��
int Numeric::Roundoff(double x) {
	int i = (int)(x + 0.5);
	if (x < -0.5)
		return i - 1;
	return i;
}