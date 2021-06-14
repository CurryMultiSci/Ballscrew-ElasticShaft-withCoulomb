#include "Numeric.h"


// 余弦定理の公式．
double Numeric::Law_of_cosines(double a, double b, double c) {
	double cos_th = (b*b + c*c - a*a) / (2*b*c);
	return cos_th;
}

// 2乗する．
double Numeric::Square(double x) {
	return x * x;
}

// xの3乗
double Numeric::Cube(double x) {
	return x * x * x;
}

// 小数第一位を四捨五入し，int型にして出力
int Numeric::Roundoff(double x) {
	int i = (int)(x + 0.5);
	if (x < -0.5)
		return i - 1;
	return i;
}