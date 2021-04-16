#define _USE_MATH_DEFINES


#include "Function.h"

using namespace std;


double fax(double x, double r, double G, double ms) {
	double ax;
	ax = -G * ms * x / pow(r, 3);
	return ax;
}

double fay(double y, double r, double G, double ms) {
	double ay;
	ay = -G * ms * y / pow(r, 3);
	return ay;
}


double fk13(double x, double r, double G, double ms) {
	return fax(x, r, G, ms);
}

double fk14(double y, double r, double G, double ms) {
	return fay(y, r, G, ms);
}


double fk21(double k11, double dt, double k13) {
	double vx;
	vx = k11 + (dt / 2) * k13;
	return vx;
}

double fk22(double k12, double dt, double k14) {
	double vy;
	vy = k12 + (dt / 2) * k14;
	return vy;
}

double fk23(double k11, double x0, double dt, double r, double G, double ms) {
	double ax, x;
	x = x0 + (dt / 2) * k11;
	ax = fax(x, r, G, ms);
	return ax;
}

double fk24(double k12, double y0, double dt, double r, double G, double ms) {
	double ay, y;
	y = y0 + (dt / 2) * k12;
	ay = fay(y, r, G, ms);
	return ay;
}


double fk31(double k11, double k23, double dt) {
	double vx;
	vx = k11 + (dt / 2) * k23;
	return vx;
}

double fk32(double k12, double k24, double dt) {
	double vy;
	vy = k12 + (dt / 2) * k24;
	return vy;
}

double fk33(double x0, double dt, double r, double G, double ms, double k21) {
	double ax, x;
	x = x0 + (dt / 2) * k21;
	ax = fax(x, r, G, ms);
	return ax;
}

double fk34(double y0, double dt, double r, double G, double ms, double k22) {
	double ay, y;
	y = y0 + (dt / 2) * k22;
	ay = fay(y, r, G, ms);
	return ay;
}


double fk41(double k11, double k33, double dt) {
	double vx;
	vx = k11 + dt * k33;
	return vx;
}

double fk42(double k12, double k34, double dt) {
	double vy;
	vy = k12 + dt * k34;
	return vy;
}

double fk43(double x0, double dt, double r, double G, double ms, double k31) {
	double ax, x;
	x = x0 + dt * k31;
	ax = fax(x, r, G, ms);
	return ax;
}

double fk44(double y0, double dt, double r, double G, double ms, double k32) {
	double ay, y;
	y = y0 + dt * k32;
	ay = fay(y, r, G, ms);
	return ay;
}




void calculateX(vector<double>& tabX, double x0, double vx0, double dt, double r, double G, double ms) {

	double k11, k13, k21, k23, k31, k33, k41, k43;
	double x, vx;


	k11 = vx0;
	k13 = fk13(x0, r, G, ms);

	k21 = fk21(k11, dt, k13);
	k23 = fk23(k11, x0, dt, r, G, ms);

	k31 = fk31(k11, k23, dt);
	k33 = fk33(x0, dt, r, G, ms, k21);

	k41 = fk41(k11, k33, dt);
	k43 = fk43(x0, dt, r, G, ms, k31);

	x = x0 + (dt / 6) * (k11 + 2 * k21 + 2 * k31 + k41);
	tabX.push_back(x);

	vx = vx0 + (dt / 6) * (k13 + 2 * k23 + 2 * k33 + k43);
	tabX.push_back(vx);

	return;
}

void calculateY(vector<double>& tabY, double y0, double vy0, double dt, double r, double G, double ms) {

	double k12, k14, k22, k24, k32, k34, k42, k44;
	double y, vy;

	k12 = vy0;
	k14 = fk14(y0, r, G, ms);

	k22 = fk22(k12, dt, k14);
	k24 = fk24(k12, y0, dt, r, G, ms);

	k32 = fk32(k12, k24, dt);
	k34 = fk34(y0, dt, r, G, ms, k22);

	k42 = fk42(k12, k34, dt);
	k44 = fk44(y0, dt, r, G, ms, k32);

	y = y0 + (dt / 6) * (k12 + 2 * k22 + 2 * k32 + k42);
	tabY.push_back(y);

	vy = vy0 + (dt / 6) * (k14 + 2 * k24 + 2 * k34 + k44);
	tabY.push_back(vy);

	return;
}


double calculateA(double a, double x0, double y0, double x, double y) {
	double da;
	da = abs(atan(y0 / x0) - atan(y / x));

	if (da <= 0) {
		return -1.0;
	}
	return a + da;
}


double MAX(double a, double b) {
	if (a > b) {
		return a;
	}
	else {
		return b;
	}
}