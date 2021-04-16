#pragma once

#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>

using namespace std;


double fax(double x, double r, double G, double ms);
double fay(double y, double r, double G, double ms);

double fk13(double x, double r, double G, double ms);
double fk14(double y, double r, double G, double ms);

double fk21(double k11, double dt, double k13);
double fk22(double k12, double dt, double k14);
double fk23(double k11, double x0, double dt, double r, double G, double ms);
double fk24(double k12, double y0, double dt, double r, double G, double ms);

double fk31(double k11, double k23, double dt);
double fk32(double k12, double k24, double dt);
double fk33(double x0, double dt, double r, double G, double ms, double k21);
double fk34(double y0, double dt, double r, double G, double ms, double k22);

double fk41(double k11, double k33, double dt);
double fk42(double k12, double k34, double dt);
double fk43(double x0, double dt, double r, double G, double ms, double k31);
double fk44(double y0, double dt, double r, double G, double ms, double k32);

void calculateX(vector<double>& tabX, double x0, double vx0, double dt, double r, double G, double ms);
void calculateY(vector<double>& tabY, double y0, double vy0, double dt, double r, double G, double ms);

double calculateA(double a, double x0, double y0, double x, double y);

double MAX(double a, double b);
#pragma once
