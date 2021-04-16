
#define _USE_MATH_DEFINES

#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include "Function.h"

using namespace std;



int main() {

    double x, y, x0, y0, r, vx, vy, vx0, vy0, au, ms, G;

    ms = 1.989 * pow(10, 30);
    G = 6.6741 * pow(10, -11);
    au = 149597870700;

    double t, dt, a, da;
    double c, ex, ey, e, tol;
    int n;
    double k11, k12, k13, k14;
    double k21, k22, k23, k24;
    double k31, k32, k33, k34;
    double k41, k42, k43, k44;

    vector<double> tabX;
    vector<double> tabY;

    vector<double> tabX0_5;
    vector<double> tabY0_5;

    vector<double> tabX1_0;
    vector<double> tabY1_0;

    fstream Data1;

    Data1.open("Data1.txt", ios_base::in | ios::_Nocreate);

    if (Data1.is_open()) {
        remove("Data1.txt");
        Data1.close();
        Data1.open("Data1.txt", ios::app);
    }
    else {
        Data1.open("Data1.txt", ios::app);
    }


    fstream Data2;

    Data2.open("Data2.txt", ios_base::in | ios::_Nocreate);

    if (Data2.is_open()) {
        remove("Data2.txt");
        Data2.close();
        Data2.open("Data2.txt", ios::app);
    }
    else {
        Data2.open("Data2.txt", ios::app);
    }


    for (int i = 0; i < 2; i++) {

        //The initial conditions 

        x0 = 0;
        y0 = 0.586 * au;
        x = x0;
        y = y0;

        r = sqrt(x0 * x0 + y0 * y0);

        vx0 = 54600;
        vy0 = 0;
        vx = vx0;
        vy = vy0;

        a = 0;
        t = 0;
        dt = 700;

        c = 0.9;
        n = 4;

        //part 1
        if (i == 0) {
            tol = 1000;
            Data1 << t << "," << dt << "," << x0 << "," << y0 << "," << r << "," << a << "\n";
        }

        //part 2
        else if (i == 1) {
            tol = 1;
            Data2 << t << "," << dt << "," << x0 << "," << y0 << "," << r << "," << a << "\n";
        }

        else {
            cout << "Incorrect i" << endl;
            return 0;
        }


        while (a < 6 * M_PI) {


            calculateX(tabX, x0, vx0, dt, r, G, ms);
            calculateY(tabY, y0, vy0, dt, r, G, ms);

            calculateX(tabX0_5, x0, vx0, dt / 2.0, r, G, ms);
            calculateY(tabY0_5, y0, vy0, dt / 2.0, r, G, ms);

            r = sqrt(pow(tabX0_5[0], 2) + pow(tabY0_5[0], 2));

            calculateX(tabX1_0, tabX0_5[0], tabX0_5[1], dt / 2.0, r, G, ms);
            calculateY(tabY1_0, tabY0_5[0], tabY0_5[1], dt / 2.0, r, G, ms);


            ex = (tabX1_0[0] - tabX[0]) / (pow(2, n) - 1.0);
            ey = (tabY1_0[0] - tabY[0]) / (pow(2, n) - 1.0);

            e = MAX(abs(ex), abs(ey));

            if (e < tol) {

                a = calculateA(a, x0, y0, tabX1_0[0], tabY1_0[0]);

                if (a == -1) {
                    cout << "Incorrect a" << endl;
                    break;
                }


                x0 = tabX1_0[0];
                y0 = tabY1_0[0];

                vx0 = tabX1_0[1];
                vy0 = tabY1_0[1];

                r = sqrt(x0 * x0 + y0 * y0);


                //part 1
                if (i == 0) {
                    Data1 << t << "," << dt << "," << x0 << "," << y0 << "," << r << "," << a << "\n";
                }

                //part 2
                else if (i == 1) {
                    Data2 << t << "," << dt << "," << x0 << "," << y0 << "," << r << "," << a << "\n";
                }

                else {
                    cout << "Incorrect i" << endl;
                    return 0;
                }

                t = t + dt;
                dt = c * dt * pow(tol / e, 1.0 / (n + 1.0));

            }

            else {
                dt = c * dt * pow(tol / e, 1.0 / (n + 1.0));
            }

            tabX.clear();
            tabY.clear();
            tabX0_5.clear();
            tabY0_5.clear();
            tabX1_0.clear();
            tabY1_0.clear();

        }
    }

    return 0;
}