
/*-------------------------------------------
 * AUTHOR: Alexandre YAHI
 * AFFILIATION: Columbia University Medical Center
 * 2015-2018
 -------------------------------------------*/
 
/*
 * recnw.h for program RecNW
 *
 */

#include <iostream>
#include <string>
#include <algorithm>
#include <stdlib.h>
#include <stdio.h>
 #include <utility>
#include <limits>

#define DEBUG 0

using namespace std;

// max functions
float max_2(float, float);
float max_3(float, float, float);

// initialization functions
void  F_init_aff( float **, int, int, float, float, bool, bool);
void  F_init_lin( float **, int, int, float, bool, bool);
void  P_init( float **, int, int, float, float, bool);
void  Q_init( float **, int, int, float, float, bool);

// alignment functions
string recnw_affine(string, string, float, float, float, float, bool, bool, bool, bool, int, int);
string recnw_reg(string, string, float, float, float, bool, bool, bool, bool, int, int);
