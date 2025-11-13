#pragma once

#include"ode_solver.h"

matrix ff0T(matrix, matrix = NAN, matrix = NAN);
matrix ff0R(matrix, matrix = NAN, matrix = NAN);
matrix df0(double, matrix, matrix = NAN, matrix = NAN);

//lab_1 funkcja testowa
matrix ff_test(matrix, matrix = NAN, matrix = NAN);

//lab_2 zbiorniki
//

matrix ff_tanks(matrix x, matrix ud1 = NAN, matrix ud2 = NAN);

matrix ff2T(matrix x, matrix ud1 = NAN, matrix ud2 = NAN);
matrix df2(double t, matrix Y, matrix ud1 = NAN, matrix ud2 = NAN);

matrix ff2R(matrix x, matrix ud1 = NAN, matrix ud2 = NAN);

// ------------------
// lab3

matrix g1(matrix x, matrix = NAN);
matrix g2(matrix x, matrix = NAN);
matrix g3(matrix x, matrix a);