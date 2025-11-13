#include"user_funs.h"

//l0
//------------------
matrix ff0T(matrix x, matrix ud1, matrix ud2)				// funkcja celu dla przypadku testowego
{
	matrix y;												// y zawiera wartość funkcji celu
	y = pow(x(0) - ud1(0), 2) + pow(x(1) - ud1(1), 2);		// ud1 zawiera współrzędne szukanego optimum
	return y;
}

matrix ff0R(matrix x, matrix ud1, matrix ud2)				// funkcja celu dla problemu rzeczywistego
{
	matrix y;												// y zawiera wartość funkcji celu
	matrix Y0 = matrix(2, 1),								// Y0 zawiera warunki początkowe
		MT = matrix(2, new double[2] { m2d(x), 0.5 });		// MT zawiera moment siły działający na wahadło oraz czas działania
	matrix* Y = solve_ode(df0, 0, 0.1, 10, Y0, ud1, MT);	// rozwiązujemy równanie różniczkowe
	int n = get_len(Y[0]);									// długość rozwiązania
	double teta_max = Y[1](0, 0);							// szukamy maksymalnego wychylenia wahadła
	for (int i = 1; i < n; ++i)
		if (teta_max < Y[1](i, 0))
			teta_max = Y[1](i, 0);
	y = abs(teta_max - m2d(ud1));							// wartość funkcji celu (ud1 to założone maksymalne wychylenie)
	Y[0].~matrix();											// usuwamy z pamięci rozwiązanie RR
	Y[1].~matrix();
	return y;
}
matrix df0(double t, matrix Y, matrix ud1, matrix ud2)
{
	matrix dY(2, 1);										// definiujemy wektor pochodnych szukanych funkcji
	double m = 1, l = 0.5, b = 0.5, g = 9.81;				// definiujemy parametry modelu
	double I = m * pow(l, 2);
	dY(0) = Y(1);																// pochodna z położenia to prędkość
	dY(1) = ((t <= ud2(1)) * ud2(0) - m * g * l * sin(Y(0)) - b * Y(1)) / I;	// pochodna z prędkości to przyspieszenie
	return dY;
}


//test function l2
	matrix ff_test(matrix x, matrix ud1, matrix ud2)
{
	double xv = x(0, 0);
	double val = -cos(0.1 * xv) * exp(-pow(0.1 * xv - 2 * M_PI, 2))
				 + 0.002 * pow(0.1 * xv, 2);
	return matrix(val);
}



matrix ff_tanks(matrix x, matrix ud1, matrix ud2)
{

	double D_A_cm2 = x(0, 0);
	if (D_A_cm2 < 1.0 || D_A_cm2 > 100.0)
		return matrix(1e10);

	double D_A = D_A_cm2 / 10000.0;

	const double P_A = 2.0, V_A0 = 5.0, T_A0 = 95.0;
	const double P_B = 1.0, V_B0 = 1.0, T_B0 = 20.0;
	const double D_B = 36.5665 / 10000.0;
	const double F_in = 0.01, T_in = 20.0;
	const double a = 0.98, b = 0.63, g = 9.81;
	const double dt = 1.0, t_end = 2000.0;

	double V_A = V_A0, T_A = T_A0;
	double V_B = V_B0, T_B = T_B0;
	double max_T_B = T_B0;

	for (double t = 0; t <= t_end; t += dt)
	{
		double h_A = V_A / P_A;
		double h_B = V_B / P_B;

 		double flow_A_to_B = a * b * D_A * sqrt(2 * g * h_A);
		double flow_B_out  = a * b * D_B * sqrt(2 * g * h_B);

		V_A -= flow_A_to_B * dt;
		V_B += (flow_A_to_B + F_in - flow_B_out) * dt;


		if (V_B < 1e-8) V_B = 1e-8;

		double dT_dt = (flow_A_to_B * (T_A - T_B) + F_in * (T_in - T_B)) / V_B;
		T_B += dT_dt * dt;

		if (T_B > max_T_B) max_T_B = T_B;
		if (V_A < 0) V_A = 0;
		if (V_B < 0) V_B = 0;
	}

	return matrix(abs(max_T_B - 50.0));
}


// Funkcja testowa z lab2 - f(x1, x2) = x1^2 + x2^2 - cos(2.5πx1) - cos(2.5πx2) + 2

matrix ff2T(matrix x, matrix ud1, matrix ud2)
{
	return matrix(
		pow(x(0), 2) + pow(x(1), 2)
		- cos(2.5 * M_PI * x(0))
		- cos(2.5 * M_PI * x(1))
		+ 2
	);
}

matrix df2(double t, matrix Y, matrix ud1, matrix ud2)
{
	double m_r = 1.0;      // masa ramienia
	double m_c = 5.0;      // masa ciężarka
	double l = 2.0;        // długość ramienia (m)
	double b = 0.25;       // współczynnik tarcia (Nms)

	double I = (1.0 / 3.0 * m_r + m_c) * pow(l, 2);

	double k1 = m2d(ud1);
	double k2 = m2d(ud2);

	double alpha_ref = M_PI;   // pozycja docelowa (rad)
	double omega_ref = 0.0;    // prędkość docelowa (rad/s)

	double alpha_error = alpha_ref - Y(0);
	double omega_error = omega_ref - Y(1);

	double M_t = k1 * alpha_error + k2 * omega_error;

	matrix dY(2, 1);
	dY(0) = Y(1);
	dY(1) = (M_t - b * Y(1)) / I;

	return dY;
}

matrix ff2R(matrix x, matrix ud1, matrix ud2)
{
	double k1 = x(0);
	double k2 = x(1);

	double t0 = 0.0;
	double t_end = 100.0;
	double dt = 0.1;

	double alpha0 = 0.0;
	double omega0 = 0.0;
	matrix Y0 = matrix(2, new double[2]{alpha0, omega0});

	matrix* res = solve_ode(df2, t0, dt, t_end, Y0, k1, k2);

	double Q = 0.0;
	for (int i = 0; i < get_len(res[0]); i++)
	{
		double alpha = res[1](i, 0);
		double omega = res[1](i, 1);

		double alpha_error = M_PI - alpha;
		double omega_error = 0.0 - omega;

		double M_t = k1 * alpha_error + k2 * omega_error;

		Q += (10 * pow(alpha_error, 2)
			+ pow(omega_error, 2)
			+ pow(M_t, 2)) * dt;
	}

	delete[] res;

	return matrix(Q);
}

// ------------------
// lab3

matrix g1(matrix x, matrix ud1)
{
	return 1 -x(0);
}

matrix g2(matrix x, matrix ud1)
{
	return  1 - x(1) ;
}

matrix g3(matrix x, matrix a)
{
	return sqrt(pow(x(0), 2) + pow(x(1), 2)) - m2d(a);
}

//f testowa




//