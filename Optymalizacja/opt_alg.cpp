#include"opt_alg.h"

solution MC(matrix(*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	// Zmienne wejściowe:
	// ff - wskaźnik do funkcji celu
	// N - liczba zmiennych funkcji celu
	// lb, ub - dolne i górne ograniczenie
	// epslion - zakłądana dokładność rozwiązania
	// Nmax - maksymalna liczba wywołań funkcji celu
	// ud1, ud2 - user data
	try
	{
		solution Xopt;
		while (true)
		{
			Xopt = rand_mat(N);									// losujemy macierz Nx1 stosując rozkład jednostajny na przedziale [0,1]
			for (int i = 0; i < N; ++i)
				Xopt.x(i) = (ub(i) - lb(i)) * Xopt.x(i) + lb(i);// przeskalowywujemy rozwiązanie do przedziału [lb, ub]
			Xopt.fit_fun(ff, ud1, ud2);							// obliczmy wartość funkcji celu
			if (Xopt.y < epsilon)								// sprawdzmy 1. kryterium stopu
			{
				Xopt.flag = 1;									// flaga = 1 ozancza znalezienie rozwiązanie z zadaną dokładnością
				break;
			}
			if (solution::f_calls > Nmax)						// sprawdzmy 2. kryterium stopu
			{
				Xopt.flag = 0;									// flaga = 0 ozancza przekroczenie maksymalne liczby wywołań funkcji celu
				break;
			}
		}
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution MC(...):\n" + ex_info);
	}
}double* expansion(matrix(*ff)(matrix, matrix, matrix),
                  double x0, double d, double alpha, int Nmax,
                  matrix ud1, matrix ud2)
{
    try
    {
        double* p = new double[2];
        int i = 0;

        matrix x0_mat(1, 1, x0);
        matrix x1_mat(1, 1, x0 + d);

        double f0 = m2d(ff(x0_mat, ud1, ud2));
        double f1 = m2d(ff(x1_mat, ud1, ud2));

        // --- 1. Jeśli funkcja jest płaska ---
        if (f1 == f0)
        {
            p[0] = std::min(x0, x0 + d);
            p[1] = std::max(x0, x0 + d);
            return p;
        }

        // --- 2. Jeśli funkcja rośnie, zmień kierunek kroku ---
        if (f1 > f0)
        {
            d = -d;
            x1_mat(0, 0) = x0 + d;
            f1 = m2d(ff(x1_mat, ud1, ud2));

            // Jeśli nadal rośnie – minimum w pobliżu x0
            if (f1 > f0)
            {
                p[0] = std::min(x0 - fabs(d), x0 + fabs(d));
                p[1] = std::max(x0 - fabs(d), x0 + fabs(d));
                return p;
            }
        }

        // --- 3. Ekspansja w kierunku malejącej funkcji ---
        while (true)
        {
            if (++i > Nmax)
                throw string("expansion(...): przekroczono limit wywołań funkcji celu");

            double x_prev = x0 + pow(alpha, i - 1) * d;
            double x_curr = x0 + pow(alpha, i) * d;

            matrix x_prev_mat(1, 1, x_prev);
            matrix x_curr_mat(1, 1, x_curr);

            double f_prev = m2d(ff(x_prev_mat, ud1, ud2));
            double f_curr = m2d(ff(x_curr_mat, ud1, ud2));

            // znaleziono koniec przedziału
            if (f_prev <= f_curr)
            {
                if (d > 0)
                {
                    p[0] = x_prev;
                    p[1] = x_curr;
                }
                else
                {
                    p[0] = x_curr;
                    p[1] = x_prev;
                }

                // upewnij się, że p[0] < p[1]
                if (p[0] > p[1])
                    std::swap(p[0], p[1]);

                return p;
            }
        }
    }
    catch (string ex_info)
    {
        throw ("double* expansion(...):\n" + ex_info);
    }
}


solution fib(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;

		// 1. Znajdź najmniejsze k takie, że F(k) > (b - a) / epsilon
		vector<int> F = {1, 1};
		while (F.back() <= (b - a) / epsilon)
			F.push_back(F[F.size() - 1] + F[F.size() - 2]);
		int k = F.size() - 1;

		double a_k = a;
		double b_k = b;

		double c_k = b_k - (double)F[k - 1] / F[k] * (b_k - a_k);
		double d_k = a_k + b_k - c_k;

		matrix mc(1, 1, c_k);
		matrix md(1, 1, d_k);
		double fc = m2d(ff(mc, ud1, ud2));
		double fd = m2d(ff(md, ud1, ud2));

		for (int i = 0; i <= k - 3; ++i)
		{
			if (solution::f_calls > 100000)
				throw string("fib(...): przekroczono limit wywołań funkcji celu");

			if (fc < fd)
			{
				b_k = d_k;
				d_k = c_k;
				fd = fc;
				c_k = b_k - (double)F[k - i - 2] / F[k - i - 1] * (b_k - a_k);
				mc(0, 0) = c_k;
				fc = m2d(ff(mc, ud1, ud2));
			}
			else
			{
				a_k = c_k;
				c_k = d_k;
				fc = fd;
				d_k = a_k + b_k - c_k;
				md(0, 0) = d_k;
				fd = m2d(ff(md, ud1, ud2));
			}
		}

		double x_star = (c_k + d_k) / 2.0;
		matrix mx(1, 1, x_star);
		Xopt.x = mx;
		Xopt.y = ff(mx, ud1, ud2);
		Xopt.flag = 0;

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution fib(...):\n" + ex_info);
	}
}

solution lag(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, double gamma, int Nmax, matrix ud1, matrix ud2)
{
    try
    {
        solution Xopt;

        double a_k = a;
        double b_k = b;
        double c_k = (a + b) / 2.0;
        double d_prev = 0.0;
        int i = 0;

        while ((b_k - a_k) > epsilon)
        {
            if (solution::f_calls > Nmax)
                throw string("lag(...): przekroczono limit wywołań funkcji celu");

            matrix ma(1, 1, a_k);
            matrix mb(1, 1, b_k);
            matrix mc(1, 1, c_k);

            double fa = m2d(ff(ma, ud1, ud2));
            double fb = m2d(ff(mb, ud1, ud2));
            double fc = m2d(ff(mc, ud1, ud2));

            double l = fa * (pow(b_k, 2) - pow(c_k, 2)) +
                       fb * (pow(c_k, 2) - pow(a_k, 2)) +
                       fc * (pow(a_k, 2) - pow(b_k, 2));

            double m = fa * (b_k - c_k) +
                       fb * (c_k - a_k) +
                       fc * (a_k - b_k);

            if (m <= 0)
                throw string("lag(...): dzielnik m <= 0, nie można kontynuować interpolacji");

            double d_k = 0.5 * l / m;

            if (fabs(d_k - d_prev) < gamma)
                break;

            d_prev = d_k;
            matrix md(1, 1, d_k);
            double fd = m2d(ff(md, ud1, ud2));

            if (a_k < d_k && d_k < c_k)
            {
                if (fd < fc)
                {
                    b_k = c_k;
                    c_k = d_k;
                }
                else
                {
                    a_k = d_k;
                }
            }
            else if (c_k < d_k && d_k < b_k)
            {
                if (fd < fc)
                {
                    a_k = c_k;
                    c_k = d_k;
                }
                else
                {
                    b_k = d_k;
                }
            }
            else
            {
                throw string("lag(...): punkt d_k poza przedziałem interpolacji");
            }

            i++;
        }

        matrix mx(1, 1, d_prev);
        Xopt.x = mx;
        Xopt.y = ff(mx, ud1, ud2);
        Xopt.flag = 0;

        return Xopt;
    }
    catch (string ex_info)
    {
        throw ("solution lag(...):\n" + ex_info);
    }
}



solution HJ(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution HJ(...):\n" + ex_info);
	}
}

solution HJ_trial(matrix(*ff)(matrix, matrix, matrix), solution XB, double s, matrix ud1, matrix ud2)
{
	try
	{
		//Tu wpisz kod funkcji

		return XB;
	}
	catch (string ex_info)
	{
		throw ("solution HJ_trial(...):\n" + ex_info);
	}
}

solution Rosen(matrix(*ff)(matrix, matrix, matrix), matrix x0, matrix s0, double alpha, double beta, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Rosen(...):\n" + ex_info);
	}
}

solution pen(matrix(*ff)(matrix, matrix, matrix), matrix x0, double c, double dc, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try {
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution pen(...):\n" + ex_info);
	}
}

solution sym_NM(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double beta, double gamma, double delta, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution sym_NM(...):\n" + ex_info);
	}
}

solution SD(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution SD(...):\n" + ex_info);
	}
}

solution CG(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution CG(...):\n" + ex_info);
	}
}

solution Newton(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix),
	matrix(*Hf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Newton(...):\n" + ex_info);
	}
}

solution golden(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution golden(...):\n" + ex_info);
	}
}

solution Powell(matrix(*ff)(matrix, matrix, matrix), matrix x0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Powell(...):\n" + ex_info);
	}
}

solution EA(matrix(*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, int mi, int lambda, matrix sigma0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution EA(...):\n" + ex_info);
	}
}
