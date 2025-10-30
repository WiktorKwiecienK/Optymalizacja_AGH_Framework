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
}



double* expansion(matrix(*ff)(matrix, matrix, matrix),
                  double x0, double d, double alpha, int Nmax,
                  matrix ud1, matrix ud2)
{
    try
    {
        double* p = new double[2];
        int i = 0;

        solution s0(matrix(1, 1, x0));
        solution s1(matrix(1, 1, x0 + d));

        double f0 = m2d(s0.fit_fun(ff, ud1, ud2));
        double f1 = m2d(s1.fit_fun(ff, ud1, ud2));

        if (f1 == f0)
        {
            p[0] = std::min(x0, x0 + d);
            p[1] = std::max(x0, x0 + d);
            return p;
        }

        if (f1 > f0)
        {
            d = -d;
            s1.x(0, 0) = x0 + d;
            f1 = m2d(s1.fit_fun(ff, ud1, ud2));

            if (f1 > f0)
            {
                p[0] = std::min(x0 - fabs(d), x0 + fabs(d));
                p[1] = std::max(x0 - fabs(d), x0 + fabs(d));
                return p;
            }
        }

        while (true)
        {
            if (solution::f_calls > Nmax)
                throw string("expansion(...): przekroczono limit wywołań funkcji celu");

            double x_prev = x0 + pow(alpha, i - 1) * d;
            double x_curr = x0 + pow(alpha, i) * d;

            solution s_prev(matrix(1, 1, x_prev));
            solution s_curr(matrix(1, 1, x_curr));

            double f_prev = m2d(s_prev.fit_fun(ff, ud1, ud2));
            double f_curr = m2d(s_curr.fit_fun(ff, ud1, ud2));

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

                if (p[0] > p[1])
                    std::swap(p[0], p[1]);

                return p;
            }

            i++;
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
        vector<int> F = {1, 1};
        while (F.back() <= (b - a) / epsilon)
            F.push_back(F[F.size() - 1] + F[F.size() - 2]);
        int k = F.size() - 1;

        double a_k = a;
        double b_k = b;

        double c_k = b_k - (double)F[k - 1] / F[k] * (b_k - a_k);
        double d_k = a_k + b_k - c_k;

        solution sc(matrix(1, 1, c_k));
        solution sd(matrix(1, 1, d_k));

        double fc = m2d(sc.fit_fun(ff, ud1, ud2));
        double fd = m2d(sd.fit_fun(ff, ud1, ud2));

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
                sc.x(0, 0) = c_k;
                fc = m2d(sc.fit_fun(ff, ud1, ud2));
            }
            else
            {
                a_k = c_k;
                c_k = d_k;
                fc = fd;
                d_k = a_k + b_k - c_k;
                sd.x(0, 0) = d_k;
                fd = m2d(sd.fit_fun(ff, ud1, ud2));
            }
        }

        double x_star = (c_k + d_k) / 2.0;
        Xopt.x = matrix(1, 1, x_star);
        Xopt.fit_fun(ff, ud1, ud2);
        Xopt.flag = 0;

        return Xopt;
    }
    catch (string ex_info)
    {
        throw ("solution fib(...):\n" + ex_info);
    }
}

solution lag(matrix(*ff)(matrix, matrix, matrix),
              double a, double b, double epsilon, double gamma,
              int Nmax, matrix ud1, matrix ud2)
{
    solution Xopt;

    double a_k = a;
    double b_k = b;
    double c_k = (a + b) / 2.0;
    double d_k = 0.0, d_prev = 0.0;

    int iter = 0;

    while ((b_k - a_k) > epsilon)
    {
        if (solution::f_calls > Nmax)
        {
            Xopt.x = matrix(1, 1, d_k);
            Xopt.fit_fun(ff, ud1, ud2);
            Xopt.flag = 0;
            return Xopt;
        }

        solution Sa(matrix(1, 1, a_k));
        solution Sb(matrix(1, 1, b_k));
        solution Sc(matrix(1, 1, c_k));

        double fa = m2d(Sa.fit_fun(ff, ud1, ud2));
        double fb = m2d(Sb.fit_fun(ff, ud1, ud2));
        double fc = m2d(Sc.fit_fun(ff, ud1, ud2));

        double l = fa * (pow(b_k, 2) - pow(c_k, 2)) +
                   fb * (pow(c_k, 2) - pow(a_k, 2)) +
                   fc * (pow(a_k, 2) - pow(b_k, 2));

        double m = fa * (b_k - c_k) +
                   fb * (c_k - a_k) +
                   fc * (a_k - b_k);

        if (fabs(m) < 1e-12)
            d_k = (a_k + b_k) / 2.0;
        else
            d_k = l / (2.0 * m);

        // zabezpieczenie na wypadek degeneracji
        if (d_k <= a_k || d_k >= b_k)
            d_k = (a_k + b_k) / 2.0;

        if (fabs(d_k - d_prev) < gamma)
            break;

        d_prev = d_k;

        solution Sd(matrix(1, 1, d_k));
        double fd = m2d(Sd.fit_fun(ff, ud1, ud2));

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
        else
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

        // aktualizacja środka dla stabilności
        c_k = (a_k + b_k) / 2.0;
        iter++;
    }

    Xopt.x = matrix(1, 1, d_k);
    Xopt.fit_fun(ff, ud1, ud2);
    Xopt.flag = 0;
    return Xopt;
}


solution fib_zmodyfikowany(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, matrix ud1, matrix ud2)
{
    try
    {
        solution Xopt;

        // --- Przygotowanie zapisu do pliku ---
        std::ofstream logFile("fib_iterations.csv");
        logFile << "Iteracja;a_k;b_k;Dlugosc_przedzialu\n";

        // 1. Znajdź najmniejsze k takie, że F(k) > (b - a) / epsilon
        std::vector<int> F = {1, 1};
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

        // Zapis stanu początkowego
        logFile << 0 << ";" << a_k << ";" << b_k << ";" << fabs(b_k - a_k) << "\n";

        // --- Iteracje Fibonacciego ---
        for (int i = 0; i <= k - 3; ++i)
        {
            if (solution::f_calls > 100000)
                throw std::string("fib(...): przekroczono limit wywołań funkcji celu");

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

            // 🔹 Zapis bieżącego stanu (do wykresu)
            logFile << (i + 1) << ";" << a_k << ";" << b_k << ";" << fabs(b_k - a_k) << "\n";
        }

        logFile.close(); // Zamknij plik

        // --- Wyznacz wynik końcowy ---
        double x_star = (c_k + d_k) / 2.0;
        matrix mx(1, 1, x_star);
        Xopt.x = mx;
        Xopt.y = ff(mx, ud1, ud2);
        Xopt.flag = 0;

        return Xopt;
    }
    catch (std::string ex_info)
    {
        throw ("solution fib(...):\n" + ex_info);
    }
}
#include <fstream>
#include <cmath>
#include "solution.h"
#include "matrix.h"
solution lag_zmodyfikowany(matrix(*ff)(matrix, matrix, matrix),
                            double a, double b, double epsilon, double gamma,
                            int Nmax, matrix ud1, matrix ud2)
{
    solution Xopt;

    double a_k = a;
    double b_k = b;
    double c_k = (a + b) / 2.0;
    double d_prev = 0.0;
    int i = 0;

    // 🔹 Plik CSV do logowania iteracji
    std::ofstream logFile("lag_iterations.csv");
    logFile << "Iteracja;a_k;b_k;c_k;d_k;Dlugosc_przedzialu\n";

    // --- zapis stanu początkowego ---
    logFile << i << ";" << a_k << ";" << b_k << ";" << c_k << ";" << 0.0 << ";" << fabs(b_k - a_k) << "\n";

    while ((b_k - a_k) > epsilon)
    {
        if (solution::f_calls > Nmax)
        {
            Xopt.x = matrix(1, 1, 0.0);
            Xopt.y = matrix(1, 1, 1e10);
            Xopt.flag = -1;
            logFile.close();
            return Xopt;
        }

        solution Sa(matrix(1, 1, a_k));
        solution Sb(matrix(1, 1, b_k));
        solution Sc(matrix(1, 1, c_k));

        double fa = m2d(Sa.fit_fun(ff, ud1, ud2));
        double fb = m2d(Sb.fit_fun(ff, ud1, ud2));
        double fc = m2d(Sc.fit_fun(ff, ud1, ud2));

        double l = fa * (pow(b_k, 2) - pow(c_k, 2)) +
                   fb * (pow(c_k, 2) - pow(a_k, 2)) +
                   fc * (pow(a_k, 2) - pow(b_k, 2));

        double m = fa * (b_k - c_k) +
                   fb * (c_k - a_k) +
                   fc * (a_k - b_k);

        if (m <= 0)
        {
            Xopt.x = matrix(1, 1, 0.0);
            Xopt.y = matrix(1, 1, 1e10);
            Xopt.flag = -2;
            logFile.close();
            return Xopt;
        }

        double d_k = 0.5 * l / m;

        if (fabs(d_k - d_prev) < gamma)
            break;

        d_prev = d_k;
        solution Sd(matrix(1, 1, d_k));
        double fd = m2d(Sd.fit_fun(ff, ud1, ud2));

        // 🔹 Logowanie iteracji
        logFile << i + 1 << ";" << a_k << ";" << b_k << ";" << c_k << ";" << d_k << ";" << fabs(b_k - a_k) << "\n";

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
            Xopt.x = matrix(1, 1, 0.0);
            Xopt.y = matrix(1, 1, 1e10);
            Xopt.flag = -3;
            logFile.close();
            return Xopt;
        }

        i++;
    }

    logFile.close();

    solution Sopt(matrix(1, 1, d_prev));
    Sopt.fit_fun(ff, ud1, ud2);
    Sopt.flag = 0;
    return Sopt;
}



solution HJ(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
    try
    {
        solution Xopt;
        int n = get_len(x0);  // liczba zmiennych

        solution XB(x0);  // punkt bazowy
        XB.fit_fun(ff, ud1, ud2);

        solution XB_old(x0);  // poprzedni punkt bazowy

        while (true)
        {
            solution X = HJ_trial(ff, XB, s, ud1, ud2);

            if (solution::f_calls > Nmax)
            {
                Xopt = XB;
                Xopt.flag = 0;
                return Xopt;
            }

            if (X.y(0, 0) < XB.y(0, 0))
            {
                while (true)
                {
                    XB_old = XB;
                    XB = X;

                    matrix x_new = 2.0 * XB.x - XB_old.x;
                    solution X_temp(x_new);
                    X_temp.fit_fun(ff, ud1, ud2);

                    X = HJ_trial(ff, X_temp, s, ud1, ud2);

                    if (solution::f_calls > Nmax)
                    {
                        Xopt = XB;
                        Xopt.flag = 0;
                        return Xopt;
                    }

                    if (X.y(0, 0) >= XB.y(0, 0))
                        break;
                }

                XB = X;
            }
            else
            {
                s = alpha * s;
            }

            if (s < epsilon)
                break;
        }

        Xopt = XB;
        Xopt.flag = 0;
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
        int n = get_len(XB.x);
        matrix e_j(n, 1);

        solution X = XB;

        for (int j = 0; j < n; ++j)
        {
            for (int i = 0; i < n; ++i)
                e_j(i) = (i == j) ? 1.0 : 0.0;

            solution X_plus(X.x + s * e_j);
            X_plus.fit_fun(ff, ud1, ud2);

            if (X_plus.y(0, 0) < X.y(0, 0))
            {
                X = X_plus;
            }
            else
            {
                // Próba kroku w kierunku -e_j
                solution X_minus(X.x - s * e_j);
                X_minus.fit_fun(ff, ud1, ud2);

                if (X_minus.y(0, 0) < X.y(0, 0))
                {
                    X = X_minus;
                }
            }
        }

        return X;
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
        int n = get_len(x0);

        int i = 0;
        matrix d(n, n);
        for (int j = 0; j < n; ++j)
            for (int k = 0; k < n; ++k)
                d(j, k) = (j == k) ? 1.0 : 0.0;

        matrix lambda(n, 1);
        for (int j = 0; j < n; ++j)
            lambda(j) = 0.0;

        matrix p(n, 1);
        for (int j = 0; j < n; ++j)
            p(j) = 0.0;

        matrix s = s0;
        solution XB(x0);
        XB.fit_fun(ff, ud1, ud2);

        while (true)
        {
            for (int j = 0; j < n; ++j)
            {
                // Wektor kierunku d_j
                matrix d_j = get_col(d, j);


                solution X_plus(XB.x + s(j) * d_j);
                X_plus.fit_fun(ff, ud1, ud2);

                if (X_plus.y(0, 0) < XB.y(0, 0))
                {
                    XB = X_plus;
                    lambda(j) = lambda(j) + s(j);
                    s(j) = s(j) * alpha;
                }
                else
                {

                    solution X_minus(XB.x - s(j) * d_j);
                    X_minus.fit_fun(ff, ud1, ud2);

                    if (X_minus.y(0, 0) < XB.y(0, 0))
                    {
                        XB = X_minus;
                        lambda(j) = lambda(j) - s(j);
                        s(j) = s(j) * alpha;
                    }
                    else
                    {

                        s(j) = s(j) * beta;
                        p(j) = p(j) + 1;
                    }
                }
            }

            i = i + 1;

            if (solution::f_calls > Nmax)
            {
                Xopt = XB;
                Xopt.flag = 0;
                return Xopt;
            }


            bool change_basis = true;
            for (int j = 0; j < n; ++j)
            {
                if (lambda(j) == 0.0 && p(j) == 0.0)
                {
                    change_basis = false;
                    break;
                }
            }

            if (change_basis)
            {

                matrix v(n, n);

                for (int j = 0; j < n; ++j)
                {
                    matrix v_j(n, 1);
                    for (int k = 0; k <= j; ++k)
                    {
                        matrix d_k = get_col(d, k);
                        v_j = v_j + lambda(k) * d_k;
                    }
                    v.set_col(v_j, j);
                }

                for (int j = 0; j < n; ++j)
                {
                    matrix v_j = get_col(v, j);
                    double norm_v = norm(v_j);

                    if (norm_v > 1e-12)
                    {
                        matrix d_j_new = v_j * (1.0 / norm_v);
                        d.set_col(d_j_new, j);
                    }
                }

                // Reset lambda i p
                for (int j = 0; j < n; ++j)
                {
                    lambda(j) = 0.0;
                    p(j) = 0.0;
                }


                s = s0;
            }

//warunek stopu
            double max_s = 0.0;
            for (int j = 0; j < n; ++j)
            {
                if (fabs(s(j)) > max_s)
                    max_s = fabs(s(j));
            }

            if (max_s < epsilon)
                break;
        }

        Xopt = XB;
        Xopt.flag = 0;
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
