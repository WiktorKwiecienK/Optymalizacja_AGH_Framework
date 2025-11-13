#include"opt_alg.h"

#include <fstream>

solution HJ_CSV(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double epsilon, int Nmax,
	matrix ud1, matrix ud2) {
	{
		try
		{
			vector<double> trajectoryX1;
			vector<double> trajectoryX2;

			cout << "\nHJ_CSV function\n";

			solution Xopt, xs(x0);
			matrix ud(x0);

			xs.fit_fun(ff, ud1, ud2);
			Xopt.flag = 0;

			do {
				solution xB = xs;

				xs = HJ_trial(ff, xB, s, ud1, ud2);
				xs.fit_fun(ff, ud1, ud2);

				// Zapis aktualnego punktu
				//	if (xs.x.size() >= 2) {
				trajectoryX1.push_back(xs.x(0, 0));
				trajectoryX2.push_back(xs.x(1, 0));
				//	}

				if (m2d(xs.y) < m2d(xB.y)) {
					do {
						solution _xB = xB;
						xB = xs;
						xs = (xB.x * 2.0) - _xB.x;
						xs = HJ_trial(ff, xs, s, ud1, ud2);

						// Zapis aktualnego punktu
						//if (xs.x.size() >= 2) {
						trajectoryX1.push_back(xs.x(0, 0));
						trajectoryX2.push_back(xs.x(1, 0));
						//}

						if (solution::f_calls > Nmax) {
							Xopt.flag = -1;
							break;
						}
					} while (m2d(xs.y) < m2d(xB.y));
					xs = xB;
				}
				else {
					s = alpha * s;
				}

				if (solution::f_calls > Nmax) {
					if (!Xopt.flag)
						Xopt.flag = -2;
					break;
				}
				ud.add_col(xs.x);
			} while (s >= epsilon);

			Xopt = xs;
			Xopt.ud = ud;

			// --- ZAPIS TRAJEKTORII DO PLIKU CSV ---
			ofstream trajFile("trajectory_HJ.csv");
			if (trajFile.is_open()) {
				trajFile << "x1;x2\n";  // nagłówek kolumn
				for (size_t i = 0; i < trajectoryX1.size(); ++i) {
					trajFile << trajectoryX1[i] << ";" << trajectoryX2[i] << "\n";
				}
				trajFile.close();
				cout << "Trajektoria zapisana do pliku trajectory_HJ.csv\n";
			} else {
				cerr << "Błąd: nie udało się otworzyć pliku trajectory.csv do zapisu.\n";
			}

			return Xopt;
		}
		catch (string ex_info) {
			throw ("solution HJ(...):\n");
		}
	}
}

#include <fstream>  // do zapisu pliku CSV

solution Rosen_CSV(matrix(*ff)(matrix, matrix, matrix), matrix x0, matrix s0, double alpha, double beta, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		const int DIM = 2;
		matrix d(DIM, DIM);
		matrix ud(x0);

		for (int w = 0; w < DIM; w++)
			for (int k = 0; k < DIM; k++)
				d(w, k) = (w == k) ? 1 : 0;

		matrix l(DIM, 1, 0.0);
		matrix p(DIM, 1, 0.0);
		matrix s(s0);

		solution xB(x0);
		xB.fit_fun(ff, ud1, ud2);
		solution Xopt(xB);
		Xopt.flag = 0;
		int max_s;

		// --- wektory do zapisu trajektorii ---
		vector<double> trajectoryX1;
		vector<double> trajectoryX2;
		vector<double> trajectoryF;

		// Zapis punktu początkowego
		trajectoryX1.push_back(xB.x(0, 0));
		trajectoryX2.push_back(xB.x(1, 0));
		trajectoryF.push_back(m2d(xB.y));

		do
		{
			for (int j = 0; j < DIM; j++)
			{
				solution _x(xB.x + s(j) * d[j]);
				if (_x.fit_fun(ff, ud1, ud2) < xB.y)
				{
					xB = _x;
					l(j) = l(j) + s(j);
					s(j) = s(j) * alpha;
				}
				else
				{
					s(j) = -s(j) * beta;
					p(j) = p(j) + 1;
				}

				// Zapis aktualnego punktu po każdej zmianie
				trajectoryX1.push_back(xB.x(0, 0));
				trajectoryX2.push_back(xB.x(1, 0));
				trajectoryF.push_back(m2d(xB.y));
			}

			Xopt = xB;

			bool zero = false;
			for (int j = 0; j < DIM; j++)
			{
				if (p(j) == 0 || abs(l(j)) < epsilon) {
					zero = true;
					break;
				}
			}

			if (!zero)
			{
				matrix _D(d);
				matrix _lQ(DIM, DIM);

				for (int i = 0; i < DIM; i++)
					for (int j = 0; j < DIM; j++)
						_lQ(i, j) = (i >= j) ? l(i) : 0.0;

				_lQ = _D * _lQ;

				matrix v(DIM, DIM);
				v.set_col(_lQ[0] / (norm(_lQ[0])), 0);

				for (int _j = 1; _j < DIM; _j++)
				{
					matrix sigma(DIM, 1);
					matrix t_lQ(trans(_lQ[_j]));

					for (int k = 0; k < _j; k++)
					{
						sigma.set_col(
							sigma[0] + (t_lQ * d[k]) * d[k],
							0);
					}
					matrix pk = _lQ[_j] - sigma[0];
					v.set_col(pk / norm(pk), _j);
				}

				d = v;

				l = matrix(DIM, 1, 0.0);
				p = matrix(DIM, 1, 0.0);
				s = s0;
			}
			if (solution::f_calls > Nmax)
			{
				Xopt.flag = -2;
				break;
			}
			ud.add_col(Xopt.x);
			max_s = 0;
			for (int j = 1; j < DIM; j++)
			{
				if (abs(s(max_s)) < abs(s(j)))
				{
					max_s = j;
				}
			}

		} while (abs(s(max_s)) >= epsilon);

		Xopt.ud = ud;

		// --- ZAPIS TRAJEKTORII DO PLIKU CSV ---
		ofstream trajFile("trajectory_rosen.csv");
		if (trajFile.is_open()) {

			trajFile << "x1;x2;f\n";

			for (size_t i = 0; i < trajectoryX1.size(); ++i) {
				trajFile << trajectoryX1[i] << ";" << trajectoryX2[i] << ";" << trajectoryF[i] << "\n";
			}

			trajFile.close();
			cout << "Trajektoria zapisana do pliku trajectory_rosen.csv\n";

		}
		else {
			cerr << "Błąd: nie udało się otworzyć pliku trajectory_ROSEN.csv do zapisu.\n";
		}
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Rosen_CSV(...):\n" + ex_info);
	}
}


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
		solution Xopt, xs(x0);
		matrix ud(x0);

		xs.fit_fun(ff, ud1, ud2);
		Xopt.flag = 0;
		do{
			solution xB = xs;

			xs = HJ_trial(ff, xB, s, ud1, ud2);
			xs.fit_fun(ff, ud1, ud2);
			if (m2d(xs.y) < m2d(xB.y)){
				do{
					solution _xB = xB;
					xB = xs;
					xs = (xB.x * 2.0) - _xB.x;
					xs = HJ_trial(ff, xs, s, ud1, ud2);
					if (solution::f_calls > Nmax){
						Xopt.flag = -1;
						break;
					}
				} while (m2d(xs.y) < m2d(xB.y));
				xs = xB;
			}else{
				s = alpha * s;
			}
			if (solution::f_calls > Nmax){
				if(!Xopt.flag)
					Xopt.flag = -2;
				break;
			}
			ud.add_col(xs.x);
		} while (s >= epsilon);

		Xopt = xs;
		Xopt.ud = ud;
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution HJ(...):\n" + ex_info);
	}
}

solution HJ_trial(matrix(*ff)(matrix, matrix, matrix), solution XB, double s, matrix ud1, matrix ud2)
{
	try{
		int Wym = 2;
        matrix dj(Wym, Wym);

		for(int k = 0; k < Wym; k++) {
			for (int l = 0; l < Wym; l++) {
				if (l == k) {
					dj(k, l) = 1;
				}else {
					dj(k, l) = 0;
				}
			}
		}

		for (int j = 0; j < Wym; j++){
			solution xj = XB.x + (s * dj[j]);
			xj.fit_fun(ff, ud1, ud2);
			if ( xj.y < XB.y ){
				XB = xj;
			}else{
				xj = XB.x - (s * dj[j]);
				xj.fit_fun(ff, ud1, ud2);
				if (xj.y < XB.y){
					XB = xj;
				}
			}
		}
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

		const int DIM = 2;
		matrix d(DIM, DIM);
		matrix ud(x0);

		for (int w = 0; w < DIM; w++)
			for (int k = 0; k < DIM; k++)
				d(w, k) = (w == k) ? 1 : 0;

		matrix l(DIM, 1, 0.0);
		matrix p(DIM, 1, 0.0);
		matrix s(s0);

		solution xB(x0);
		xB.fit_fun(ff, ud1, ud2);
		solution Xopt(xB);
		Xopt.flag = 0;
		int max_s;
		do
		{
			for (int j = 0; j < DIM; j++)
			{
				solution _x(xB.x + s(j) * d[j]);
				if (_x.fit_fun(ff, ud1, ud2) <  xB.y)
				{
					xB = _x;
					l(j) = l(j) + s(j);
					s(j) = s(j) * alpha;
				}
				else
				{
					s(j) = -s(j) * beta;
					p(j) = p(j) + 1;
				}
			}
			Xopt = xB;

			bool zero = false;
			for (int j = 0; j < DIM; j++)
			{
				if (p(j) == 0 || abs(l(j)) < epsilon ) {
					zero = true;
					break;
				}
			}

			if (!zero)
			{
				matrix _D(d);
				matrix _lQ(DIM, DIM);

				for (int i = 0; i < DIM; i++)
					for (int j = 0; j < DIM; j++)
						_lQ(i, j) = (i >= j) ? l(i) : 0.0;

				_lQ = _D * _lQ;

				matrix v(DIM,DIM);
				v.set_col(_lQ[0]/(norm(_lQ[0])), 0);

				for (int _j = 1; _j < DIM; _j++)
				{
					matrix sigma(DIM,1);
					matrix t_lQ(trans(_lQ[_j]));

					for (int k = 0; k < _j; k++)
					{
						sigma.set_col(
							sigma[0] + (t_lQ * d[k]) * d[k],
							0);
					}
					matrix pk = _lQ[_j] - sigma[0];
					v.set_col(pk/norm(pk), _j);
				}

				d = v;

				l = matrix(DIM, 1, 0.0);
				p = matrix(DIM, 1, 0.0);
				s = s0;
			}
			if (solution::f_calls > Nmax)
			{
				Xopt.flag = -2;
				break;
			}
			ud.add_col(Xopt.x);
			max_s = 0;
			for (int j = 1; j < DIM; j++)
			{
				if (abs(s(max_s)) < abs(s(j)))
				{
					max_s = j;
				}
			}

		} while (abs(s(max_s)) >= epsilon);
		Xopt.ud = ud;
		return Xopt;

	}
	catch (string ex_info)
	{
		throw ("solution Rosen(...):\n" + ex_info);
	}
}

//penalty

solution pen(matrix(*ff)(matrix, matrix, matrix), matrix x0, double c, double dc, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try {

		solution xi(x0), x_i;
		xi.flag = 0;
		matrix init_v_S(2, 1);
		init_v_S(0) = c;
		init_v_S(1) = ud2(0);
		double nm = 0;
		//double tmp1 = m2d(ff3T(xi.x)), tmp2;
		do
		{
			x_i = xi;

			//tmp2 = tmp1;

			xi = sym_NM(ff,xi.x,ud1(0),ud1(1),ud1(2),ud1(3), ud1(4), ud1(5), Nmax, init_v_S);

			//tmp1 = m2d(ff3T(xi.x));
//std::c                                                   out << "PEN:\nx:" << x_i.x(0) << " " << x_i.x(1) << " y: " << xi.y << "\nx:" << xi.x(0) << " " << xi.x(1) << " y: " << xi.y << "\n";

			init_v_S(0) = init_v_S(0) * dc;

			if (solution::f_calls > Nmax)
			{
				xi.flag = -2;
				break;
			}
			if (dc < 1.0)
			{
				double sum = 0.0;
				sum += 1 / m2d(g1(xi.x, NAN));
				sum += 1 / m2d(g2(xi.x, NAN));
				sum += 1 / m2d(g3(xi.x, init_v_S(1)));
				if (c * fabs(sum) < epsilon)
					break;
			}

			nm = norm(xi.x - x_i.x);

		} while (nm >= epsilon);

		return xi;
	}
	catch (string ex_info)
	{
		throw ("matrix pen(...):\n" + ex_info);
	}
}

solution sym_NM(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double beta, double gamma, double delta, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{

		solution Xopt;
		int g_min = 0;
		Xopt.flag = 0;
		// set function's dimension
		const int DIM = 2;
		// set initial matrixes
		matrix p(DIM, 1 + DIM); // DIM (points' dimension) x DIM+1 (points' amount) square point matrix
		matrix e = ident_mat(DIM); // DIM x DIM square identity matrix

		p.set_col(x0, 0); // p0 = x0
		for (int i = 1; i <= DIM; i++)
			p.set_col(p[0] + e[i-1] * s, i); //pi = p0 + ei*s

		solution::f_calls += 1 + DIM;
		matrix p_f(DIM+1, 1);
		for (int i = 0; i <= DIM; i++)
			p_f(i) = m2d(ff(p[i], ud1, NAN)); // returns matrix 1x1

		double max_norm;
		do {
			max_norm = 0.0;
			int p_max = 0, p_min = 0;
			for (int i = 1; i <= DIM; i++) {
				if (p_f(p_max) < p_f(i)) p_max = i;
				if (p_f(p_min) > p_f(i)) p_min = i;
			}
			if (p_max == p_min)
				p_max = (p_max+1)%(DIM+1);

			matrix p_s(DIM,1); // _p
			for (int i = 0; i <= DIM; i++)
			{
				if (i == p_max) continue;
				p_s.set_col(p_s[0] + p[i], 0); // _p = E(i!=max) pi
			}

			p_s.set_col(p_s[0] / DIM, 0);  // _p /= n
			matrix p_odb = p_s[0] + (p_s[0] - p[p_max]) * alpha; // p_odb = _p + a(_p - p_max)

			solution::f_calls++;
			double p_odb_f = m2d(ff(p_odb, ud1, NAN)); // returns matrix 1x1

			if (m2d(p_odb_f) < p_f(p_max))
			{
				matrix p_e = p_s + (p_odb[0] - p_s[0]) * gamma;

				solution::f_calls++;
				double p_e_f = m2d(ff(p_e, ud1, NAN));

				if (ff(p_e, ud1, NAN) < p_odb_f)
				{
					p.set_col(p_e[0], p_max);
					p_f(p_max) = p_e_f;
				}
				else
				{
					p.set_col(p_odb[0], p_max);
					p_f(p_max) = p_odb_f;
				}
			}
			else
			{
				if (p_f(p_min) <= p_odb_f && p_odb_f < p_f(p_max))
				{
					p.set_col(p_odb[0], p_max);
					p_f(p_max) = p_odb(0);
				}
				else
				{
					matrix p_z = p_s[0] + (p[p_max] - p_s[0])*beta;

					solution::f_calls++;
					double p_z_f = m2d(ff(p_z, ud1, NAN));

					if (p_z_f >= p_f(p_max))
					{
						for (int i = 0; i <= DIM; i++)
						{
							if (i == p_min) continue;
							p.set_col((p[i] + p[p_min]) * delta, i);
							solution::f_calls++;
							p_f(i) = m2d(ff(p[i], ud1, NAN));
						}
					}
					else
					{
						p.set_col(p_z[0], p_max);
						p_f(p_max) = p_z_f;
					}
				}
			}
			g_min = p_min;
			if (solution::f_calls > Nmax)
			{
				Xopt.flag = -2;
				break;
				//throw("Nie znaleziono przedzialu po Nmax probach (f_calls > Nmax)\n");
			}
			for (int i = 0; i <= DIM; i++)
			{
				if (i == p_min) continue;
				double i_norm = norm(p[p_min] - p[i]);
				if (i_norm > max_norm)
					max_norm = i_norm;
			}
			//std::cout << max_norm << std::endl;
		} while (max_norm >= epsilon);
		Xopt.x = p[g_min];
		Xopt.y = p_f(g_min);
//std::cout << "NM: " << Xopt << std::endl;
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
