/*********************************************
Kod stanowi uzupe?nienie materia??w do ?wicze?
w ramach przedmiotu metody optymalizacji.
Kod udost?pniony na licencji CC BY-SA 3.0
Autor: dr in?. ?ukasz Sztangret
Katedra Informatyki Stosowanej i Modelowania
Akademia G?rniczo-Hutnicza
Data ostatniej modyfikacji: 30.09.2025
*********************************************/

#include <iomanip>

#include"opt_alg.h"

void lab0();
void lab1();
void lab2();
void lab3();
void lab4();
void lab5();
void lab6();

//void lab1_test_function();

void test_DA(double da);

int main()
{
	try
	{
		lab2();
	}
	catch (string EX_INFO)
	{
		cerr << "ERROR:\n";
		cerr << EX_INFO << endl << endl;
	}
	return 0;
}

void lab0()
{
	//Funkcja testowa
	double epsilon = 1e-2;									// dok?adno??
	int Nmax = 10000;										// maksymalna liczba wywo?a? funkcji celu
	matrix lb(2, 1, -5), ub(2, 1, 5),						// dolne oraz g?rne ograniczenie
		a(2, 1);											// dok?adne rozwi?zanie optymalne
	solution opt;											// rozwi?zanie optymalne znalezione przez algorytm
	a(0) = -1;
	a(1) = 2;
	opt = MC(ff0T, 2, lb, ub, epsilon, Nmax, a);			// wywo?anie procedury optymalizacji
	cout << opt << endl << endl;							// wypisanie wyniku
	solution::clear_calls();								// wyzerowanie licznik?w

	//Wahadlo
	Nmax = 1000;											// dok?adno??
	epsilon = 1e-2;											// maksymalna liczba wywo?a? funkcji celu
	lb = 0, ub = 5;											// dolne oraz g?rne ograniczenie
	double teta_opt = 1;									// maksymalne wychylenie wahad?a
	opt = MC(ff0R, 1, lb, ub, epsilon, Nmax, teta_opt);		// wywo?anie procedury optymalizacji
	cout << opt << endl << endl;							// wypisanie wyniku
	solution::clear_calls();								// wyzerowanie licznik?w

	//Zapis symulacji do pliku csv
	matrix Y0 = matrix(2, 1),								// Y0 zawiera warunki pocz?tkowe
		MT = matrix(2, new double[2] { m2d(opt.x), 0.5 });	// MT zawiera moment si?y dzia?aj?cy na wahad?o oraz czas dzia?ania
	matrix* Y = solve_ode(df0, 0, 0.1, 10, Y0, NAN, MT);	// rozwi?zujemy r?wnanie r?niczkowe
	ofstream Sout("symulacja_lab0.csv");					// definiujemy strumie? do pliku .csv
	Sout << hcat(Y[0], Y[1]);								// zapisyjemy wyniki w pliku
	Sout.close();											// zamykamy strumie?
	Y[0].~matrix();											// usuwamy z pami?ci rozwi?zanie RR
	Y[1].~matrix();
}

void lab1()
{
	try
	{
		srand(time(NULL));

		int N = 100;
		double eps = 1e-5;
		int Nmax = 200;
		double d = 5.0;
		double lag_lambda = 5.0;
		int lag_iter = 100;

		double alphas[3] = {1.2, 1.5, 2.0};

		ofstream fout("tabela1.csv");
		fout << "alpha;run;x0;a;b;FIB_x;FIB_y;LAG_x;LAG_y\n";

		for (int a = 0; a < 3; ++a)
		{
			double alpha = alphas[a];

			for (int i = 0; i < N; ++i)
			{
				double x0 = 1.0 + (rand() / (double)RAND_MAX) * 99.0;

				double* interval = expansion(ff_test, x0, d, alpha, Nmax, NAN, NAN);
				double a_exp = interval[0];
				double b_exp = interval[1];

				solution SeFib = fib(ff_test, a_exp, b_exp, eps, NAN, NAN);
				solution SeLag = lag(ff_test, a_exp, b_exp, eps, lag_lambda, lag_iter, NAN, NAN);

				fout << fixed << setprecision(3)
					 << alpha << ";" << i + 1 << ";" << x0 << ";"
					 << a_exp << ";" << b_exp << ";"
					 << SeFib.x(0,0) << ";" << SeFib.y(0,0) << ";"
					 << SeLag.x(0,0) << ";" << SeLag.y(0,0) << "\n";

				delete[] interval;
			}
		}

		fout.close();
		cout << "Zapisano dane do pliku: tabela1.csv" << endl;

		cout << "=====>>>>> FUNKCJA TESTOWA <<<<<<=====" << endl;


		cout << "===== METODA FIBONACCIEGO =====" << endl;
		solution Sfib = fib(ff_test, 0, 100, 1e-5, NAN, NAN);
		cout << "Minimum (Fibonacci):" << endl;
		cout << Sfib << endl;

			cout << "\n===== METODA LAGRANGE'A =====" << endl;
			solution Slag;

			//	for (int i = 0; i < 101; ++i)
			Slag = lag(ff_test, 0, 100, 1e-5, 1.0, 100, NAN, NAN);

			cout << "Minimum (Lagrange):" << endl;
			cout << Slag << endl;



		cout << "\n===== METODA EKSPANSJI =====" << endl;
		double x0 = 50.0;
		 d = 7.0;
		double alpha = 1.5;
		Nmax = 1000;

		double* interval = expansion(ff_test, x0, d, alpha, Nmax, NAN, NAN);
		cout << "Przedzia? po ekspansji: [" << interval[0] << ", " << interval[1] << "]" << endl;


		//lag() i fib() dla przedzia?u znalezionego przez metode ekspansji
		solution SeFib = fib(ff_test, interval[0], interval[1], 1e-5, NAN, NAN);
		cout << "Minimum (Fibonacci po ekspansji):" << endl;
		cout << SeFib << endl;

		solution Selag = lag(ff_test, interval[0], interval[1], 1e-5, 1.0, 100, NAN, NAN);
		cout << "Minimum (Lagrange'a po ekspansji):" << endl;
		cout << Selag << endl;

		// ------------------ZBIORNIKI----------------------
		cout << "\n\n=====>>>>> ZADANIE ZE ZBIORNIKAMI <<<<<<=====" << endl;

		cout << "\n===== TEST MAKSYMALNEJ TEMPERATURY dla DA = 50cm^2 =====" << endl;
		test_DA(50.0);

		cout << "\n===== METODA EKSPANSJI =====" << endl;

		// USTALONY PUNKT STARTOWY
		x0 = 60, d = 10.0, alpha = 1.8;
		Nmax = 1000;

		interval = expansion(ff_tanks, x0, d, alpha, Nmax, NAN, NAN);
		cout << "Wyliczony przedzial metoda ekspansji: [" << interval[0] << ", " << interval[1] << "]" << endl;

		cout << "\n===== METODA FIBONACCIEGO =====" << endl;
		Sfib = fib(ff_tanks, interval[0], interval[1], 1e-3, NAN, NAN);
		//solution Sfib = fib(ff_tanks, 0, 20, 1e-3, NAN, NAN);
		cout << Sfib << endl;

		cout << "\n===== METODA LAGRANGE'A =====" << endl;
		Slag = lag(ff_tanks, interval[0], interval[1], 1e-3, 1.0, 100, NAN, NAN);
		cout << Slag << endl;

		delete[] interval;



		  eps = 1e-5;
		 double gamma = 1.0;
		 Nmax = 100;

		// === OPTIMALIZACJA ===
		 Sfib = fib(ff_tanks, 1.0, 100.0, eps, NAN, NAN);
		 Slag = lag(ff_tanks, 1.0, 100.0, eps, gamma, Nmax, NAN, NAN);

		// === SYMULACJA DLA OPTYMALNEGO D_A ===
		double D_A_cm2 = Sfib.x(0,0);

		const double P_A = 2.0, V_A0 = 5.0, T_A0 = 95.0;
		const double P_B = 1.0, V_B0 = 1.0, T_B0 = 20.0;
		const double D_B = 36.5665 / 10000.0;
		const double F_in = 0.01, T_in = 20.0;
		const double a = 0.98, b = 0.63, g = 9.81;
		const double dt = 1.0, t_end = 2000.0;

		double D_A = D_A_cm2 / 10000.0;
		double V_A = V_A0, T_A = T_A0;
		double V_B = V_B0, T_B = T_B0;

		ofstream fsim("symulacja.csv");
		fsim << "t;T_B_FIb;V_A_Fib;V_B_Fib\n";

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

			if (V_A < 0) V_A = 0;
			if (V_B < 0) V_B = 0;

			if ((int)t%10 == 0)
				fsim << t << ";" << T_B << ";" << V_A << ";" << V_B << "\n";
		}

		fsim.close();
		cout << "Zapisano przebieg symulacji do pliku: symulacja.csv" << endl;

	}
	catch (string ex_info)
	{
		cerr << "B??d: " << ex_info << endl;
	}
}

/*
 * Funkcja zapisuj?ca wynik symulacji dla okre?lonego przekroju
 */

void test_DA(double Da)
{

	//double D_A_cm2 = 12.5;
	double D_A_cm2 = Da;

	//Zbiornik A -> zmienny przekroj
	double D_A = D_A_cm2 / 10000.0; // cm² ? m²
	const double P_A = 2.0, V_A0 = 5.0, T_A0 = 95.0;
	//Zbiornik B
	const double P_B = 1.0, V_B0 = 1.0, T_B0 = 20.0;
	const double D_B = 36.5665 / 10000.0;
	const double F_in = 0.01, T_in = 20.0;

	const double a = 0.98, b = 0.63, g = 9.81;
	const double dt = 1.0, t_end = 2000.0;

	double V_A = V_A0, T_A = T_A0;
	double V_B = V_B0, T_B = T_B0;
	double max_T_B = T_B0;

	ofstream fout("symulacja_DA_50.csv");
	fout << "t;T_B;V_B" << endl;

	for (double t = 0; t <= t_end; t += dt)
	{
		double h_A = V_A / P_A;
		double h_B = V_B / P_B;

		double flow_A_to_B = a * b * D_A * sqrt(2 * g * h_A);
		double flow_B_out  = a * b * D_B * sqrt(2 * g * h_B);

		V_A -= flow_A_to_B * dt;
		V_B += (flow_A_to_B + F_in - flow_B_out) * dt;

		//Wersja kacprowa
		//double dT_dt = (flow_A_to_B * T_A + F_in * T_in - flow_B_out * T_B) / V_B;
		//T_B += dT_dt * dt;

		//Wersja wiktorowa
		double dT_dt = (flow_A_to_B * (T_A - T_B) + F_in * (T_in - T_B)) / V_B;
		T_B += dT_dt * dt;

		if (T_B > max_T_B) max_T_B = T_B;
		if (V_A < 0) V_A = 0;
		if (V_B < 0) V_B = 0;

		fout << t << ";" << T_B << ";" << V_B << endl;
	}

	fout.close();

	cout << "Zapisano symulacj? do pliku: symulacja_DA_50.csv" << endl;
	cout << "Maksymalna temperatura w zbiorniku B: " << max_T_B << " °C" << endl;

}


void lab2()
{
    try
    {
        srand(time(NULL));

        int n = 2;
        double epsilon = 1e-3;
        int Nmax = 10000;

        int num_tests = 100;

        cout << "\n===== TESTY DLA RÓ?NYCH D?UGO?CI KROKU =====" << endl;
        double step_sizes[3] = {0.5, 0.1, 0.01};

    	 ofstream fout2("lab2_tabela1.csv");
        fout2 << "Step_size;i;x1_start;x2_start;x1_HJ;x2_HJ;y_HS;f_calls;lokalne/globalne;x1_ROS;x2_ROS;y_ROS;f_calls;globalne/lokalne\n";

        for (int s_idx = 0; s_idx < 3; ++s_idx)
        {
            double step = step_sizes[s_idx];

            cout << "\nD?ugo?? kroku s = " << step << endl;

            for (int i = 0; i < num_tests; ++i)
            {
                double x1_start = -1.0 + (rand() / (double)RAND_MAX) * 2.0;
                double x2_start = -1.0 + (rand() / (double)RAND_MAX) * 2.0;
                matrix x0(2, 1);
                x0(0) = x1_start;
                x0(1) = x2_start;

                // HOOKE-JEEVES
                solution::clear_calls();
                solution result_HJ = HJ(ff2T, x0, step, 0.5, epsilon, Nmax, NAN, NAN);

                // ROSENBROCK
                solution::clear_calls();
                matrix s0(2, 1);
                s0(0) = step;
                s0(1) = step;
                solution result_Rosen = Rosen(ff2T, x0, s0, 2.0, 0.5, epsilon, Nmax, NAN, NAN);

			//	cout<<"RP = ("<<x1_start<<";"<<x2_start<<")\n";
            //	cout<<"HJ = ("<<result_HJ.x<<")\n";
            //	cout<<"R = ("<<result_Rosen.x<<")\n\n";
				// cout<<get_row(result_Rosen.x, 0)<<";"<<get_row(result_Rosen.x,1)<<endl<<endl;
				// fout2<<step<<";"<<i<<";"<<get_row(result_Rosen.x,0)<<get_row(result_Rosen.x,1)<<endl;


				double x1_HJ = m2d(get_row(result_HJ.x, 0));
            	double x2_HJ = m2d(get_row(result_HJ.x, 1));

            	double x1_R = m2d(get_row(result_Rosen.x, 0));
            	double x2_R = m2d(get_row(result_Rosen.x, 1));

            	bool isGlobalHJ = false;
            	bool isGlobalRosen = false;

            	if (fabs(x1_HJ) < 0.25 && fabs(x2_HJ) < 0.25)
            		isGlobalHJ = true;

            	if (fabs(x1_R) < 0.25 && fabs(x2_R) < 0.25)
            		isGlobalRosen = true;

            	string minHJ = isGlobalHJ ? "globalne" : "lokalne";
            	string minRosen = isGlobalRosen ? "globalne" : "lokalne";


				fout2<<step<<";"<<i<<";"
				<<x1_start<<";"<<x2_start<<";"
            	<<x1_HJ<<";"<<x2_HJ<<";"<< (double)m2d(result_HJ.y) <<";"
            	<<result_HJ.f_calls<<";"<<minHJ<<";"<<x1_R<<";"<<x2_R
            	<<";"<<(double)m2d(result_Rosen.y)<<";"<<result_Rosen.f_calls<< ";"<<minRosen<<endl;
            }
        }
	fout2.close();
    }
    catch (string ex_info)
    {
        cerr << "B??d: " << ex_info << endl;
    }
}

void lab3()
{

}

void lab4()
{

}

void lab5()
{

}

void lab6()
{

}
