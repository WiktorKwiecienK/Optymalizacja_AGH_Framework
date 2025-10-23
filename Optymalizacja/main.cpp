/*********************************************
Kod stanowi uzupe?nienie materia??w do ?wicze?
w ramach przedmiotu metody optymalizacji.
Kod udost?pniony na licencji CC BY-SA 3.0
Autor: dr in?. ?ukasz Sztangret
Katedra Informatyki Stosowanej i Modelowania
Akademia G?rniczo-Hutnicza
Data ostatniej modyfikacji: 30.09.2025
*********************************************/

#include"opt_alg.h"

void lab0();
void lab1();

	void zadanie_5a();
	void zadanie_5b();

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

		//lab1();

	//	double test_Da = 50.0;
	//	cout<< "\n\nTest DA = " << test_Da <<"\n\n"; <sprawdza poprawno?? modelu (dla Da = 50 T_max ~= 62.5)
	//	test_DA(test_Da);


		// ???
	//zadanie_5a();

		//archiwalna
		//lab1_kacper();

		zadanie_5b();

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
//ZADANIE PROBLEM RZECZYWISTY


void zadanie_5a()
{
	srand(time(NULL)); // losowo?? punktów startowych

	const int N = 100;
	const double eps = 1e-5;
	const int Nmax = 1000;
	const double d = 10.0;
	const double lag_lambda = 1.0;
	const int lag_iter = 100;

	double alphas[3] = {1.2, 1.5, 2.0};

	ofstream fout("wyniki_zadanie_5a.csv");
	fout << "alpha;run;x0;FIB_x;FIB_y;FIB_calls;LAG_x;LAG_y;LAG_calls;FIBexp_x;FIBexp_y;FIBexp_calls;LAGexp_x;LAGexp_y;LAGexp_calls" << endl;

	for (int a = 0; a < 3; ++a)
	{
		double alpha = alphas[a];

		for (int i = 0; i < N; ++i)
		{
			double x0 = 1.0 + (rand() / (double)RAND_MAX) * 99.0; // losowy punkt startowy z [1,100]

			// Bez ekspansji
			solution Sfib = fib(ff_test, 0, 100, eps, NAN, NAN);
			solution Slag = lag(ff_test, 0, 100, eps, lag_lambda, lag_iter, NAN, NAN);
			cout<<"Punkt startowy losowy = "<<x0<<endl;
			// Z ekspansj?
			double* interval = expansion(ff_test, x0, d, alpha, Nmax, NAN, NAN);

			cout<<"Ekspansja :"<<interval[0]<<" - "<<interval[1]<<" Dla alpha = "<<alpha<<endl;
			solution SeFib = fib(ff_test, interval[0], interval[1], eps, NAN, NAN);
			cout<<SeFib<<endl<<endl;

			//solution SeLag = lag(ff_test, interval[0], interval[1], eps, lag_lambda, lag_iter, NAN, NAN);

			//cout<<"LAG => "<<SeLag<<endl;
			// Zapis do pliku
			// fout << fixed << setprecision(6)
			// << alpha << ";" << i + 1 << ";" << x0 << ";"
			// << Sfib.x << ";" << Sfib.y << ";" << Sfib.f_calls << ";"
			// << Slag.x << ";" << Slag.y << ";" << Slag.f_calls << ";"
			// << SeFib.x << ";" << SeFib.y << ";" << SeFib.f_calls << ";"
			// << SeLag.x << ";" << SeLag.y << ";" << SeLag.f_calls << endl;

			delete[] interval;
		}
	}

	fout.close();
	cout << "Zapisano wyniki do pliku: wyniki_zadanie_5a.csv" << endl;
}



void zadanie_5b()
{
	const double eps = 1e-5;
	const double gamma = 1.0;
	const int Nmax = 100;

	// === OPTIMALIZACJA ===
	solution Sfib = fib(ff_tanks, 1.0, 100.0, eps, NAN, NAN);
	solution Slag = lag(ff_tanks, 1.0, 100.0, eps, gamma, Nmax, NAN, NAN);

	// === ZAPIS DO TABELI 3 ===
	ofstream fout("tabela3.csv");
	fout << "Metoda;x;y;f_calls\n";
	fout << fixed << setprecision(6);
	fout << "Fibonacci;" << Sfib.x(0,0) << ";" << Sfib.y(0,0) << ";" << Sfib.f_calls << "\n";
	fout << "Lagrange;" << Slag.x(0,0) << ";" << Slag.y(0,0) << ";" << Slag.f_calls << "\n";
	fout.close();
	cout << "Zapisano wyniki optymalizacji do pliku: tabela3.csv" << endl;

	// === SYMULACJA DLA OPTYMALNEGO D_A ===
	double D_A_cm2 = Sfib.x(0,0); // lub Slag.x(0,0) je?li chcesz porówna?

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
	fsim << "t;T_B;V_A;V_B\n";

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

		fsim << t << ";" << T_B << ";" << V_A << ";" << V_B << "\n";
	}

	fsim.close();
	cout << "Zapisano przebieg symulacji do pliku: symulacja.csv" << endl;
}


//ZADANIE Z FUNKCJA TESTOWA



void lab1()
{
	try
	{

		cout << "=====>>>>> FUNKCJA TESTOWA <<<<<<=====" << endl;

		//Lag() i Fib() Dla przedzialu 0 - 100
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
		double d = 7.0;
		double alpha = 1.5;
		int Nmax = 1000;

		double* interval = expansion(ff_test, x0, d, alpha, Nmax, NAN, NAN);
		cout << "Przedzia? po ekspansji: [" << interval[0] << ", " << interval[1] << "]" << endl;


		//lag() i fib() dla przedzia?u znalezionego przez metode ekspansji
		solution SeFib = fib(ff_test, interval[0], interval[1], 1e-5, NAN, NAN);
		cout << "Minimum (Fibonacci po ekspansji):" << endl;
		cout << SeFib << endl;

		solution Selag = lag(ff_test, interval[0], interval[1], 1e-5, 1.0, 100, NAN, NAN);
		cout << "Minimum (Lagrange'a po ekspansji):" << endl;
		cout << Selag << endl;


		cout<<"<><><><><> ZADANIE Z FUNKCJA TESTOWA <><><><><><><>\n\n";
		zadanie_5a();

		//delete[] interval;

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
