// Kamil Paluszewski 180194
#define  _CRT_SECURE_NO_WARNINGS
#include <math.h>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>

using namespace std;

#define N_PROBEK 512

typedef struct Mat {
	double** vals;
	int size;
} Mat;

void InitMatrix(Mat* mat, int size)
{
	mat->size = size;
	mat->vals = (double**)malloc(size * sizeof(double*));

	for (int i = 0; i < size; i++)
	{
		mat->vals[i] = (double*)malloc(size * sizeof(double));
	}

	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			mat->vals[i][j] = 0;
		}
	}
}

void ForwardSubstitution(Mat* mat, double* y, double* b, int size)
{
	for (int i = 0; i < size; i++)
	{
		double sum = 0.0;
		for (int j = 0; j <= i - 1; j++)
		{
			sum += mat->vals[i][j] * y[j];
		}

		y[i] = (b[i] - sum) / mat->vals[i][i];
	}
}

void BackwardSubstitution(Mat* mat, double* x, double* y, int size)
{
	for (int i = size - 1; i >= 0; i--)
	{
		double sum = 0.0;
		for (int j = i + 1; j < size; j++)
		{
			sum += mat->vals[i][j] * x[j];
		}

		x[i] = (y[i] - sum) / mat->vals[i][i];
	}
}

void CopyMatrix(Mat* from, Mat* to) // deep copy macierzy
{
	to->size = from->size;

	for (int i = 0; i < from->size; i++)
	{
		for (int j = 0; j < from->size; j++)
		{
			to->vals[i][j] = from->vals[i][j];
		}
	}
}

double* CrossProduct(Mat* mat, double* vector)
{
	int size = mat->size;
	double* result = (double*)malloc(size * sizeof(double));

	for (int i = 0; i < size; i++)
	{
		double value = 0.0;

		for (int j = 0; j < size; j++)
		{
			value += mat->vals[i][j] * vector[j];
		}
		result[i] = value;
	}

	return result;
}

void Pivoting(Mat* L, Mat* U, Mat* P, int j)
{
	double pivot = abs(U->vals[j][j]);
	int index = j; // indeks pivota

	for (int i = j + 1; i < U->size; i++)
	{
		if (abs(U->vals[i][j]) > pivot)
		{
			pivot = abs(U->vals[i][j]);
			index = i;
		}
	}

	if (index != j)
	{
		for (int i = 0; i < U->size; i++)
		{
			if (i < j)
			{
				double temp = L->vals[j][i];
				L->vals[j][i] = L->vals[index][i];
				L->vals[index][i] = temp;
			}
			else
			{
				double temp = U->vals[j][i];
				U->vals[j][i] = U->vals[index][i];
				U->vals[index][i] = temp;
			}

			double temp = P->vals[j][i];
			P->vals[j][i] = P->vals[index][i];
			P->vals[index][i] = temp;

		}
	}

}

void LU(Mat* mat, double* b, double* x)
{
	int size = mat->size;

	Mat L, U, P; // P = macierz permutacji
	InitMatrix(&L, size);

	InitMatrix(&U, size);
	CopyMatrix(mat, &U);

	InitMatrix(&P, size);

	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			if (i == j)
			{
				L.vals[i][j] = 1;
				P.vals[i][j] = 1;
			}
		}
	}

	double* y = (double*)malloc(size * sizeof(double));

	for (int i = 0; i < size - 1; i++)
	{
		Pivoting(&L, &U, &P, i);

		for (int j = i + 1; j < size; j++)
		{
			L.vals[j][i] = U.vals[j][i] / U.vals[i][i];

			for (int k = i; k < size; k++)
			{
				U.vals[j][k] = U.vals[j][k] - L.vals[j][i] * U.vals[i][k];
			}
		}
	}

	b = CrossProduct(&P, b);

	ForwardSubstitution(&L, y, b, size);

	BackwardSubstitution(&U, x, y, size);

	free(y);
	for (int j = 0; j < size; j++)
	{
		free(L.vals[j]);
		free(U.vals[j]);
		free(P.vals[j]);
	}
	free(L.vals);
	free(U.vals);
	free(P.vals);
}

int Wczytaj(const char* pobrane, double** probki, char separator)
{
	ifstream inFile;
	string nazwa = "dane/";
	nazwa += pobrane;

	inFile.open(nazwa, inFile.in);
	if (!inFile.is_open())
	{
		cout << "Nie bylo mozliwe poprawne otworzenie pliku";
		return 1;
	}
	string s;

	for (int i = 0; i < N_PROBEK; i++)
	{

		getline(inFile, s, separator);
		probki[i][0] = stod(s);
		getline(inFile, s);
		probki[i][1] = stod(s);
	}
	return 0;
}

void Lagrange(double** wezly, int n_wezlow, const char* pobrane, const char* rozmieszczenie, double** probki) {

	string plik_wyjsciowy = "output/Lagrange/";

	plik_wyjsciowy += to_string(n_wezlow);
	plik_wyjsciowy += "_";
	plik_wyjsciowy += rozmieszczenie;
	plik_wyjsciowy += "_";
	plik_wyjsciowy += pobrane;
	plik_wyjsciowy += ".csv";

	cout << "Lagrange rozmieszczenie " << rozmieszczenie << " " << n_wezlow << " wezlow" << endl;

	FILE* out;
	out = fopen(plik_wyjsciowy.c_str(), "w");

	double RMSD = 0.0;

	for (int i = 0; i < N_PROBEK; i++)
	{
		double wysokosc = 0.0;

		for (int j = 0; j < n_wezlow; j++)
		{
			double a = 1.0;

			for (int k = 0; k < j; k++)
			{
				a = a * (probki[i][0] - wezly[k][0]) / (wezly[j][0] - wezly[k][0]);
			}

			for (int k = j + 1; k < n_wezlow; k++)
			{
				a = a * (probki[i][0] - wezly[k][0]) / (wezly[j][0] - wezly[k][0]);
			}

			wysokosc = wysokosc + a * wezly[j][1];
		}

		double delta = wysokosc - probki[i][1];

		RMSD += pow(delta, 2);

		char result[200];
		// w pliku csv: odleg³oœæ (x), wysokosc rzeczywista(zmierzona), wysokosc ESTYMOWANA w trzeciej kolumnie
		//(w celu przyspieszenia tworzenia wykresu)

		sprintf_s(result, 200, "%f,%f,%f \n", probki[i][0], probki[i][1], wysokosc);
		fputs(result, out);
	}

	RMSD = RMSD / N_PROBEK;
	RMSD = sqrt(RMSD);

	cout << "BLAD SREDIOKWADRATOWY: " << RMSD << endl << endl;

	fclose(out);
}

void SplinesMethod(double** wezly, int n_wezlow, const char* pobrane, const char* rozmieszczenie, double** probki)
{
	string plik_wyjsciowy = "output/splajny/";

	plik_wyjsciowy += to_string(n_wezlow);
	plik_wyjsciowy += "_";
	plik_wyjsciowy += rozmieszczenie;
	plik_wyjsciowy += "_";
	plik_wyjsciowy += pobrane;
	plik_wyjsciowy += ".csv";

	cout << "Splajny rozmieszczenie " << rozmieszczenie << " " << n_wezlow << " wezlow" << endl;

	FILE* out;
	out = fopen(plik_wyjsciowy.c_str(), "w");

	double RMSD = 0.0;

	int size = (n_wezlow - 1) * 4;
	Mat A;
	InitMatrix(&A, size);

	double* b = (double*)malloc(size * sizeof(double));
	double* x = (double*)malloc(size * sizeof(double));

	for (int j = 0; j < size; j++)
	{
		x[j] = 1;
		b[j] = 0;
	}

	//tworzenie macierzy

	//(1)  S0(x0) = f(x0) => a0=f(x0)
	A.vals[0][0] = 1;
	b[0] = wezly[0][1];

	//(2) S0(x1) = f(x1) => a0 + b0*h + c0*h^2 + d0*h^3 = f(x1)
	double h = wezly[1][0] - wezly[0][0];
	A.vals[1][0] = 1;
	A.vals[1][1] = h;
	A.vals[1][2] = h * h;
	A.vals[1][3] = h * h * h;
	b[1] = wezly[1][1];


	for (int j = 1; j < n_wezlow - 1; j++) {

		// rozmieszczenie nie musi byæ równomierne, dlatego rozdzielam warunek x1-x0=x2-x1=h na h1 i h2
		// (poprawia to jakoœæ interpolacji z losowo wybranymi wêz³ami)
		double h1 = wezly[j][0] - wezly[j - 1][0];
		double h2 = wezly[j + 1][0] - wezly[j][0];

		//(3) Sj(xj) = f(xj) => aj = f(xj)
		A.vals[4 * j - 2][4 * j] = 1;
		b[4 * j - 2] = wezly[j][1];

		//(4) Sj(xj+1) = f(xj+1) => aj + bj*h + cj*h^2 + dj*h^3 = f(xj+1)
		A.vals[4 * j - 1][4 * j] = 1;
		A.vals[4 * j - 1][4 * j + 1] = h2;
		A.vals[4 * j - 1][4 * j + 2] = h2 * h2;
		A.vals[4 * j - 1][4 * j + 3] = h2 * h2 * h2;
		b[4 * j - 1] = wezly[j + 1][1];

		//(5) Sj-1'(xj) = Sj'(xj) => bj-1 + 2*(cj-1)*h + 3*(dj-1)*h^2 -bj = 0
		A.vals[4 * j][4 * (j - 1) + 1] = 1;
		A.vals[4 * j][4 * (j - 1) + 2] = 2 * h1;
		A.vals[4 * j][4 * (j - 1) + 3] = 3 * h1 * h1;
		A.vals[4 * j][4 * j + 1] = -1;
		b[4 * j] = 0;

		//(6) Sj-1''(xj) = Sj''(xj) => 2*(cj-1) + 6*(dj-1)*h - 2*cj = 0
		A.vals[4 * j + 1][4 * (j - 1) + 2] = 2;
		A.vals[4 * j + 1][4 * (j - 1) + 3] = 6 * h1;
		A.vals[4 * j + 1][4 * j + 2] = -2;
		b[4 * j + 1] = 0;
	}

	//(7) S0''(x0) = 0 => c0=0
	A.vals[4 * (n_wezlow - 1) - 2][2] = 1;
	b[4 * (n_wezlow - 1) - 2] = 0;

	h = wezly[n_wezlow - 1][0] - wezly[n_wezlow - 2][0];

	//(8) Sn-1''(xn) = 0  => 2*cn-1 + 6*dn-1 = 0
	A.vals[4 * (n_wezlow - 1) - 1][4 * (n_wezlow - 2) + 2] = 2;
	A.vals[4 * (n_wezlow - 1) - 1][4 * (n_wezlow - 2) + 3] = 6 * h;
	b[4 * (n_wezlow - 1) - 1] = 0;

	LU(&A, b, x);

	for (int i = 0; i < N_PROBEK; i++)
	{
		double wysokosc = 0.0;
		for (int j = 0; j < n_wezlow; j++)
		{
			if (probki[i][0] > wezly[j][0] && probki[i][0] < wezly[j + 1][0])
			{
				for (int k = 0; k < 4; k++)
				{
					double h = probki[i][0] - wezly[j][0];
					wysokosc += x[4 * j + k] * pow(h, k);

				}

				char result[200];
				/* w pliku csv: odleg³oœæ (x), wysokosc rzeczywista(zmierzona), wysokosc ESTYMOWANA w trzeciej kolumnie
					(w celu przyspieszenia tworzenia wykresu)*/
				sprintf_s(result, 200, "%f,%f,%f \n", probki[i][0], probki[i][1], wysokosc);
				fputs(result, out);
			}
			else if (probki[i][0] == wezly[j][0])
			{
				// nie interpoluj wartosci w samych wezlach (wpisz ich rzeczywista zmierzona wartosc)
				wysokosc = wezly[j][1];
				char result[200];
				// w pliku csv: odleg³oœæ (x), wysokosc rzeczywista(zmierzona), wysokosc ESTYMOWANA w trzeciej kolumnie
				//(aczkolwiek w wezlach wysokoscia po estymacji jest wysokosc rzeczywista w tym wezle)
				sprintf_s(result, 200, "%f,%f,%f \n", probki[i][0], probki[i][1], wysokosc);
				fputs(result, out);
			}
		}
		double delta = wysokosc - probki[i][1];
		RMSD += pow(delta, 2);
	}

	RMSD = RMSD / N_PROBEK;
	RMSD = sqrt(RMSD);

	cout << "BLAD SREDIOKWADRATOWY: " << RMSD << endl << endl;

	fclose(out);
	free(b);
	free(x);
}


int main()
{
	const char* pobrane = "genoa_rapallo.txt";
	char separator = ' '; // jakim znakiem rozdzielone sa dane w pliku wejœciowym

	double** probki;

	cout << pobrane << endl;

	probki = (double**)malloc(sizeof(double*) * N_PROBEK);

	for (int i = 0; i < N_PROBEK; i++)
	{
		probki[i] = (double*)malloc(sizeof(double) * 2);
	}

	Wczytaj(pobrane, probki, separator);

	int n_wezlow[3] = { 9,17,33 }; // najlepiej: dzielnik liczby próbek + 1

	for (int i = 0; i < 3; i++)
	{
		/* ROZMIESZCZENIE RÓWNOMIERNE */

		int co_ile = ceil((double)N_PROBEK / (n_wezlow[i] - 1));


		double** wezly = (double**)malloc(sizeof(double*) * n_wezlow[i]);

		for (int j = 0; j < n_wezlow[i]; j++)
		{
			wezly[j] = (double*)malloc(sizeof(double) * 2);
		}

		int idx = 0;
		for (int j = 0; j < N_PROBEK; j++)
		{
			if (j % co_ile == 0)
			{
				wezly[idx][0] = probki[j][0];
				wezly[idx][1] = probki[j][1];
				idx++;
			}
		}

		//ostatnia probka zostanie ostatnim wezlem
		wezly[n_wezlow[i] - 1][0] = probki[N_PROBEK - 1][0];
		wezly[n_wezlow[i] - 1][1] = probki[N_PROBEK - 1][1];

		Lagrange(wezly, n_wezlow[i], pobrane, "rownomierne", probki);
		SplinesMethod(wezly, n_wezlow[i], pobrane, "rownomierne", probki);

		/* ROZMIESZCZENIE LOSOWE */

		srand((unsigned int)time(NULL));

		int* wylosowane = (int*)malloc(sizeof(int) * n_wezlow[i]);

		//Losujê n-2 punktów -> jako pozosta³e 2 punkty okreœlam z góry krañce zakresu probek

		for (int j = 0; j < n_wezlow[i] - 2; j++)
		{
			bool unique = false;
			int wylosowana_probka = 0;

			while (unique == false)
			{
				unique = true;
				wylosowana_probka = rand() % (N_PROBEK - 2) + 1;

				for (int k = 0; k < j; k++)
				{
					if (wylosowana_probka == wylosowane[k])
					{
						unique = false;
					}
				}
			}
			wylosowane[j] = wylosowana_probka;
		}

		wylosowane[n_wezlow[i] - 1] = N_PROBEK - 1;
		wylosowane[n_wezlow[i] - 2] = 0;

		for (int j = n_wezlow[i] - 1; j >= 0; j--)
		{
			int probki_idx = -1; // indeks w tablicy próbek
			int idx = -1; // indeks w tablicy wylosowanych

			for (int k = 0; k < n_wezlow[i]; k++)
			{
				if (wylosowane[k] > probki_idx)
				{
					probki_idx = wylosowane[k];
					idx = k;
				}
			}

			wezly[j][0] = probki[probki_idx][0];
			wezly[j][1] = probki[probki_idx][1];

			wylosowane[idx] = -1;

		}

		Lagrange(wezly, n_wezlow[i], pobrane, "losowe", probki);
		SplinesMethod(wezly, n_wezlow[i], pobrane, "losowe", probki);

		free(wylosowane);

		for (int j = 0; j < n_wezlow[i]; j++)
		{
			free(wezly[j]);
		}
		free(wezly);
	}

	for (int j = 0; j < N_PROBEK; j++)
	{
		free(probki[j]);
	}
	free(probki);

	return 0;
}
