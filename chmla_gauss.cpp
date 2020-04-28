#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
using namespace std;


// Макроси для вектора і матриці типу double і прототипи функцій
typedef vector<double> rvector;
typedef vector<rvector> rmatrix; 
rvector back_substitution(const rmatrix& A, const rvector& b);
rvector get_result(rmatrix& A, rvector& b);
bool gauss_elimination(rmatrix& A, rvector& b);
int find_max(const rmatrix& A, int k);
void check(const rmatrix& A, const rvector& x);
void print(const rmatrix& A, const rvector& b);
void print(const rvector& b);
void test();


// Маємо рівняння Ax=b
// Тут A - це nxn матриця коефіцієнтів, b - це n-елементний вектор
// Дана функція повертає x - n-елементний вектор розв'язків
// Якщо A - вироджена, то повертає вектор нульового розміру
rvector get_result(rmatrix& A, rvector& b) {
	if (gauss_elimination(A, b)) return rvector(0);
	return back_substitution(A, b);
}


// Функція виконує зворотній процес Гауса
rvector back_substitution(const rmatrix& A, const rvector& b) {
	// Створюємо вектор розв'язків
	int n = A.size();
	rvector x(n); 

	// Обчислюємо х від x[n-1] до x[0]
	for (int i = n - 1; i >= 0; --i) {
		// Змінна s для збереження сум добутків і-тих значень з матриці А на відповідні ікси
		double s = 0;
		for (int j = i + 1; j < n; ++j) s += A[i][j] * x[j]; 
		// Обчислюємо ікси
		x[i] = (b[i] - s) / A[i][i];
	}
	return x;
}


// Функція, що знаходить макс елементи для постовпцевого вибору
int find_max(const rmatrix& A, int k){
	int n = A.size();
	int imax = k; // Індекс рядка для знаходження макс елемента
	double max_pivot = abs(A[k][k]);

	for (int i = k + 1; i < n; ++i) {
		double a = abs(A[i][k]);
		if (a > max_pivot) {
			max_pivot = a;
			imax = i;
		}
	}
	return imax;
}


// Прямий хід методу Гауса
bool gauss_elimination(rmatrix& A, rvector& b) {
	int n = A.size();
	double det = 1; // Визначник
	int swap_count=0;
	
	cout << "A before elimination" << endl;
	print(A, b);

	for (int k = 0; k < n; ++k) {
		int imax = find_max(A, k); // Знаходимо номер рядка з найбільшим по модулю елементом
		if (A[imax][k] == 0) {
			cout << endl << "Det A = 0" << endl;
			return true;
		} // Якщо матриця вироджена
		swap(A[k], A[imax]); swap(b[k], b[imax]); // Перестановка відповідних елементів в матриці А та векторі b
		if (imax != k) ++swap_count; 
			
		for (int i = k + 1; i < n; ++i) {
			double c = -A[i][k] / A[k][k]; // Коефіцієнт для зведення  

			for (int j = k; j < n; ++j) {
				A[i][j] += c * A[k][j];
			} // Зводимо елементи рядка матриці А
			b[i] += c * b[k]; // Зводимо відповідний елемент b[i]
		}
		det *= A[k][k]; // Обчислюємо визначник
	}

	if (!A.empty()) {
		cout << endl << "A after elimination" << endl;
		print(A, b);
	}
	else cout << "A is empty" << endl;
	
	if (swap_count % 2) det *= -1;
	cout << endl << "Det A = " << det << endl;
	return false; // Якщо невироджена
}


// Друкуємо матрицю А і вектор b  
void print(const rmatrix& A, const rvector& b) {
	cout.precision(10); // 10 значущих знаків після коми
	for (size_t i = 0; i < A.size(); ++i) {
		for (size_t j = 0; j < A.size(); ++j) {
			cout << left << setw(14) << A[i][j]; // Форматуємо вивід, вирівнюємо вліво
		}
		cout << left << '|' << " " << b[i] << endl;
	}
}


// Друкуємо вектор
void print(const rvector& b) {
	cout.precision(10);
	for (const auto& i : b) {
		cout << i << ' ';
	}
	cout << endl;
}


// Перевірка правильності нашої реалізації методу Гауса
void check(const rmatrix& A, const rvector& x) {
	double val;
	for (size_t i = 0; i < A.size(); ++i) {
		val = 0;
		for (size_t j = 0; j < A.size(); ++j) {
			val += A[i][j]*x[j];
		}
		cout.precision(10);
		cout << val << ' ';
	}
	cout << endl;
}


// Функція для вводу з файла
void file_input(ifstream& in, rmatrix& A, rvector& b) {
	if (!in.is_open()) { 
		cout << "Error: File was not open" << endl;
		return; 
	}
	// В файлі повинні бути дані в такому форматі
	// Число n-розмір матриці
	// (n+1)*n чисел: кожні n-1 чисел матриці A і 1 число вектора b
	// a b c --- a, b -> A, c->b
	// d e f --- d, e -> A, e->b
	int temp;
	in >> temp;
	A.resize(temp);
	for (auto& i : A) i.resize(temp);
	b.resize(temp);
	
	for (size_t i = 0; i < A.size(); ++i) {
		for (size_t j = 0; j < A.size(); ++j) {
			in >> A[i][j];
		}
		in >> b[i];
	}
}


// Тестуємо нашу реалізацію методу Гауса
void test() {
	// Простестуємо: спочатку ввід з файлу 
	ifstream in("data.txt"); rmatrix a; rvector b;
	for (int i = 0; i < 4; ++i) {
		cout << "=======================================Test " << i + 1 << "=======================================" << endl;
		file_input(in, a, b);
		rvector x = get_result(a, b);
		if (!x.empty()) {
			cout << endl << "Vector of x" << endl;
			print(x);
			cout << endl << "Check our results (Ax=b)" << endl;
			check(a, x);
		}
		else if (x.empty() && in.is_open()) {
			cout << "A is singular, it has no solutions" << endl;
		}
		cout << endl << "Press ENTER to continue" << endl;
		cin.get();
		cout << endl << endl;
	}
}


int main()
{
	test();


	system("pause");
	return 0;
}
