#include<iostream>
#include <fstream>
#include <cmath>
#include <time.h>
#include <sstream>
#include <string>
#include <vector>
#include "Header.h"
#include <chrono>
using namespace std::chrono;
using namespace std;

template <class T>
void hoanVi(T& a, T& b)
{
	T x = a;
	a = b;
	b = x;
}

//-------------------------------------------------

// Hàm phát sinh mảng dữ liệu ngẫu nhiên
void generateRandomData(int a[], int n)
{
	srand((unsigned int)time(NULL));

	for (int i = 0; i < n; i++)
	{
		a[i] = rand() % n;
	}
}

// Hàm phát sinh mảng dữ liệu có thứ tự tăng dần
void generateSortedData(int a[], int n)
{
	for (int i = 0; i < n; i++)
	{
		a[i] = i;
	}
}

// Hàm phát sinh mảng dữ liệu có thứ tự ngược (giảm dần)
void generateReverseData(int a[], int n)
{
	for (int i = 0; i < n; i++)
	{
		a[i] = n - 1 - i;
	}
}

// Hàm phát sinh mảng dữ liệu gần như có thứ tự
void generateNearlySortedData(int a[], int n)
{
	for (int i = 0; i < n; i++)
	{
		a[i] = i;
	}
	srand((unsigned int)time(NULL));
	for (int i = 0; i < 10; i++)
	{
		int r1 = rand() % n;
		int r2 = rand() % n;
		hoanVi(a[r1], a[r2]);
	}
}

void generateData(int a[], int n, int dataType)
{
	switch (dataType)
	{
	case 0:	// ngẫu nhiên
		generateRandomData(a, n);
		break;
	case 1:	// có thứ tự
		generateSortedData(a, n);
		break;
	case 2:	// có thứ tự ngược
		generateReverseData(a, n);
		break;
	case 3:	// gần như có thứ tự
		generateNearlySortedData(a, n);
		break;
	default:
		printf("Error: unknown data type!\n");
	}
}



void Swap(int& a, int& b) {
	int temp = a;
	a = b;
	b = temp;
}

//Selection Sort

void selectionSort(int* a, int n) {

	int minIndex;

	for (int i = 0; i < n - 1; i++)
	{
		minIndex = i;
		for (int j = i + 1; j < n; j++)
		{
			if (a[j] < a[minIndex])
				minIndex = j;
		}

		Swap(a[i], a[minIndex]);
	}


}
unsigned long long countSelectionSort(int* a, int n) {

	unsigned long long comparisons = 0;

	int minIndex;
	for (int i = 0; (++comparisons) && (i < n - 1); i++)
	{
		minIndex = i;
		for (int j = i + 1; (++comparisons) && (j < n); j++)
		{
			if ((++comparisons) && (a[j] < a[minIndex]))
				minIndex = j;
		}

		Swap(a[i], a[minIndex]);
	}

	return comparisons;
}

long double timeSelectionSort(int* a, int n) {

	auto t1 = high_resolution_clock::now();
	selectionSort(a, n);
	auto t2 = high_resolution_clock::now();

	duration<long double, std::milli> ms_double = t2 - t1;

	return ms_double.count();
}

//Quick Sort

int Partition(int* a, int left, int right) {
	
	int randomIndex = left + rand() % (right - left + 1);
	Swap(a[left], a[randomIndex]); 
	int pivot = a[left];

	
	int count = 0;
	for (int i = left + 1; i <= right; i++) {
		if (a[i] <= pivot)
			count++;
	}

	
	int pivotIndex = left + count;
	Swap(a[pivotIndex], a[left]);

	
	int i = left, j = right;
	while (i < pivotIndex && j > pivotIndex) {
		while (a[i] <= pivot) i++;
		while (a[j] > pivot) j--;
		if (i < pivotIndex && j > pivotIndex)
			Swap(a[i++], a[j--]);
	}

	return pivotIndex;
}

void quickSort(int* a, int left, int right) {
	if (left < right)
	{
		int p = Partition(a, left, right);
		quickSort(a, left, p - 1);
		quickSort(a, p + 1, right);
	}
}

int countPartition(int* a, int left, int right, unsigned long long& countCompQuickSort) {
	int randomIndex = left + rand() % (right - left + 1);
	Swap(a[left], a[randomIndex]);
	int pivot = a[left];

	int count = 0;
	for (int i = left + 1; (++countCompQuickSort) && (i <= right); i++) {
		if ((++countCompQuickSort) && (a[i] <= pivot))
			count++;
	}

	int pivotIndex = left + count;
	Swap(a[pivotIndex], a[left]);

	int i = left, j = right;
	while (((++countCompQuickSort) && (i < pivotIndex)) && ((++countCompQuickSort) && (j > pivotIndex))) {
		while ((++countCompQuickSort) && (a[i] <= pivot)) i++;
		while ((++countCompQuickSort) && (a[j] > pivot)) j--;
		if (((++countCompQuickSort) && (i < pivotIndex)) && ((++countCompQuickSort) && (j > pivotIndex)))
			Swap(a[i++], a[j--]);
	}

	return pivotIndex;
}

unsigned long long countQuickSort(int* a, int left, int right) {
	unsigned long long countCompQuickSort = 0;
	if (left < right) {
		int p = countPartition(a, left, right, countCompQuickSort);
		countCompQuickSort += countQuickSort(a, left, p - 1);
		countCompQuickSort += countQuickSort(a, p + 1, right);
	}
	return countCompQuickSort;
}

long double timeQuickSort(int* a, int left, int right) {

	auto t1 = high_resolution_clock::now();
	quickSort(a, left, right);
	auto t2 = high_resolution_clock::now();

	duration<long double, std::milli> ms_double = t2 - t1;

	return ms_double.count();

}

//----------------------------------------------------------------------

// --------------------INSERTION SORT--------------------
unsigned long long insertionSortCMP(int* a, int n)
{
	unsigned long long cmp = 0;
	int key;
	int j;
	for (int i = 1; ++cmp && i < n; i++)
	{
		key = a[i];
		j = i - 1;

		while (++cmp && j >= 0 && ++cmp && a[j] > key)
		{
			a[j + 1] = a[j];
			j = j - 1;
		}
		a[j + 1] = key;
	}

	return cmp;
}
double insertionSortRT(int* a, int n)
{
	// get time before sort
	auto t1 = high_resolution_clock::now();

	// sort
	int key;
	int j;
	for (int i = 1; i < n; i++)
	{
		key = a[i];
		j = i - 1;

		while (j >= 0 && a[j] > key)
		{
			a[j + 1] = a[j];
			j = j - 1;
		}
		a[j + 1] = key;
	}

	//get time after sorted
	auto t2 = high_resolution_clock::now();
	// caculate time (ms)
	std::chrono::duration<double, std::milli> runTime = t2 - t1;

	return runTime.count();
}

// ----------------------MERGE SORT----------------------
unsigned long long mergeCMP(int* a, int left, int right)
{
	unsigned long long cmp = 0;
	int mid = left + (right - left) / 2;



	// divide into tempL and tempR
	// ps: tempL contains mid
	int sizeL = mid - left + 1;
	int sizeR = right - mid;
	int* tempL = new int[sizeL];
	int* tempR = new int[sizeR];
	for (int i = 0; ++cmp && i < sizeL; i++)
		tempL[i] = a[left + i];
	for (int j = 0; ++cmp && j < sizeR; j++)
		tempR[j] = a[mid + 1 + j];



	// use i,j,k for 3 while below
	int i = 0; // index of tempL
	int j = 0; // index of tempR
	int k = left; // index of merged arr




	// Merge tempL tempR into arr
	while (++cmp && i < sizeL && ++cmp && j < sizeR) {
		if (++cmp && tempL[i] <= tempR[j]) {
			a[k] = tempL[i];
			i++;
		}
		else {
			a[k] = tempR[j];
			j++;
		}
		k++;
	}
	// if sizeL != sizeR -> copy remainings to arr
	while (++cmp && i < sizeL) {
		a[k] = tempL[i];
		i++;
		k++;
	}
	while (++cmp && j < sizeR) {
		a[k] = tempR[j];
		j++;
		k++;
	}




	delete[] tempL;
	delete[] tempR;
	return cmp;
}
int mergeSortCMP(int* a, int begin, int end)
{
	int cmp = 0;

	// stop case
	if (++cmp && begin >= end)
		return cmp;

	int mid = begin + (end - begin) / 2;
	cmp += mergeSortCMP(a, begin, mid);
	cmp += mergeSortCMP(a, mid + 1, end);
	cmp += mergeCMP(a, begin, end);

	return cmp;
}

void merge(int* a, int left, int right)
{

	// divide into tempL and tempR
	// ps: tempL contains mid
	int mid = left + (right - left) / 2;
	int sizeL = mid - left + 1;
	int sizeR = right - mid;
	int* tempL = new int[sizeL];
	int* tempR = new int[sizeR];
	for (int i = 0; i < sizeL; i++)
		tempL[i] = a[left + i];
	for (int j = 0; j < sizeR; j++)
		tempR[j] = a[mid + 1 + j];



	// use i,j,k for 3 while below
	int i = 0; // index of tempL
	int j = 0; // index of tempR
	int k = left; // index of merged arr




	// Merge tempL tempR into arr
	while (i < sizeL && j < sizeR) {
		if (tempL[i] <= tempR[j]) {
			a[k] = tempL[i];
			i++;
		}
		else {
			a[k] = tempR[j];
			j++;
		}
		k++;
	}
	// if sizeL != sizeR -> copy remainings to arr
	while (i < sizeL) {
		a[k] = tempL[i];
		i++;
		k++;
	}
	while (j < sizeR) {
		a[k] = tempR[j];
		j++;
		k++;
	}



	delete[] tempL;
	delete[] tempR;
}
void mergeSort(int* a, int begin, int end)
{
	// stop case
	if (begin >= end)
		return;

	int mid = begin + (end - begin) / 2;
	mergeSort(a, begin, mid);
	mergeSort(a, mid + 1, end);
	merge(a, begin, end);
}
double mergeSortRT(int* a, int n)
{
	int begin = 0;
	int end = begin + n - 1;

	// get time before sort
	auto t1 = high_resolution_clock::now();

	// stop case
	if (begin >= end)
	{
		//get time after sort
		auto t2 = high_resolution_clock::now();
		// caculate time (ms)
		std::chrono::duration<double, std::milli> runTime = t2 - t1;

		return runTime.count();
	}

	mergeSort(a, begin, end);

	//get time after sort
	auto t2 = high_resolution_clock::now();
	// caculate time (ms)
	std::chrono::duration<double, std::milli> runTime = t2 - t1;

	return runTime.count();
}

//-----------------------------------
int* Copy_Array(int* arr, int n)
{
	int* temp = new int[n];
	for (int i = 0; i < n; i++)
	{
		temp[i] = arr[i];
	}
	return temp;
}

void Heap(int* arr, int N, int i)
{
	int langest = i;
	int l = 2 * i + 1;
	int r = 2 * i + 2;
	if (l<N && arr[l]>arr[i])
		langest = l;
	if (r<N && arr[r]>arr[i])
	{
		if (arr[r] > arr[l])
			langest = r;
	}
	if (langest != i)
	{
		hoanVi(arr[i], arr[langest]);
		Heap(arr, N, langest);
	}
}

//runtime heap_sort
double Heap_sort(int* arr, int n)
{
	auto start = chrono::steady_clock::now();

	for (int i = n / 2 - 1; i >= 0; i--)
	{
		Heap(arr, n, i);
	}
	for (int i = n - 1; i > 0; i--)
	{
		hoanVi(arr[i], arr[0]);
		Heap(arr, i, 0);
	}
	auto end = chrono::steady_clock::now();
	auto time_run = (end - start);
	return chrono::duration<double, milli>(time_run).count();
}
// Tinh phep so sanh
void heap_Comparison(int* arr, int N, int i, unsigned long long& cmpr)
{
	int langest = i;
	int l = 2 * i + 1;
	int r = 2 * i + 2;
	if ((++cmpr && l < N) && (++cmpr && arr[l] > arr[i]))
		langest = l;
	if ((++cmpr && r < N) && (++cmpr && arr[r] > arr[i]))
	{
		if (++cmpr && arr[r] > arr[l])
			langest = r;
	}
	if (++cmpr && langest != i)
	{
		hoanVi(arr[i], arr[langest]);
		heap_Comparison(arr, N, langest, cmpr);
	}
}
unsigned long long   Heap_sort_Comparison(int* arr, int n)
{
	unsigned long long  cmpr = 0;
	for (int i = n / 2 - 1; ++cmpr && i >= 0; i--)
	{
		heap_Comparison(arr, n, i, cmpr);
	}
	for (int i = n - 1; ++cmpr && i > 0; i--)
	{
		hoanVi(arr[i], arr[0]);
		heap_Comparison(arr, i, 0, cmpr);
	}
	return cmpr;
}
// Bubble Sort
double BubbleSort(int* arr, int n)
{
	auto start = chrono::steady_clock::now();
	int i, j;
	for (i = 0; i < n - 1; i++)
	{

		for (j = 0; j < n - i - 1; j++)
		{
			if (arr[j] > arr[j + 1])
				hoanVi(arr[j], arr[j + 1]);
		}
	}
	auto end = chrono::steady_clock::now();
	auto time_run = (end - start);

	return chrono::duration<double, milli>(time_run).count();
}

// Tinh Phep So Sanh
unsigned long long  BubbleSort_Comparison(int* arr, int n)
{
	unsigned long long   cmpr = 0;
	int i, j;
	for (i = 0; (++cmpr && i < n - 1); i++)
	{
		for (j = 0; (++cmpr && j < n - i - 1); j++)
		{
			if (++cmpr && arr[j] > arr[j + 1])
				hoanVi(arr[j], arr[j + 1]);
		}
	}
	return cmpr;
}

//--------------------
int getMax(int arr[], int n)
{
	int max = arr[0];
	for (int i = 1; i < n; i++)
		if (arr[i] > max)
			max = arr[i];
	return max;
}
void countSort(int* arr, int n, int place)
{

	int* output = new int[n];
	int i, count[10] = { 0 };


	for (i = 0; i < n; i++)
		count[(arr[i] / place) % 10]++;


	for (i = 1; i < 10; i++)
		count[i] += count[i - 1];


	for (i = n - 1; i >= 0; i--) {
		output[count[(arr[i] / place) % 10] - 1] = arr[i];
		count[(arr[i] / place) % 10]--;
	}


	for (i = 0; i < n; i++)
		arr[i] = output[i];
}
void countSortToCountComparisons(int* arr, int n, int place, int& Comparisons)
{

	int* output = new int[n];
	int i, count[10] = { 0 };


	for (i = 0; (++Comparisons) && (i < n); i++)
		count[(arr[i] / place) % 10]++;


	for (i = 1; (++Comparisons) && (i < 10); i++)
		count[i] += count[i - 1];


	for (i = n - 1; (++Comparisons) && (i >= 0); i--) {
		output[count[(arr[i] / place) % 10] - 1] = arr[i];
		count[(arr[i] / place) % 10]--;
	}


	for (i = 0; (++Comparisons) && (i < n); i++)
		arr[i] = output[i];
}
void radixSort(int* arr, int n)
{
	int max = getMax(arr, n);
	for (int place = 1; max / place > 0; place *= 10)
		countSort(arr, n, place);
}
double getTime(int* arr, int n) {
	auto start = std::chrono::high_resolution_clock::now();
	radixSort(arr, n);
	auto finish = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double, std::milli> timeusing = finish - start;

	return timeusing.count();
}
int Counting(int* arr, int n) {
	int Comparison = 0;

	int m = getMax(arr, n);
	for (int place = 1; (++Comparison) && (m / place > 0); place *= 10)
		countSortToCountComparisons(arr, n, place, Comparison);
	return Comparison;

}
// -------------Shell sort-------------------
unsigned long long shellSortCMP(int* a, int n)
{
	unsigned long long cmp = 0;

	// Dynamic gap sequence
	int gap = n / 2;
	while (++cmp && gap > 0)
	{
		for (int i = gap; ++cmp && i < n; i++)
		{
			int temp = a[i];
			int j = i;
			while (++cmp && j >= gap && ++cmp && a[j - gap] > temp)
			{
				a[j] = a[j - gap];
				j -= gap;
			}
			a[j] = temp;
		}
		gap /= 2;
	}
	return cmp;
}

double shellSortRT(int* a, int n)
{
	auto t1 = std::chrono::high_resolution_clock::now();

	// Perform Shell Sort
	int gap = n / 2;
	while (gap > 0)
	{
		for (int i = gap; i < n; i++)
		{
			int temp = a[i];
			int j = i;
			while (j >= gap && a[j - gap] > temp)
			{
				a[j] = a[j - gap];
				j -= gap;
			}
			a[j] = temp;
		}
		gap /= 2;
	}

	auto t2 = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double, std::milli> runTime = t2 - t1;

	return runTime.count();

}
void ShellSort(int* a, int n)
{
	int gap = n / 2;
	while (gap > 0)
	{
		for (int i = gap; i < n; i++)
		{
			int temp = a[i];
			int j = i;
			while (j >= gap && a[j - gap] > temp)
			{
				a[j] = a[j - gap];
				j -= gap;
			}
			a[j] = temp;
		}
		gap /= 2;
	}
}
// -------------------- Count_Sort-------------------
unsigned long long countingSortCMP(int* arr, int n) {
	unsigned long long comparisons = 0;

	int maxElement = *max_element(arr, arr + n);
	int minElement = *min_element(arr, arr + n);
	int range = maxElement - minElement + 1;

	vector<int> count(range, 0);
	vector<int> output(n);

	for (int i = 0; (++comparisons) && i < n; ++i)
		count[arr[i] - minElement]++;

	for (int i = 1; (++comparisons) && i < range; ++i)
		count[i] += count[i - 1];

	for (int i = n - 1; (++comparisons) && i >= 0; --i) {
		output[count[arr[i] - minElement] - 1] = arr[i];
		count[arr[i] - minElement]--;
	}

	return comparisons;
}


double countingSortRT(int* arr, int n) {
    auto start = chrono::high_resolution_clock::now();

    int maxElement = *max_element(arr, arr + n);
    int minElement = *min_element(arr, arr + n);
    int range = maxElement - minElement + 1;

    vector<int> count(range, 0);
    vector<int> output(n);

    for (int i = 0; i < n; ++i)
        count[arr[i] - minElement]++;

    for (int i = 1; i < range; ++i)
        count[i] += count[i - 1];

    for (int i = n - 1; i >= 0; --i) {
        output[count[arr[i] - minElement] - 1] = arr[i];
        count[arr[i] - minElement]--;
    }

    for (int i = 0; i < n; ++i)
        arr[i] = output[i];

    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double, milli> time = end - start;
    return time.count();
}
void Print_file_Input(string File_in, int* arr, int n)
{
	ofstream Fout;
	Fout.open(File_in);
	Fout << n << endl;
	for (int i = 0; i < n; i++)
	{
		Fout << arr[i] << " ";
	}
	Fout.close();
}
char Xac_Dinh_Thuat_Toan(char* algorithm)
{
	char* k = NULL;
	k = strstr(algorithm, "selection");
	if (k != NULL)
	{
		return 's';
	}
	k = strstr(algorithm, "insertion");
	if (k != NULL)
	{
		return 'i';
	}
	k = strstr(algorithm, "bubble");
	if (k != NULL)
	{
		return 'b';
	}
	k = strstr(algorithm, "heap");
	if (k != NULL)
	{
		return 'h';
	}
	k = strstr(algorithm, "merge");
	if (k != NULL)
	{
		return 'm';
	}
	k = strstr(algorithm, "quick");
	if (k != NULL)
	{
		return 'q';
	}
	k = strstr(algorithm, "radix");
	if (k != NULL)
	{
		return 'r';
	}
	k= strstr(algorithm, "shell");
	if (k != NULL)
	{
		return 'sh';
	}
	k = strstr(algorithm, "count");
	if (k != NULL)
	{
		return 'c';
	}
	return 'e';
}
int Xac_dinh_parameter(char* para)
{
	switch (para[1])
	{
	case 't':
	{
		return 1;
		break;
	}
	case 'c':
	{
		return 2;
		break;
	}
	case 'b':
	{
		return 3;
		break;
	}
	}
	return 0;
}
void Cout_consolve(char algorithm, int para, int* arr, int n)
{
	switch (algorithm)
	{
	case 's':
	{
		switch (para)
		{
		case 1:
		{
			cout << "Running Time: " << timeSelectionSort(arr, n) << " ms\n";
			break;
		}
		case 2:
		{
			cout << "Comparisions: " << countSelectionSort(arr, n) << endl;
			break;
		}
		case 3:
		{
			int* arr2 = Copy_Array(arr, n);

			cout << "Running Time: " << timeSelectionSort(arr, n) << " ms\n";
			cout << "Comparisions: " << countSelectionSort(arr2, n) << endl;

			delete[] arr2;
			break;
		}
		}
		break;
	}
	case 'i':
	{
		switch (para)
		{
		case 1:
		{
			cout << "Running Time: " << insertionSortRT(arr, n) << " ms\n";
			break;
		}
		case 2:
		{
			cout << "Comparisions: " << insertionSortCMP(arr, n) << endl;
			break;
		}
		case 3:
		{
			int* arr2 = Copy_Array(arr, n);

			cout << "Running Time: " << insertionSortRT(arr, n) << " ms\n";
			cout << "Comparisions: " << insertionSortCMP(arr2, n) << endl;

			delete[] arr2;
			break;
		}
		}
		break;
	}
	case 'b':
	{
		switch (para)
		{
		case 1:
		{
			cout << "Running Time: " << BubbleSort(arr, n) << " ms\n";
			break;
		}
		case 2:
		{
			cout << "Comparisions: " << BubbleSort_Comparison(arr, n) << endl;
			break;
		}
		case 3:
		{
			int* arr2 = Copy_Array(arr, n);

			cout << "Running Time: " << BubbleSort(arr, n) << " ms\n";
			cout << "Comparisions: " << BubbleSort_Comparison(arr2, n) << endl;

			delete[] arr2;
			break;
		}
		}
		break;
	}
	case 'h':
	{
		switch (para)
		{
		case 1:
		{
			cout << "Running Time: " << Heap_sort(arr, n) << " ms\n";
			break;
		}
		case 2:
		{
			cout << "Comparisions: " << Heap_sort_Comparison(arr, n) << endl;
			break;
		}
		case 3:
		{
			int* arr2 = Copy_Array(arr, n);

			cout << "Running Time: " << Heap_sort(arr, n) << " ms\n";
			cout << "Comparisions: " << Heap_sort_Comparison(arr2, n) << endl;

			delete[] arr2;
			break;
		}
		}
		break;
	}
	case 'm':
	{
		switch (para)
		{
		case 1:
		{
			cout << "Running Time: " << mergeSortRT(arr, n) << " ms\n";
			break;
		}
		case 2:
		{
			cout << "Comparisions: " << mergeSortCMP(arr, 0, n - 1) << endl;
			break;
		}
		case 3:
		{
			int* arr2 = Copy_Array(arr, n);

			cout << "Running Time: " << mergeSortRT(arr, n) << " ms\n";
			cout << "Comparisions: " << mergeSortCMP(arr2, 0, n - 1) << endl;

			delete[] arr2;
			break;
		}
		}
		break;
	}
	case 'q':
	{
		switch (para)
		{
		case 1:
		{
			cout << "Running Time: " << timeQuickSort(arr, 0, n - 1) << " ms\n";
			break;
		}
		case 2:
		{
			cout << "Comparisions: " << countQuickSort(arr, 0, n - 1) << endl;
			break;
		}
		case 3:
		{
			int* arr2 = Copy_Array(arr, n);

			cout << "Running Time: " << timeQuickSort(arr, 0, n - 1) << " ms\n";
			cout << "Comparisions: " << countQuickSort(arr2, 0, n - 1) << endl;

			delete[] arr2;
			break;
		}
		}
		break;
	}
	case 'r':
	{
		switch (para)
		{
		case 1:
		{
			cout << "Running Time: " << getTime(arr, n) << " ms\n";
			break;
		}
		case 2:
		{
			cout << "Comparisions: " << Counting(arr, n) << endl;
			break;
		}
		case 3:
		{
			int* arr2 = Copy_Array(arr, n);

			cout << "Running Time: " << getTime(arr, n) << " ms\n";
			cout << "Comparisions: " << Counting(arr2, n) << endl;

			delete[] arr2;
			break;
		}
		}
		break;
	}
	case 'sh':
	{
		switch (para)
		{
		case 1:
			{
			cout << "Running Time: " << shellSortRT(arr, n) << "ms " << endl;
			break;
			}
		case 2:
		{
			cout << "Comparisions: " << shellSortCMP( arr,  n)<< endl;
		}
		case 3:
		{
			int* arr2 = Copy_Array(arr, n);

			cout << "Running Time: " << shellSortRT(arr, n) << " ms\n";
			cout << "Comparisions: " << shellSortCMP(arr2, n) << endl;
			delete[] arr2;
			break;
		}
		}

	}
	case 'c':
	{
		switch (para)
		{
		case 1:
		{
			cout << "Running Time: " << countingSortRT(arr, n) << "ms " << endl;
			break;
		}
		case 2:
		{
			cout << "Comparisions: " << countingSortCMP(arr, n) << endl;
		}
		case 3:
		{
			int* arr2 = Copy_Array(arr, n);

			cout << "Running Time: " << countingSortRT(arr, n) << " ms\n";
			cout << "Comparisions: " << countingSortCMP(arr, n) << endl;
			delete[] arr2;
			break;
		}
		}
	}
	case 'e':
	{
		cout << "false";
		return;
	}
	}
}

//--------------------
int strcompare(const string a, const string b)
{
	int s = a.size();
	if (b.size() != s)
		return 0;
	for (int i = 0; i < s; ++i)
		if (tolower(a[i]) != tolower(b[i]))
			return 0;
	return 1;
}
void printInfo1(string ss1, string ss2, string ss3, int* a, int n) {
	string selection_sort = "selection-sort", quick_sort = "quick-sort", insertion_sort = "insertion-sort", merge_sort = "merge-sort",
		heap_sort = "heap-sort", bubble_sort = "bubble-sort", radix_sort = "radix-sort", shell_sort = "shell-sort", count_sort="count-sort";
	string time = "-time", comparsion = "-comp", both = "-both";

	if (strcompare(ss1, selection_sort) == 1) {
		cout << "Algorithm: Selection Sort" << endl;
		cout << "Input file: " << ss2 << endl;
		cout << "Input size: " << n << endl;
		cout << "----------------------------------------" << endl;
		if (strcompare(ss3, time)) {
			cout << "Running time: " << timeSelectionSort(a, n) << " ms" << endl;
		}
		else if (strcompare(ss3, comparsion)) {
			cout << "Comparisons: " << countSelectionSort(a, n) << endl;
		}
		else if (strcompare(ss3, both)) {
			cout << "Running time: " << timeSelectionSort(a, n) << " ms" << endl;
			cout << "Comparisons: " << countSelectionSort(a, n) << endl;
		}
		else cout << "Wrong input!!" << endl;

	}
	else if (strcompare(ss1, insertion_sort) == 1) {
		cout << "Algorithm: Insertion Sort" << endl;
		cout << "Input file: " << ss2 << endl;
		cout << "Input size: " << n << endl;
		cout << "----------------------------------------" << endl;
		if (strcompare(ss3, time)) {
			cout << "Running time: " << insertionSortRT(a, n) << " ms" << endl;
		}
		else if (strcompare(ss3, comparsion)) {
			cout << "Comparisons: " << insertionSortCMP(a, n) << endl;
		}
		else if (strcompare(ss3, both)) {
			cout << "Running time: " << insertionSortRT(a, n) << " ms" << endl;
			cout << "Comparisons: " << insertionSortCMP(a, n) << endl;
		}
		else cout << "Wrong input!!" << endl;
	}
	else if (strcompare(ss1, merge_sort) == 1) {
		cout << "Algorithm: Merge Sort" << endl;
		cout << "Input file: " << ss2 << endl;
		cout << "Input size: " << n << endl;
		cout << "----------------------------------------" << endl;
		if (strcompare(ss3, time)) {
			cout << "Running time: " << mergeSortRT(a, n) << " ms" << endl;
		}
		else if (strcompare(ss3, comparsion)) {
			cout << "Comparisons: " << mergeSortCMP(a, 0, n - 1) << endl;
		}
		else if (strcompare(ss3, both)) {
			cout << "Running time: " << mergeSortRT(a, n) << " ms" << endl;
			cout << "Comparisons: " << mergeSortCMP(a, 0, n - 1) << endl;
		}
		else cout << "Wrong input!!" << endl;
	}
	else if (strcompare(ss1, heap_sort) == 1) {
		cout << "Algorithm: Heap Sort" << endl;
		cout << "Input file: " << ss3 << endl;
		cout << "Input size: " << n << endl;
		cout << "----------------------------------------" << endl;
		if (strcompare(ss2, time)) {
			cout << "Running time: " << Heap_sort(a, n) << " ms" << endl;
		}
		else if (strcompare(ss2, comparsion)) {
			cout << "Comparisons: " << Heap_sort_Comparison(a, n) << endl;
		}
		else if (strcompare(ss2, both)) {
			cout << "Running time: " << Heap_sort(a, n) << " ms" << endl;
			cout << "Comparisons: " << Heap_sort_Comparison(a, n) << endl;
		}
		else cout << "Wrong input!!" << endl;
	}
	else if (strcompare(ss1, bubble_sort) == 1) {
		cout << "Algorithm: Bubble Sort" << endl;
		cout << "Input file: " << ss2 << endl;
		cout << "Input size: " << n << endl;
		cout << "----------------------------------------" << endl;
		if (strcompare(ss3, time)) {
			cout << "Running time: " << BubbleSort(a, n) << " ms" << endl;
		}
		else if (strcompare(ss3, comparsion)) {
			cout << "Comparisons: " << BubbleSort_Comparison(a, n) << endl;
		}
		else if (strcompare(ss3, both)) {
			cout << "Running time: " << BubbleSort(a, n) << " ms" << endl;
			cout << "Comparisons: " << BubbleSort_Comparison(a, n) << endl;
		}
		else cout << "Wrong input!!" << endl;
	}
	else if (strcompare(ss1, radix_sort) == 1) {
		cout << "Algorithm: Radix Sort" << endl;
		cout << "Input file: " << ss2 << endl;
		cout << "Input size: " << n << endl;
		cout << "----------------------------------------" << endl;
		if (strcompare(ss3, time)) {
			cout << "Running time: " << getTime(a, n) << " ms" << endl;
		}
		else if (strcompare(ss3, comparsion)) {
			cout << "Comparisons: " << Counting(a, n) << endl;
		}
		else if (strcompare(ss3, both)) {
			cout << "Running time: " << getTime(a, n) << " ms" << endl;
			cout << "Comparisons: " << Counting(a, n) << endl;
		}
		else cout << "Wrong input!!" << endl;
	}
	else if (strcompare(ss1, quick_sort) == 1) {
		cout << "Algorithm: Quick Sort" << endl;
		cout << "Input file: " << ss2 << endl;
		cout << "Input size: " << n << endl;
		cout << "----------------------------------------" << endl;
		if (strcompare(ss3, time)) {
			cout << "Running time: " << timeQuickSort(a, 0, n - 1) << " ms" << endl;
		}
		else if (strcompare(ss3, comparsion)) {
			cout << "Comparisons: " << countQuickSort(a, 0, n - 1) << endl;
		}
		else if (strcompare(ss3, both)) {
			cout << "Running time: " << timeQuickSort(a, 0, n - 1) << " ms" << endl;
			cout << "Comparisons: " << countQuickSort(a, 0, n - 1) << endl;
		}
		else cout << "Wrong input!!" << endl;
	}
	else if (strcompare(ss1, shell_sort) == 1)
	{
		cout << "Algorithm: Shell sort" << endl;
		cout << "Input file: " << ss2 << endl;
		cout << "Input size: " << n << endl;
		cout << "----------------------------------------" << endl;
		if (strcompare(ss3, time)) {
			cout << "Running time: " << shellSortRT(a, n ) << " ms" << endl;
		}
		else if (strcompare(ss3, comparsion)) {
			cout << "Comparisons: " << shellSortCMP(a, n ) << endl;
		}
		else if (strcompare(ss3, both)) {
			cout << "Running time: " << shellSortRT(a, n) << " ms" << endl;
			cout << "Comparisons: " << shellSortCMP(a, n) << endl;
		}
		else cout << "Wrong input!!" << endl;
	}
	else if (strcompare(ss1, count_sort) == 1) {
		cout << "Algorithm: Counting Sort" << endl;
		cout << "Input file: " << ss2 << endl;
		cout << "Input size: " << n << endl;
		cout << "----------------------------------------" << endl;
		if (strcompare(ss3, time)) {
			cout << "Running time: " << countingSortRT(a, n) << " ms" << endl;
		}
		else if (strcompare(ss3, comparsion)) {
			cout << "Comparisons: " << countingSortCMP(a, n) << endl;
		}
		else if (strcompare(ss3, both)) {
			cout << "Running time: " << countingSortRT(a, n) << " ms" << endl;
			cout << "Comparisons: " << countingSortCMP(a, n) << endl;
		}
		else cout << "Wrong input!!" << endl;
		}
	else cout << "Wrong input!!" << endl;
}
int command1(int argc, char* argv[]) {


	if (argc != 5) {
		cout << "Please enter 5 arguments" << endl;
		return 0;
	}
	fstream fin(argv[3], ios::in);
	if (!fin.is_open()) {
		cout << "can not open file" << endl;
		return 0;
	}

	int elements = 0;
	fin >> elements;
	if (elements > 1000000) {
		cout << "Please enter lower number" << endl;
		return 0;
	}
	int* a = new int[elements];
	int i = 0;
	while (i < elements) {
		fin >> a[i];
		i++;
	}
	printInfo1(argv[2], argv[3], argv[4], a, elements);
	fin.close();
	fstream fout("output.txt", ios::out);
	fout << elements << endl;
	for (int i = 0; i < elements; i++) {
		fout << a[i] << " ";
	}
	fout.close();
	system("pause");
	return 1;
}
void printInfo2(string ss1, string ss2, string ss3, int* a, int n) {
	string selection_sort = "selection-sort", quick_sort = "quick-sort", insertion_sort = "insertion-sort", merge_sort = "merge-sort",
		heap_sort = "heap-sort", bubble_sort = "bubble-sort", radix_sort = "radix-sort", shell_sort = "shell-sort",count_sort="count-sort";
	string time = "-time", comparsion = "-comp", both = "-both";
	string order;
	if (strcompare(ss2, "-rand")) {
		order = "Randomized data";
	}
	else if (strcompare(ss2, "-nsorted")) {
		order = " Nearly sorted data";
	}
	else if (strcompare(ss2, "-sorted")) {
		order = " Sorted data";
	}
	else if (strcompare(ss2, "-rev")) {
		order = "Reverse sorted data";
	}

	if (strcompare(ss1, selection_sort) == 1) {
		cout << "Algorithm: Selection Sort" << endl;
		cout << "Input size: " << n << endl;
		cout << "Input order: " << order << endl;
		cout << "----------------------------------------" << endl;
		if (strcompare(ss3, time)) {
			cout << "Running time: " << timeSelectionSort(a, n) << " ms" << endl;
		}
		else if (strcompare(ss3, comparsion)) {
			cout << "Comparisons: " << countSelectionSort(a, n) << endl;
		}
		else if (strcompare(ss3, both)) {
			cout << "Running time: " << timeSelectionSort(a, n) << " ms" << endl;
			cout << "Comparisons: " << countSelectionSort(a, n) << endl;
		}
		else cout << "Wrong input!!" << endl;

	}
	else if (strcompare(ss1, insertion_sort) == 1) {
		cout << "Algorithm: Isertion Sort" << endl;
		cout << "Input size: " << n << endl;
		cout << "Input order: " << order << endl;
		cout << "----------------------------------------" << endl;
		if (strcompare(ss3, time)) {
			cout << "Running time: " << insertionSortRT(a, n) << " ms" << endl;
		}
		else if (strcompare(ss3, comparsion)) {
			cout << "Comparisons: " << insertionSortCMP(a, n) << endl;
		}
		else if (strcompare(ss3, both)) {
			cout << "Running time: " << insertionSortRT(a, n) << " ms" << endl;
			cout << "Comparisons: " << insertionSortCMP(a, n) << endl;
		}
		else cout << "Wrong input!!" << endl;
	}
	else if (strcompare(ss1, merge_sort) == 1) {
		cout << "Algorithm: Merge Sort" << endl;
		cout << "Input size: " << n << endl;
		cout << "Input order: " << order << endl;
		cout << "----------------------------------------" << endl;
		if (strcompare(ss3, time)) {
			cout << "Running time: " << mergeSortRT(a, n) << " ms" << endl;
		}
		else if (strcompare(ss3, comparsion)) {
			cout << "Comparisons: " << mergeSortCMP(a, 0, n - 1) << endl;
		}
		else if (strcompare(ss3, both)) {
			cout << "Running time: " << mergeSortRT(a, n) << " ms" << endl;
			cout << "Comparisons: " << mergeSortCMP(a, 0, n - 1) << endl;
		}
		else cout << "Wrong input!!" << endl;
	}
	else if (strcompare(ss1, heap_sort) == 1) {
		cout << "Algorithm: Heap Sort" << endl;
		cout << "Input size: " << n << endl;
		cout << "Input order: " << order << endl;
		cout << "----------------------------------------" << endl;
		if (strcompare(ss2, time)) {
			cout << "Running time: " << Heap_sort(a, n) << " ms" << endl;
		}
		else if (strcompare(ss2, comparsion)) {
			cout << "Comparisons: " << Heap_sort_Comparison(a, n) << endl;
		}
		else if (strcompare(ss2, both)) {
			cout << "Running time: " << Heap_sort(a, n) << " ms" << endl;
			cout << "Comparisons: " << Heap_sort_Comparison(a, n) << endl;
		}
		else cout << "Wrong input!!" << endl;
	}
	else if (strcompare(ss1, bubble_sort) == 1) {
		cout << "Algorithm: Bubble Sort" << endl;
		cout << "Input size: " << n << endl;
		cout << "Input order: " << order << endl;
		cout << "----------------------------------------" << endl;
		if (strcompare(ss3, time)) {
			cout << "Running time: " << BubbleSort(a, n) << " ms" << endl;
		}
		else if (strcompare(ss3, comparsion)) {
			cout << "Comparisons: " << BubbleSort_Comparison(a, n) << endl;
		}
		else if (strcompare(ss3, both)) {
			cout << "Running time: " << BubbleSort(a, n) << " ms" << endl;
			cout << "Comparisons: " << BubbleSort_Comparison(a, n) << endl;
		}
		else cout << "Wrong input!!" << endl;
	}
	else if (strcompare(ss1, radix_sort) == 1) {
		cout << "Algorithm: Radix Sort" << endl;
		cout << "Input size: " << n << endl;
		cout << "Input order: " << order << endl;
		cout << "----------------------------------------" << endl;
		if (strcompare(ss3, time)) {
			cout << "Running time: " << getTime(a, n) << " ms" << endl;
		}
		else if (strcompare(ss3, comparsion)) {
			cout << "Comparisons: " << Counting(a, n) << endl;
		}
		else if (strcompare(ss3, both)) {
			cout << "Running time: " << getTime(a, n) << " ms" << endl;
			cout << "Comparisons: " << Counting(a, n) << endl;
		}
		else cout << "Wrong input!!" << endl;
	}
	else if (strcompare(ss1, quick_sort) == 1) {
		cout << "Algorithm: Quick Sort" << endl;
		cout << "Input size: " << n << endl;
		cout << "Input order: " << order << endl;
		cout << "----------------------------------------" << endl;
		if (strcompare(ss3, time)) {
			cout << "Running time: " << timeQuickSort(a, 0, n - 1) << " ms" << endl;
		}
		else if (strcompare(ss3, comparsion)) {
			cout << "Comparisons: " << countQuickSort(a, 0, n - 1) << endl;
		}
		else if (strcompare(ss3, both)) {
			cout << "Running time: " << timeQuickSort(a, 0, n - 1) << " ms" << endl;
			cout << "Comparisons: " << countQuickSort(a, 0, n - 1) << endl;
		}
		else cout << "Wrong input!!" << endl;
	}
	else if (strcompare(ss1, shell_sort) == 1) {
		cout << "Algorithm: Shell sort" << endl;
		cout << "Input size: " << n << endl;
		cout << "Input order: " << order << endl;
		cout << "----------------------------------------" << endl;
		if (strcompare(ss3, time)) {
			cout << "Running time: " << shellSortRT(a, n) << " ms" << endl;
		}
		else if (strcompare(ss3, comparsion)) {
			cout << "Comparisons: " << shellSortCMP(a, n) << endl;
		}
		else if (strcompare(ss3, both)) {
			cout << "Running time: " << shellSortRT(a, n) << " ms" << endl;
			cout << "Comparisons: " << shellSortCMP(a, n) << endl;
		}
		else cout << "Wrong input!!" << endl;
		}
	else if (strcompare(ss1, count_sort) == 1) {
		cout << "Algorithm: Counting Sort" << endl;
		cout << "Input size: " << n << endl;
		cout << "Input order: " << order << endl;
		cout << "----------------------------------------" << endl;
		if (strcompare(ss3, time)) {
			cout << "Running time: " << countingSortRT(a, n) << " ms" << endl;
		}
		else if (strcompare(ss3, comparsion)) {
			cout << "Comparisons: " << countingSortCMP(a, n) << endl;
		}
		else if (strcompare(ss3, both)) {
			cout << "Running time: " << countingSortRT(a, n) << " ms" << endl;
			cout << "Comparisons: " << countingSortCMP(a, n) << endl;
		}
		else cout << "Wrong input!!" << endl;
		}
	else cout << "Wrong input!!" << endl;
}
int command2(int argc, char* argv[]) {


	if (argc != 6) {
		cout << "Please enter 6 arguments" << endl;
		return 0;
	}
	stringstream s(argv[3]);
	int elements = 0;
	s >> elements;
	if (elements > 1000000) {
		cout << "Please enter lower number" << endl;
		return 0;
	}
	int* a = new int[elements];

	fstream finput("input.txt", ios::out);

	if (strcompare(argv[4], "-rand")) {
		generateRandomData(a, elements);
		finput << elements << endl;
		for (int i = 0; i < elements; i++) {
			finput << a[i] << " ";
		}
		printInfo2(argv[2], argv[4], argv[5], a, elements);
	}
	else if (strcompare(argv[4], "-nsorted")) {
		generateNearlySortedData(a, elements);
		finput << elements << endl;
		for (int i = 0; i < elements; i++) {
			finput << a[i] << " ";
		}
		printInfo2(argv[2], argv[4], argv[5], a, elements);
	}
	else if (strcompare(argv[4], "-sorted")) {
		generateSortedData(a, elements);
		finput << elements << endl;
		for (int i = 0; i < elements; i++) {
			finput << a[i] << " ";
		}
		printInfo2(argv[2], argv[4], argv[5], a, elements);
	}
	else if (strcompare(argv[4], "-rev")) {
		generateReverseData(a, elements);
		finput << elements << endl;
		for (int i = 0; i < elements; i++) {
			finput << a[i] << " ";
		}
		printInfo2(argv[2], argv[4], argv[5], a, elements);
	}
	else {
		cout << "Wrong size of data" << endl;
		return 0;
	}
	finput.close();

	fstream fout("output.txt", ios::out);
	fout << elements << endl;
	for (int i = 0; i < elements; i++) {
		fout << a[i] << " ";
	}
	fout.close();
	system("pause");
	return 1;
}