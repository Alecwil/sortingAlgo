#pragma once
using namespace std;

void generateData(int a[], int n, int dataType);

unsigned long long countSelectionSort(int* a, int n);
long double timeSelectionSort(int* a, int n);

unsigned long long countQuickSort(int* a, int left, int right);
long double timeQuickSort(int* a, int left, int right);

//-----------
int mergeSortCMP(int* a, int begin, int end);
double mergeSortRT(int* a, int n);
unsigned long long insertionSortCMP(int* a, int n);
double insertionSortRT(int* a, int n);

//-------------
double Heap_sort(int* arr, int n); //tgian heapsort
unsigned long long   Heap_sort_Comparison(int* arr, int n);
double BubbleSort(int* arr, int n); //tgian bubblesort
unsigned long long  BubbleSort_Comparison(int* arr, int n);
//------------------

double countingSortRT(int* arr, int n);
unsigned long long countingSortCMP(int* arr, int n);


//------------------
double getTime(int* arr, int n); //tgian radixsort
int Counting(int* arr, int n); //phep sosanh radixsort

//------------------
void Print_file_Input(string File_in, int* arr, int n);
char Xac_Dinh_Thuat_Toan(char* algorithm);
int Xac_dinh_parameter(char* para);
void Cout_consolve(char algorithm, int para, int* arr, int n);

//-----------------
int command1(int argc, char* argv[]);
int command2(int argc, char* argv[]);
