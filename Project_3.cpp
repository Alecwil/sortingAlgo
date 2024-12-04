
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include "Header.h"
#include <string.h>
using namespace std;


int main(int argc, char* argv[]) {


	if (strcmp(argv[1], "-a") == 0) {

		cout << "ALGORITHM MODE" << endl;
		if (argc == 5 && (strstr(argv[3], ".") != NULL)) { //command 1	
			command1(argc, argv);
		}

		else if (argc == 6) { //command 2
			command2(argc, argv);
		}

		else if (argc == 5 && (strstr(argv[3], ".") == NULL)) { //command 3
			int n = atoi(argv[3]);
			char* k = NULL;
			k = strstr(argv[2], "sort");
			if (k != NULL)
			{
				k = k - 1;
				k[0] = '\0';
			}
			else
			{
				cout << "\nfalse\n";
				return 0;
			}
			cout << "Algorithm: " << argv[2] << " sort" << endl;
			cout << "Input size: " << n << endl << endl;
			int* arr = new int[n];
			char al = Xac_Dinh_Thuat_Toan(argv[2]);
			int para = Xac_dinh_parameter(argv[4]);

			cout << "Input order: Randomize\n ";
			cout << "----------------------------" << endl;
			generateData(arr, n, 0);
			Print_file_Input("input_1.txt", arr, n);
			Cout_consolve(al, para, arr, n);


			cout << "\nInput order: Nearly Sorted\n ";
			cout << "----------------------------" << endl;
			generateData(arr, n, 3);
			Print_file_Input("input_2.txt", arr, n);
			Cout_consolve(al, para, arr, n);


			cout << "\nInput order: Sorted\n ";
			cout << "----------------------------" << endl;
			generateData(arr, n, 1);
			Print_file_Input("input_3.txt", arr, n);
			Cout_consolve(al, para, arr, n);


			cout << "\nInput order: Reversed\n ";
			cout << "----------------------------" << endl;
			generateData(arr, n, 2);
			Print_file_Input("input_4.txt", arr, n);
			Cout_consolve(al, para, arr, n);


			delete[] arr;
		}

	}

	else if (strcmp(argv[1], "-c") == 0) {
		cout << "COMPARE MODE" << endl;

		if (argc == 5) { //command 4

			// fix algorithm name from cmd for printing
			string al1(argv[2]);
			al1[0] = toupper(al1[0]);
			string al2(argv[3]);
			al2[0] = toupper(al2[0]);
			for (int i = 1; i < al1.size(); i++)
				if (al1[i] == '-') {
					al1[i] = ' ';
					al1[i + 1] = toupper(al1[i + 1]);
				}
			for (int i = 1; i < al2.size(); i++)
				if (al2[i] == '-') {
					al2[i] = ' ';
					al2[i + 1] = toupper(al2[i + 1]);
				}





			// read input file
			ifstream f;
			f.open(argv[4], ios::binary);
			if (!f)
			{
				cout << "Cannot open file to read data" << endl;
				return 0;
			}
			// input size
			string line;
			getline(f, line);
			int inputSize = stoi(line);

			if (inputSize > 1000000) {
				cout << "Input size must be less than or equal to 1,000,000!" << endl;
				return 0;
			}

			// data
			int* a = new int[inputSize];
			string element;
			int i = 0;
			while (getline(f, line))
			{
				stringstream ss(line);
				i = 0;
				while (getline(ss, element, ' ') && i < inputSize)
				{
					a[i] = stoi(element);
					i++;
				}
			}
			f.close();

			// cout input algorithm + input file + input size
			cout << "Algorithm: " << al1 << " | " << al2 << endl;
			cout << "Input file: " << string(argv[4]) << endl;
			cout << "Input size: " << inputSize << endl;

			cout << "---------------------------" << endl;

			int* aTemp = new int[inputSize];
			cout << "Running time: ";
			for (int i = 0; i < 2; i++)
			{

				for (int j = 0; j < inputSize; j++)
				{
					aTemp[j] = a[j];
				}
				if (strcmp(argv[i + 2], "selection-sort") == 0)
					cout << fixed << setprecision(4) << timeSelectionSort(aTemp, inputSize) << " ms";
				else if (strcmp(argv[i + 2], "insertion-sort") == 0)
					cout << fixed << setprecision(4) << insertionSortRT(aTemp, inputSize) << " ms";
				else if (strcmp(argv[i + 2], "bubble-sort") == 0)
					cout << fixed << setprecision(4) << BubbleSort(aTemp, inputSize) << " ms";
				else if (strcmp(argv[i + 2], "heap-sort") == 0)
					cout << fixed << setprecision(4) << Heap_sort(aTemp, inputSize) << " ms";
				else if (strcmp(argv[i + 2], "merge-sort") == 0)
					cout << fixed << setprecision(4) << mergeSortRT(aTemp, inputSize) << " ms";
				else if (strcmp(argv[i + 2], "quick-sort") == 0)
					cout << fixed << setprecision(4) << timeQuickSort(aTemp, 0, inputSize - 1) << " ms";
				else if (strcmp(argv[i + 2], "radix-sort") == 0)
					cout << fixed << setprecision(4) << getTime(aTemp, inputSize) << " ms";
				else if (strcmp(argv[i + 2], "counting-sort") == 0)
					cout << fixed << setprecision(4) << countingSortRT(aTemp, inputSize) << " ms";
				else
					cout << "Algorithm doesn't exist";

				if (i == 0)
					cout << " | ";
			}
			cout << endl;
			cout << "Comparisons: ";
			for (int i = 0; i < 2; i++)
			{

				for (int j = 0; j < inputSize; j++)
				{
					aTemp[j] = a[j];
				}
				if (strcmp(argv[i + 2], "selection-sort") == 0)
					cout << countSelectionSort(aTemp, inputSize);
				else if (strcmp(argv[i + 2], "insertion-sort") == 0)
					cout << insertionSortCMP(aTemp, inputSize);
				else if (strcmp(argv[i + 2], "bubble-sort") == 0)
					cout << BubbleSort_Comparison(aTemp, inputSize);
				else if (strcmp(argv[i + 2], "heap-sort") == 0)
					cout << Heap_sort_Comparison(aTemp, inputSize);
				else if (strcmp(argv[i + 2], "merge-sort") == 0)
					cout << mergeSortCMP(aTemp, 0, inputSize - 1);
				else if (strcmp(argv[i + 2], "quick-sort") == 0)
					cout << countQuickSort(aTemp, 0, inputSize - 1);
				else if (strcmp(argv[i + 2], "radix-sort") == 0)
					cout << Counting(aTemp, inputSize);
				else if (strcmp(argv[i + 2], "counting-sort") == 0)
					cout << countingSortCMP(aTemp, inputSize);
				else
					cout << "Algorithm doesn't exist";

				if (i == 0)
					cout << " | ";
			}

			delete[] aTemp;
			delete[] a;

		}

		else if (argc == 6) { //command 5

			int inputSize = atoi(argv[4]);

			if (inputSize > 1000000) {
				cout << "Input size must be less than or equal to 1,000,000!" << endl;
				return 0;
			}

			string al1(argv[2]);
			al1[0] = toupper(al1[0]);
			string al2(argv[3]);
			al2[0] = toupper(al2[0]);

			for (int i = 1; i < al1.size(); i++)
				if (al1[i] == '-') {
					al1[i] = ' ';
					al1[i + 1] = toupper(al1[i + 1]);
				}

			for (int i = 1; i < al2.size(); i++)
				if (al2[i] == '-') {
					al2[i] = ' ';
					al2[i + 1] = toupper(al2[i + 1]);
				}

			int dataOrder = -1;

			if (strcmp(argv[5], "-rand") == 0) dataOrder = 0;
			else if (strcmp(argv[5], "-sorted") == 0) dataOrder = 1;
			else if (strcmp(argv[5], "-rev") == 0) dataOrder = 2;
			else if (strcmp(argv[5], "-nsorted") == 0) dataOrder = 3;
			else {
				cout << "Wrong data order!" << endl;
				return 0;
			}

			cout << "Algorithm: " << al1 << " | " << al2 << endl;
			cout << "Input size: " << inputSize << endl;
			cout << "Input order: ";

			if (dataOrder == 0) cout << "Randomized data" << endl;
			else if (dataOrder == 1) cout << "Sorted data" << endl;
			else if (dataOrder == 2) cout << "Reverse sorted data" << endl;
			else if (dataOrder == 3) cout << "Nearly sorted data" << endl;

			int* a = new int[inputSize];
			generateData(a, inputSize, dataOrder);

			ofstream f("input.txt");
			if (!f) {
				cout << "Cannot open file for writing data." << endl;
				return 0;
			}
			f << inputSize << endl;
			for (int i = 0; i < inputSize; i++)
			{
				f << a[i] << " ";
			}
			f.close();


			cout << "---------------------------" << endl;

			int* aTemp = new int[inputSize];

			cout << "Running time: ";
			for (int i = 0; i < 2; i++)
			{

				for (int j = 0; j < inputSize; j++)
				{
					aTemp[j] = a[j];
				}
				if (strcmp(argv[i + 2], "selection-sort") == 0)
					cout << fixed << setprecision(4) << timeSelectionSort(aTemp, inputSize) << " ms";
				else if (strcmp(argv[i + 2], "insertion-sort") == 0)
					cout << fixed << setprecision(4) << insertionSortRT(aTemp, inputSize) << " ms";
				else if (strcmp(argv[i + 2], "bubble-sort") == 0)
					cout << fixed << setprecision(4) << BubbleSort(aTemp, inputSize) << " ms";
				else if (strcmp(argv[i + 2], "heap-sort") == 0)
					cout << fixed << setprecision(4) << Heap_sort(aTemp, inputSize) << " ms";
				else if (strcmp(argv[i + 2], "merge-sort") == 0)
					cout << fixed << setprecision(4) << mergeSortRT(aTemp, inputSize) << " ms";
				else if (strcmp(argv[i + 2], "quick-sort") == 0)
					cout << fixed << setprecision(4) << timeQuickSort(aTemp, 0, inputSize - 1) << " ms";
				else if (strcmp(argv[i + 2], "radix-sort") == 0)
					cout << fixed << setprecision(4) << getTime(aTemp, inputSize) << " ms";
				else if (strcmp(argv[i + 2], "counting-sort") == 0)
					cout << fixed << setprecision(4) << countingSortRT(aTemp, inputSize) << " ms";
				else
					cout << "Algorithm doesn't exist";

				if (i == 0)
					cout << " | ";
			}

			cout << endl;

			cout << "Comparisons: ";
			for (int i = 0; i < 2; i++)
			{

				for (int j = 0; j < inputSize; j++)
				{
					aTemp[j] = a[j];
				}
				if (strcmp(argv[i + 2], "selection-sort") == 0)
					cout << countSelectionSort(aTemp, inputSize);
				else if (strcmp(argv[i + 2], "insertion-sort") == 0)
					cout << insertionSortCMP(aTemp, inputSize);
				else if (strcmp(argv[i + 2], "bubble-sort") == 0)
					cout << BubbleSort_Comparison(aTemp, inputSize);
				else if (strcmp(argv[i + 2], "heap-sort") == 0)
					cout << Heap_sort_Comparison(aTemp, inputSize);
				else if (strcmp(argv[i + 2], "merge-sort") == 0)
					cout << mergeSortCMP(aTemp, 0, inputSize - 1);
				else if (strcmp(argv[i + 2], "quick-sort") == 0)
					cout << countQuickSort(aTemp, 0, inputSize - 1);
				else if (strcmp(argv[i + 2], "radix-sort") == 0)
					cout << Counting(aTemp, inputSize);
				else if (strcmp(argv[i + 2], "counting-sort") == 0)
					cout << countingSortCMP(aTemp, inputSize);
				else
					cout << "Algorithm doesn't exist";

				if (i == 0)
					cout << " | ";
			}



			delete[] aTemp;
			delete[] a;
		}
	}

	return 0;
}
