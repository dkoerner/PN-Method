#include <iostream>
#include <string>

#include <cuda.h>
#include <builtin_types.h>
//#include <drvapi_error_string.h>



// Forward declare the function in the .cu file
void vectorAddition(const float* a, const float* b, float* c, int n);


void printArray(const float* a, const unsigned int n)
{
	std::string s = "(";
	//unsigned int ii;
	//for (ii = 0; ii < n - 1; ++ii)
	//	s += ing::number(a[ii])).append(", ");
	//s.append(QString::number(a[ii])).append(")");

	std::cout << s << std::endl;
}

int main()
{
	std::cout << "Hello World!" << std::endl;





	int deviceCount = 0;
	int cudaDevice = 0;
	char cudaDeviceName [100];

	unsigned int N = 50;
	float *a, *b, *c;

	cuInit(0);
	cuDeviceGetCount(&deviceCount);
	cuDeviceGet(&cudaDevice, 0);
	cuDeviceGetName(cudaDeviceName, 100, cudaDevice);
	std::cout << "Number of devices: " << deviceCount << std::endl;
	std::cout << "Device name:" << cudaDeviceName << std::endl;

	a = new float [N];    b = new float [N];    c = new float [N];
	for (unsigned int ii = 0; ii < N; ++ii)
	{
		a[ii] = 1;
		b[ii] = 2;
		//a[ii] = qrand();
		//b[ii] = qrand();
	}

	// This is the function call in which the kernel is called
	vectorAddition(a, b, c, N);


	std::cout << "input a:"; printArray(a, N);
	std::cout << "input b:"; printArray(b, N);
	std::cout << "output c:"; printArray(c, N);

	if (a) delete a;
	if (b) delete b;
	if (c) delete c;

	return 0;
}

