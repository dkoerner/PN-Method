#include <PNSystem.h>

void set_system_row(PNSystem::Voxel& sys,
					PNSystem::Fields& fields)
{
	V2i vi = sys.getVoxel();
	V2d vd = sys.getVoxel().cast<double>();

	// term 0 ----------------

	// term 1 ----------------
	sys.A( 2, vi + V2i(-2,0), 1 ) += 4.16666666667;
	sys.A( 2, vi + V2i(-2,1), 1 ) += 4.16666666667;
	sys.A( 2, vi + V2i(-1,0), 1 ) += 4.16666666667;
	sys.A( 2, vi + V2i(-1,1), 1 ) += 4.16666666667;
	sys.A( 2, vi + V2i(-1,0), 1 ) += -4.16666666667;
	sys.A( 2, vi + V2i(-1,1), 1 ) += -4.16666666667;
	sys.A( 2, vi + V2i(0,0), 1 ) += -4.16666666667;
	sys.A( 2, vi + V2i(0,1), 1 ) += -4.16666666667;
	sys.A( 2, vi + V2i(-1,0), 1 ) += -4.16666666667;
	sys.A( 2, vi + V2i(-1,1), 1 ) += -4.16666666667;
	sys.A( 2, vi + V2i(0,0), 1 ) += -4.16666666667;
	sys.A( 2, vi + V2i(0,1), 1 ) += -4.16666666667;
	sys.A( 2, vi + V2i(0,0), 1 ) += 4.16666666667;
	sys.A( 2, vi + V2i(0,1), 1 ) += 4.16666666667;
	sys.A( 2, vi + V2i(1,0), 1 ) += 4.16666666667;
	sys.A( 2, vi + V2i(1,1), 1 ) += 4.16666666667;

	// term 2 ----------------

	// term 3 ----------------
	sys.A( 2, vi + V2i(-1,0), 2 ) += -16.6666666667;
	sys.A( 2, vi + V2i(0,0), 2 ) += 16.6666666667;
	sys.A( 2, vi + V2i(0,0), 2 ) += 16.6666666667;
	sys.A( 2, vi + V2i(1,0), 2 ) += -16.6666666667;

	// term 4 ----------------

	// term 5 ----------------
	sys.A( 2, vi + V2i(-1,0), 1 ) += -std::complex<double>(0.0, 16.666666666666664);
	sys.A( 2, vi + V2i(-1,1), 1 ) += std::complex<double>(0.0, 16.666666666666664);
	sys.A( 2, vi + V2i(0,0), 1 ) += std::complex<double>(0.0, 16.666666666666664);
	sys.A( 2, vi + V2i(0,1), 1 ) += -std::complex<double>(0.0, 16.666666666666664);

	// term 6 ----------------

	// term 7 ----------------
	sys.A( 2, vi + V2i(-1,-1), 2 ) += -std::complex<double>(0.0, 4.166666666666666);
	sys.A( 2, vi + V2i(-1,0), 2 ) += -std::complex<double>(0.0, 4.166666666666666);
	sys.A( 2, vi + V2i(0,-1), 2 ) += -std::complex<double>(0.0, 4.166666666666666);
	sys.A( 2, vi + V2i(0,0), 2 ) += -std::complex<double>(0.0, 4.166666666666666);
	sys.A( 2, vi + V2i(-1,0), 2 ) += std::complex<double>(0.0, 4.166666666666666);
	sys.A( 2, vi + V2i(-1,1), 2 ) += std::complex<double>(0.0, 4.166666666666666);
	sys.A( 2, vi + V2i(0,0), 2 ) += std::complex<double>(0.0, 4.166666666666666);
	sys.A( 2, vi + V2i(0,1), 2 ) += std::complex<double>(0.0, 4.166666666666666);
	sys.A( 2, vi + V2i(0,-1), 2 ) += std::complex<double>(0.0, 4.166666666666666);
	sys.A( 2, vi + V2i(0,0), 2 ) += std::complex<double>(0.0, 4.166666666666666);
	sys.A( 2, vi + V2i(1,-1), 2 ) += std::complex<double>(0.0, 4.166666666666666);
	sys.A( 2, vi + V2i(1,0), 2 ) += std::complex<double>(0.0, 4.166666666666666);
	sys.A( 2, vi + V2i(0,0), 2 ) += -std::complex<double>(0.0, 4.166666666666666);
	sys.A( 2, vi + V2i(0,1), 2 ) += -std::complex<double>(0.0, 4.166666666666666);
	sys.A( 2, vi + V2i(1,0), 2 ) += -std::complex<double>(0.0, 4.166666666666666);
	sys.A( 2, vi + V2i(1,1), 2 ) += -std::complex<double>(0.0, 4.166666666666666);

	// term 8 ----------------

	// term 9 ----------------

	// term 10 ----------------
	sys.A( 2, vi + V2i(-2,0), 1 ) += 0.833333333333;
	sys.A( 2, vi + V2i(-2,1), 1 ) += 0.833333333333;
	sys.A( 2, vi + V2i(-1,0), 1 ) += 0.833333333333;
	sys.A( 2, vi + V2i(-1,1), 1 ) += 0.833333333333;
	sys.A( 2, vi + V2i(-1,0), 1 ) += -0.833333333333;
	sys.A( 2, vi + V2i(-1,1), 1 ) += -0.833333333333;
	sys.A( 2, vi + V2i(0,0), 1 ) += -0.833333333333;
	sys.A( 2, vi + V2i(0,1), 1 ) += -0.833333333333;
	sys.A( 2, vi + V2i(-1,0), 1 ) += -0.833333333333;
	sys.A( 2, vi + V2i(-1,1), 1 ) += -0.833333333333;
	sys.A( 2, vi + V2i(0,0), 1 ) += -0.833333333333;
	sys.A( 2, vi + V2i(0,1), 1 ) += -0.833333333333;
	sys.A( 2, vi + V2i(0,0), 1 ) += 0.833333333333;
	sys.A( 2, vi + V2i(0,1), 1 ) += 0.833333333333;
	sys.A( 2, vi + V2i(1,0), 1 ) += 0.833333333333;
	sys.A( 2, vi + V2i(1,1), 1 ) += 0.833333333333;

	// term 11 ----------------

	// term 12 ----------------
	sys.A( 0, vi + V2i(-1,0), 0 ) += -16.6666666667;
	sys.A( 0, vi + V2i(0,0), 0 ) += 16.6666666667;
	sys.A( 0, vi + V2i(0,0), 0 ) += 16.6666666667;
	sys.A( 0, vi + V2i(1,0), 0 ) += -16.6666666667;
	sys.A( 1, vi + V2i(-1,0), 1 ) += -20.0;
	sys.A( 1, vi + V2i(0,0), 1 ) += 20.0;
	sys.A( 1, vi + V2i(0,0), 1 ) += 20.0;
	sys.A( 1, vi + V2i(1,0), 1 ) += -20.0;
	sys.A( 2, vi + V2i(-1,0), 2 ) += -3.33333333333;
	sys.A( 2, vi + V2i(0,0), 2 ) += 3.33333333333;
	sys.A( 2, vi + V2i(0,0), 2 ) += 3.33333333333;
	sys.A( 2, vi + V2i(1,0), 2 ) += -3.33333333333;

	// term 13 ----------------

	// term 14 ----------------
	sys.A( 2, vi + V2i(-1,0), 1 ) += -std::complex<double>(0.0, 3.333333333333333);
	sys.A( 2, vi + V2i(-1,1), 1 ) += std::complex<double>(0.0, 3.333333333333333);
	sys.A( 2, vi + V2i(0,0), 1 ) += std::complex<double>(0.0, 3.333333333333333);
	sys.A( 2, vi + V2i(0,1), 1 ) += -std::complex<double>(0.0, 3.333333333333333);

	// term 15 ----------------

	// term 16 ----------------
	sys.A( 0, vi + V2i(-1,-1), 0 ) += -std::complex<double>(0.0, 4.166666666666666);
	sys.A( 0, vi + V2i(-1,0), 0 ) += -std::complex<double>(0.0, 4.166666666666666);
	sys.A( 0, vi + V2i(0,-1), 0 ) += -std::complex<double>(0.0, 4.166666666666666);
	sys.A( 0, vi + V2i(0,0), 0 ) += -std::complex<double>(0.0, 4.166666666666666);
	sys.A( 0, vi + V2i(-1,0), 0 ) += std::complex<double>(0.0, 4.166666666666666);
	sys.A( 0, vi + V2i(-1,1), 0 ) += std::complex<double>(0.0, 4.166666666666666);
	sys.A( 0, vi + V2i(0,0), 0 ) += std::complex<double>(0.0, 4.166666666666666);
	sys.A( 0, vi + V2i(0,1), 0 ) += std::complex<double>(0.0, 4.166666666666666);
	sys.A( 0, vi + V2i(0,-1), 0 ) += std::complex<double>(0.0, 4.166666666666666);
	sys.A( 0, vi + V2i(0,0), 0 ) += std::complex<double>(0.0, 4.166666666666666);
	sys.A( 0, vi + V2i(1,-1), 0 ) += std::complex<double>(0.0, 4.166666666666666);
	sys.A( 0, vi + V2i(1,0), 0 ) += std::complex<double>(0.0, 4.166666666666666);
	sys.A( 0, vi + V2i(0,0), 0 ) += -std::complex<double>(0.0, 4.166666666666666);
	sys.A( 0, vi + V2i(0,1), 0 ) += -std::complex<double>(0.0, 4.166666666666666);
	sys.A( 0, vi + V2i(1,0), 0 ) += -std::complex<double>(0.0, 4.166666666666666);
	sys.A( 0, vi + V2i(1,1), 0 ) += -std::complex<double>(0.0, 4.166666666666666);
	sys.A( 1, vi + V2i(-1,-1), 1 ) += -std::complex<double>(0.0, 5.0);
	sys.A( 1, vi + V2i(-1,0), 1 ) += -std::complex<double>(0.0, 5.0);
	sys.A( 1, vi + V2i(0,-1), 1 ) += -std::complex<double>(0.0, 5.0);
	sys.A( 1, vi + V2i(0,0), 1 ) += -std::complex<double>(0.0, 5.0);
	sys.A( 1, vi + V2i(-1,0), 1 ) += std::complex<double>(0.0, 5.0);
	sys.A( 1, vi + V2i(-1,1), 1 ) += std::complex<double>(0.0, 5.0);
	sys.A( 1, vi + V2i(0,0), 1 ) += std::complex<double>(0.0, 5.0);
	sys.A( 1, vi + V2i(0,1), 1 ) += std::complex<double>(0.0, 5.0);
	sys.A( 1, vi + V2i(0,-1), 1 ) += std::complex<double>(0.0, 5.0);
	sys.A( 1, vi + V2i(0,0), 1 ) += std::complex<double>(0.0, 5.0);
	sys.A( 1, vi + V2i(1,-1), 1 ) += std::complex<double>(0.0, 5.0);
	sys.A( 1, vi + V2i(1,0), 1 ) += std::complex<double>(0.0, 5.0);
	sys.A( 1, vi + V2i(0,0), 1 ) += -std::complex<double>(0.0, 5.0);
	sys.A( 1, vi + V2i(0,1), 1 ) += -std::complex<double>(0.0, 5.0);
	sys.A( 1, vi + V2i(1,0), 1 ) += -std::complex<double>(0.0, 5.0);
	sys.A( 1, vi + V2i(1,1), 1 ) += -std::complex<double>(0.0, 5.0);
	sys.A( 2, vi + V2i(-1,-1), 2 ) += -std::complex<double>(0.0, 0.8333333333333333);
	sys.A( 2, vi + V2i(-1,0), 2 ) += -std::complex<double>(0.0, 0.8333333333333333);
	sys.A( 2, vi + V2i(0,-1), 2 ) += -std::complex<double>(0.0, 0.8333333333333333);
	sys.A( 2, vi + V2i(0,0), 2 ) += -std::complex<double>(0.0, 0.8333333333333333);
	sys.A( 2, vi + V2i(-1,0), 2 ) += std::complex<double>(0.0, 0.8333333333333333);
	sys.A( 2, vi + V2i(-1,1), 2 ) += std::complex<double>(0.0, 0.8333333333333333);
	sys.A( 2, vi + V2i(0,0), 2 ) += std::complex<double>(0.0, 0.8333333333333333);
	sys.A( 2, vi + V2i(0,1), 2 ) += std::complex<double>(0.0, 0.8333333333333333);
	sys.A( 2, vi + V2i(0,-1), 2 ) += std::complex<double>(0.0, 0.8333333333333333);
	sys.A( 2, vi + V2i(0,0), 2 ) += std::complex<double>(0.0, 0.8333333333333333);
	sys.A( 2, vi + V2i(1,-1), 2 ) += std::complex<double>(0.0, 0.8333333333333333);
	sys.A( 2, vi + V2i(1,0), 2 ) += std::complex<double>(0.0, 0.8333333333333333);
	sys.A( 2, vi + V2i(0,0), 2 ) += -std::complex<double>(0.0, 0.8333333333333333);
	sys.A( 2, vi + V2i(0,1), 2 ) += -std::complex<double>(0.0, 0.8333333333333333);
	sys.A( 2, vi + V2i(1,0), 2 ) += -std::complex<double>(0.0, 0.8333333333333333);
	sys.A( 2, vi + V2i(1,1), 2 ) += -std::complex<double>(0.0, 0.8333333333333333);

	// term 17 ----------------

	// term 18 ----------------

	// term 19 ----------------

	// term 20 ----------------

	// term 21 ----------------
	sys.A( 1, vi + V2i(-1,0), 1 ) += -16.6666666667;
	sys.A( 1, vi + V2i(0,0), 1 ) += 16.6666666667;
	sys.A( 1, vi + V2i(0,0), 1 ) += 16.6666666667;
	sys.A( 1, vi + V2i(1,0), 1 ) += -16.6666666667;

	// term 22 ----------------

	// term 23 ----------------
	sys.A( 1, vi + V2i(-1,-1), 2 ) += 4.16666666667;
	sys.A( 1, vi + V2i(-1,0), 2 ) += 4.16666666667;
	sys.A( 1, vi + V2i(0,-1), 2 ) += 4.16666666667;
	sys.A( 1, vi + V2i(0,0), 2 ) += 4.16666666667;
	sys.A( 1, vi + V2i(0,-1), 2 ) += -4.16666666667;
	sys.A( 1, vi + V2i(0,0), 2 ) += -4.16666666667;
	sys.A( 1, vi + V2i(1,-1), 2 ) += -4.16666666667;
	sys.A( 1, vi + V2i(1,0), 2 ) += -4.16666666667;
	sys.A( 1, vi + V2i(0,-1), 2 ) += -4.16666666667;
	sys.A( 1, vi + V2i(0,0), 2 ) += -4.16666666667;
	sys.A( 1, vi + V2i(1,-1), 2 ) += -4.16666666667;
	sys.A( 1, vi + V2i(1,0), 2 ) += -4.16666666667;
	sys.A( 1, vi + V2i(1,-1), 2 ) += 4.16666666667;
	sys.A( 1, vi + V2i(1,0), 2 ) += 4.16666666667;
	sys.A( 1, vi + V2i(2,-1), 2 ) += 4.16666666667;
	sys.A( 1, vi + V2i(2,0), 2 ) += 4.16666666667;

	// term 24 ----------------

	// term 25 ----------------
	sys.A( 1, vi + V2i(-1,-1), 1 ) += std::complex<double>(0.0, 4.166666666666666);
	sys.A( 1, vi + V2i(-1,0), 1 ) += std::complex<double>(0.0, 4.166666666666666);
	sys.A( 1, vi + V2i(0,-1), 1 ) += std::complex<double>(0.0, 4.166666666666666);
	sys.A( 1, vi + V2i(0,0), 1 ) += std::complex<double>(0.0, 4.166666666666666);
	sys.A( 1, vi + V2i(-1,0), 1 ) += -std::complex<double>(0.0, 4.166666666666666);
	sys.A( 1, vi + V2i(-1,1), 1 ) += -std::complex<double>(0.0, 4.166666666666666);
	sys.A( 1, vi + V2i(0,0), 1 ) += -std::complex<double>(0.0, 4.166666666666666);
	sys.A( 1, vi + V2i(0,1), 1 ) += -std::complex<double>(0.0, 4.166666666666666);
	sys.A( 1, vi + V2i(0,-1), 1 ) += -std::complex<double>(0.0, 4.166666666666666);
	sys.A( 1, vi + V2i(0,0), 1 ) += -std::complex<double>(0.0, 4.166666666666666);
	sys.A( 1, vi + V2i(1,-1), 1 ) += -std::complex<double>(0.0, 4.166666666666666);
	sys.A( 1, vi + V2i(1,0), 1 ) += -std::complex<double>(0.0, 4.166666666666666);
	sys.A( 1, vi + V2i(0,0), 1 ) += std::complex<double>(0.0, 4.166666666666666);
	sys.A( 1, vi + V2i(0,1), 1 ) += std::complex<double>(0.0, 4.166666666666666);
	sys.A( 1, vi + V2i(1,0), 1 ) += std::complex<double>(0.0, 4.166666666666666);
	sys.A( 1, vi + V2i(1,1), 1 ) += std::complex<double>(0.0, 4.166666666666666);

	// term 26 ----------------

	// term 27 ----------------
	sys.A( 1, vi + V2i(0,-1), 2 ) += std::complex<double>(0.0, 16.666666666666664);
	sys.A( 1, vi + V2i(0,0), 2 ) += -std::complex<double>(0.0, 16.666666666666664);
	sys.A( 1, vi + V2i(1,-1), 2 ) += -std::complex<double>(0.0, 16.666666666666664);
	sys.A( 1, vi + V2i(1,0), 2 ) += std::complex<double>(0.0, 16.666666666666664);

	// term 28 ----------------

	// term 29 ----------------

	// term 30 ----------------
	sys.A( 0, vi + V2i(-1,0), 0 ) += -16.6666666667;
	sys.A( 0, vi + V2i(0,0), 0 ) += 16.6666666667;
	sys.A( 0, vi + V2i(0,0), 0 ) += 16.6666666667;
	sys.A( 0, vi + V2i(1,0), 0 ) += -16.6666666667;
	sys.A( 1, vi + V2i(-1,0), 1 ) += -3.33333333333;
	sys.A( 1, vi + V2i(0,0), 1 ) += 3.33333333333;
	sys.A( 1, vi + V2i(0,0), 1 ) += 3.33333333333;
	sys.A( 1, vi + V2i(1,0), 1 ) += -3.33333333333;
	sys.A( 2, vi + V2i(-1,0), 2 ) += -20.0;
	sys.A( 2, vi + V2i(0,0), 2 ) += 20.0;
	sys.A( 2, vi + V2i(0,0), 2 ) += 20.0;
	sys.A( 2, vi + V2i(1,0), 2 ) += -20.0;

	// term 31 ----------------

	// term 32 ----------------
	sys.A( 1, vi + V2i(-1,-1), 2 ) += 0.833333333333;
	sys.A( 1, vi + V2i(-1,0), 2 ) += 0.833333333333;
	sys.A( 1, vi + V2i(0,-1), 2 ) += 0.833333333333;
	sys.A( 1, vi + V2i(0,0), 2 ) += 0.833333333333;
	sys.A( 1, vi + V2i(0,-1), 2 ) += -0.833333333333;
	sys.A( 1, vi + V2i(0,0), 2 ) += -0.833333333333;
	sys.A( 1, vi + V2i(1,-1), 2 ) += -0.833333333333;
	sys.A( 1, vi + V2i(1,0), 2 ) += -0.833333333333;
	sys.A( 1, vi + V2i(0,-1), 2 ) += -0.833333333333;
	sys.A( 1, vi + V2i(0,0), 2 ) += -0.833333333333;
	sys.A( 1, vi + V2i(1,-1), 2 ) += -0.833333333333;
	sys.A( 1, vi + V2i(1,0), 2 ) += -0.833333333333;
	sys.A( 1, vi + V2i(1,-1), 2 ) += 0.833333333333;
	sys.A( 1, vi + V2i(1,0), 2 ) += 0.833333333333;
	sys.A( 1, vi + V2i(2,-1), 2 ) += 0.833333333333;
	sys.A( 1, vi + V2i(2,0), 2 ) += 0.833333333333;

	// term 33 ----------------

	// term 34 ----------------
	sys.A( 0, vi + V2i(-1,-1), 0 ) += std::complex<double>(0.0, 4.166666666666666);
	sys.A( 0, vi + V2i(-1,0), 0 ) += std::complex<double>(0.0, 4.166666666666666);
	sys.A( 0, vi + V2i(0,-1), 0 ) += std::complex<double>(0.0, 4.166666666666666);
	sys.A( 0, vi + V2i(0,0), 0 ) += std::complex<double>(0.0, 4.166666666666666);
	sys.A( 0, vi + V2i(-1,0), 0 ) += -std::complex<double>(0.0, 4.166666666666666);
	sys.A( 0, vi + V2i(-1,1), 0 ) += -std::complex<double>(0.0, 4.166666666666666);
	sys.A( 0, vi + V2i(0,0), 0 ) += -std::complex<double>(0.0, 4.166666666666666);
	sys.A( 0, vi + V2i(0,1), 0 ) += -std::complex<double>(0.0, 4.166666666666666);
	sys.A( 0, vi + V2i(0,-1), 0 ) += -std::complex<double>(0.0, 4.166666666666666);
	sys.A( 0, vi + V2i(0,0), 0 ) += -std::complex<double>(0.0, 4.166666666666666);
	sys.A( 0, vi + V2i(1,-1), 0 ) += -std::complex<double>(0.0, 4.166666666666666);
	sys.A( 0, vi + V2i(1,0), 0 ) += -std::complex<double>(0.0, 4.166666666666666);
	sys.A( 0, vi + V2i(0,0), 0 ) += std::complex<double>(0.0, 4.166666666666666);
	sys.A( 0, vi + V2i(0,1), 0 ) += std::complex<double>(0.0, 4.166666666666666);
	sys.A( 0, vi + V2i(1,0), 0 ) += std::complex<double>(0.0, 4.166666666666666);
	sys.A( 0, vi + V2i(1,1), 0 ) += std::complex<double>(0.0, 4.166666666666666);
	sys.A( 1, vi + V2i(-1,-1), 1 ) += std::complex<double>(0.0, 0.8333333333333333);
	sys.A( 1, vi + V2i(-1,0), 1 ) += std::complex<double>(0.0, 0.8333333333333333);
	sys.A( 1, vi + V2i(0,-1), 1 ) += std::complex<double>(0.0, 0.8333333333333333);
	sys.A( 1, vi + V2i(0,0), 1 ) += std::complex<double>(0.0, 0.8333333333333333);
	sys.A( 1, vi + V2i(-1,0), 1 ) += -std::complex<double>(0.0, 0.8333333333333333);
	sys.A( 1, vi + V2i(-1,1), 1 ) += -std::complex<double>(0.0, 0.8333333333333333);
	sys.A( 1, vi + V2i(0,0), 1 ) += -std::complex<double>(0.0, 0.8333333333333333);
	sys.A( 1, vi + V2i(0,1), 1 ) += -std::complex<double>(0.0, 0.8333333333333333);
	sys.A( 1, vi + V2i(0,-1), 1 ) += -std::complex<double>(0.0, 0.8333333333333333);
	sys.A( 1, vi + V2i(0,0), 1 ) += -std::complex<double>(0.0, 0.8333333333333333);
	sys.A( 1, vi + V2i(1,-1), 1 ) += -std::complex<double>(0.0, 0.8333333333333333);
	sys.A( 1, vi + V2i(1,0), 1 ) += -std::complex<double>(0.0, 0.8333333333333333);
	sys.A( 1, vi + V2i(0,0), 1 ) += std::complex<double>(0.0, 0.8333333333333333);
	sys.A( 1, vi + V2i(0,1), 1 ) += std::complex<double>(0.0, 0.8333333333333333);
	sys.A( 1, vi + V2i(1,0), 1 ) += std::complex<double>(0.0, 0.8333333333333333);
	sys.A( 1, vi + V2i(1,1), 1 ) += std::complex<double>(0.0, 0.8333333333333333);
	sys.A( 2, vi + V2i(-1,-1), 2 ) += std::complex<double>(0.0, 5.0);
	sys.A( 2, vi + V2i(-1,0), 2 ) += std::complex<double>(0.0, 5.0);
	sys.A( 2, vi + V2i(0,-1), 2 ) += std::complex<double>(0.0, 5.0);
	sys.A( 2, vi + V2i(0,0), 2 ) += std::complex<double>(0.0, 5.0);
	sys.A( 2, vi + V2i(-1,0), 2 ) += -std::complex<double>(0.0, 5.0);
	sys.A( 2, vi + V2i(-1,1), 2 ) += -std::complex<double>(0.0, 5.0);
	sys.A( 2, vi + V2i(0,0), 2 ) += -std::complex<double>(0.0, 5.0);
	sys.A( 2, vi + V2i(0,1), 2 ) += -std::complex<double>(0.0, 5.0);
	sys.A( 2, vi + V2i(0,-1), 2 ) += -std::complex<double>(0.0, 5.0);
	sys.A( 2, vi + V2i(0,0), 2 ) += -std::complex<double>(0.0, 5.0);
	sys.A( 2, vi + V2i(1,-1), 2 ) += -std::complex<double>(0.0, 5.0);
	sys.A( 2, vi + V2i(1,0), 2 ) += -std::complex<double>(0.0, 5.0);
	sys.A( 2, vi + V2i(0,0), 2 ) += std::complex<double>(0.0, 5.0);
	sys.A( 2, vi + V2i(0,1), 2 ) += std::complex<double>(0.0, 5.0);
	sys.A( 2, vi + V2i(1,0), 2 ) += std::complex<double>(0.0, 5.0);
	sys.A( 2, vi + V2i(1,1), 2 ) += std::complex<double>(0.0, 5.0);

	// term 35 ----------------

	// term 36 ----------------
	sys.A( 1, vi + V2i(0,-1), 2 ) += std::complex<double>(0.0, 3.333333333333333);
	sys.A( 1, vi + V2i(0,0), 2 ) += -std::complex<double>(0.0, 3.333333333333333);
	sys.A( 1, vi + V2i(1,-1), 2 ) += -std::complex<double>(0.0, 3.333333333333333);
	sys.A( 1, vi + V2i(1,0), 2 ) += std::complex<double>(0.0, 3.333333333333333);

	// term 37 ----------------

	// term 38 ----------------

	// term 39 ----------------

	// term 40 ----------------

	// term 41 ----------------
	sys.A( 2, vi + V2i(-1,0), 1 ) += -std::complex<double>(0.0, 16.666666666666664);
	sys.A( 2, vi + V2i(0,0), 1 ) += std::complex<double>(0.0, 16.666666666666664);
	sys.A( 2, vi + V2i(-1,1), 1 ) += std::complex<double>(0.0, 16.666666666666664);
	sys.A( 2, vi + V2i(0,1), 1 ) += -std::complex<double>(0.0, 16.666666666666664);

	// term 42 ----------------

	// term 43 ----------------
	sys.A( 2, vi + V2i(-1,-1), 2 ) += std::complex<double>(0.0, 4.166666666666666);
	sys.A( 2, vi + V2i(-1,0), 2 ) += std::complex<double>(0.0, 4.166666666666666);
	sys.A( 2, vi + V2i(0,-1), 2 ) += std::complex<double>(0.0, 4.166666666666666);
	sys.A( 2, vi + V2i(0,0), 2 ) += std::complex<double>(0.0, 4.166666666666666);
	sys.A( 2, vi + V2i(0,-1), 2 ) += -std::complex<double>(0.0, 4.166666666666666);
	sys.A( 2, vi + V2i(0,0), 2 ) += -std::complex<double>(0.0, 4.166666666666666);
	sys.A( 2, vi + V2i(1,-1), 2 ) += -std::complex<double>(0.0, 4.166666666666666);
	sys.A( 2, vi + V2i(1,0), 2 ) += -std::complex<double>(0.0, 4.166666666666666);
	sys.A( 2, vi + V2i(-1,0), 2 ) += -std::complex<double>(0.0, 4.166666666666666);
	sys.A( 2, vi + V2i(-1,1), 2 ) += -std::complex<double>(0.0, 4.166666666666666);
	sys.A( 2, vi + V2i(0,0), 2 ) += -std::complex<double>(0.0, 4.166666666666666);
	sys.A( 2, vi + V2i(0,1), 2 ) += -std::complex<double>(0.0, 4.166666666666666);
	sys.A( 2, vi + V2i(0,0), 2 ) += std::complex<double>(0.0, 4.166666666666666);
	sys.A( 2, vi + V2i(0,1), 2 ) += std::complex<double>(0.0, 4.166666666666666);
	sys.A( 2, vi + V2i(1,0), 2 ) += std::complex<double>(0.0, 4.166666666666666);
	sys.A( 2, vi + V2i(1,1), 2 ) += std::complex<double>(0.0, 4.166666666666666);

	// term 44 ----------------

	// term 45 ----------------
	sys.A( 2, vi + V2i(-1,-1), 1 ) += -4.16666666667;
	sys.A( 2, vi + V2i(-1,0), 1 ) += -4.16666666667;
	sys.A( 2, vi + V2i(0,-1), 1 ) += -4.16666666667;
	sys.A( 2, vi + V2i(0,0), 1 ) += -4.16666666667;
	sys.A( 2, vi + V2i(-1,0), 1 ) += 4.16666666667;
	sys.A( 2, vi + V2i(-1,1), 1 ) += 4.16666666667;
	sys.A( 2, vi + V2i(0,0), 1 ) += 4.16666666667;
	sys.A( 2, vi + V2i(0,1), 1 ) += 4.16666666667;
	sys.A( 2, vi + V2i(-1,0), 1 ) += 4.16666666667;
	sys.A( 2, vi + V2i(-1,1), 1 ) += 4.16666666667;
	sys.A( 2, vi + V2i(0,0), 1 ) += 4.16666666667;
	sys.A( 2, vi + V2i(0,1), 1 ) += 4.16666666667;
	sys.A( 2, vi + V2i(-1,1), 1 ) += -4.16666666667;
	sys.A( 2, vi + V2i(-1,2), 1 ) += -4.16666666667;
	sys.A( 2, vi + V2i(0,1), 1 ) += -4.16666666667;
	sys.A( 2, vi + V2i(0,2), 1 ) += -4.16666666667;

	// term 46 ----------------

	// term 47 ----------------
	sys.A( 2, vi + V2i(0,-1), 2 ) += -16.6666666667;
	sys.A( 2, vi + V2i(0,0), 2 ) += 16.6666666667;
	sys.A( 2, vi + V2i(0,0), 2 ) += 16.6666666667;
	sys.A( 2, vi + V2i(0,1), 2 ) += -16.6666666667;

	// term 48 ----------------

	// term 49 ----------------

	// term 50 ----------------
	sys.A( 2, vi + V2i(-1,0), 1 ) += -std::complex<double>(0.0, 3.333333333333333);
	sys.A( 2, vi + V2i(0,0), 1 ) += std::complex<double>(0.0, 3.333333333333333);
	sys.A( 2, vi + V2i(-1,1), 1 ) += std::complex<double>(0.0, 3.333333333333333);
	sys.A( 2, vi + V2i(0,1), 1 ) += -std::complex<double>(0.0, 3.333333333333333);

	// term 51 ----------------

	// term 52 ----------------
	sys.A( 0, vi + V2i(-1,-1), 0 ) += std::complex<double>(0.0, 4.166666666666666);
	sys.A( 0, vi + V2i(-1,0), 0 ) += std::complex<double>(0.0, 4.166666666666666);
	sys.A( 0, vi + V2i(0,-1), 0 ) += std::complex<double>(0.0, 4.166666666666666);
	sys.A( 0, vi + V2i(0,0), 0 ) += std::complex<double>(0.0, 4.166666666666666);
	sys.A( 0, vi + V2i(0,-1), 0 ) += -std::complex<double>(0.0, 4.166666666666666);
	sys.A( 0, vi + V2i(0,0), 0 ) += -std::complex<double>(0.0, 4.166666666666666);
	sys.A( 0, vi + V2i(1,-1), 0 ) += -std::complex<double>(0.0, 4.166666666666666);
	sys.A( 0, vi + V2i(1,0), 0 ) += -std::complex<double>(0.0, 4.166666666666666);
	sys.A( 0, vi + V2i(-1,0), 0 ) += -std::complex<double>(0.0, 4.166666666666666);
	sys.A( 0, vi + V2i(-1,1), 0 ) += -std::complex<double>(0.0, 4.166666666666666);
	sys.A( 0, vi + V2i(0,0), 0 ) += -std::complex<double>(0.0, 4.166666666666666);
	sys.A( 0, vi + V2i(0,1), 0 ) += -std::complex<double>(0.0, 4.166666666666666);
	sys.A( 0, vi + V2i(0,0), 0 ) += std::complex<double>(0.0, 4.166666666666666);
	sys.A( 0, vi + V2i(0,1), 0 ) += std::complex<double>(0.0, 4.166666666666666);
	sys.A( 0, vi + V2i(1,0), 0 ) += std::complex<double>(0.0, 4.166666666666666);
	sys.A( 0, vi + V2i(1,1), 0 ) += std::complex<double>(0.0, 4.166666666666666);
	sys.A( 1, vi + V2i(-1,-1), 1 ) += std::complex<double>(0.0, 5.0);
	sys.A( 1, vi + V2i(-1,0), 1 ) += std::complex<double>(0.0, 5.0);
	sys.A( 1, vi + V2i(0,-1), 1 ) += std::complex<double>(0.0, 5.0);
	sys.A( 1, vi + V2i(0,0), 1 ) += std::complex<double>(0.0, 5.0);
	sys.A( 1, vi + V2i(0,-1), 1 ) += -std::complex<double>(0.0, 5.0);
	sys.A( 1, vi + V2i(0,0), 1 ) += -std::complex<double>(0.0, 5.0);
	sys.A( 1, vi + V2i(1,-1), 1 ) += -std::complex<double>(0.0, 5.0);
	sys.A( 1, vi + V2i(1,0), 1 ) += -std::complex<double>(0.0, 5.0);
	sys.A( 1, vi + V2i(-1,0), 1 ) += -std::complex<double>(0.0, 5.0);
	sys.A( 1, vi + V2i(-1,1), 1 ) += -std::complex<double>(0.0, 5.0);
	sys.A( 1, vi + V2i(0,0), 1 ) += -std::complex<double>(0.0, 5.0);
	sys.A( 1, vi + V2i(0,1), 1 ) += -std::complex<double>(0.0, 5.0);
	sys.A( 1, vi + V2i(0,0), 1 ) += std::complex<double>(0.0, 5.0);
	sys.A( 1, vi + V2i(0,1), 1 ) += std::complex<double>(0.0, 5.0);
	sys.A( 1, vi + V2i(1,0), 1 ) += std::complex<double>(0.0, 5.0);
	sys.A( 1, vi + V2i(1,1), 1 ) += std::complex<double>(0.0, 5.0);
	sys.A( 2, vi + V2i(-1,-1), 2 ) += std::complex<double>(0.0, 0.8333333333333333);
	sys.A( 2, vi + V2i(-1,0), 2 ) += std::complex<double>(0.0, 0.8333333333333333);
	sys.A( 2, vi + V2i(0,-1), 2 ) += std::complex<double>(0.0, 0.8333333333333333);
	sys.A( 2, vi + V2i(0,0), 2 ) += std::complex<double>(0.0, 0.8333333333333333);
	sys.A( 2, vi + V2i(0,-1), 2 ) += -std::complex<double>(0.0, 0.8333333333333333);
	sys.A( 2, vi + V2i(0,0), 2 ) += -std::complex<double>(0.0, 0.8333333333333333);
	sys.A( 2, vi + V2i(1,-1), 2 ) += -std::complex<double>(0.0, 0.8333333333333333);
	sys.A( 2, vi + V2i(1,0), 2 ) += -std::complex<double>(0.0, 0.8333333333333333);
	sys.A( 2, vi + V2i(-1,0), 2 ) += -std::complex<double>(0.0, 0.8333333333333333);
	sys.A( 2, vi + V2i(-1,1), 2 ) += -std::complex<double>(0.0, 0.8333333333333333);
	sys.A( 2, vi + V2i(0,0), 2 ) += -std::complex<double>(0.0, 0.8333333333333333);
	sys.A( 2, vi + V2i(0,1), 2 ) += -std::complex<double>(0.0, 0.8333333333333333);
	sys.A( 2, vi + V2i(0,0), 2 ) += std::complex<double>(0.0, 0.8333333333333333);
	sys.A( 2, vi + V2i(0,1), 2 ) += std::complex<double>(0.0, 0.8333333333333333);
	sys.A( 2, vi + V2i(1,0), 2 ) += std::complex<double>(0.0, 0.8333333333333333);
	sys.A( 2, vi + V2i(1,1), 2 ) += std::complex<double>(0.0, 0.8333333333333333);

	// term 53 ----------------

	// term 54 ----------------
	sys.A( 2, vi + V2i(-1,-1), 1 ) += -0.833333333333;
	sys.A( 2, vi + V2i(-1,0), 1 ) += -0.833333333333;
	sys.A( 2, vi + V2i(0,-1), 1 ) += -0.833333333333;
	sys.A( 2, vi + V2i(0,0), 1 ) += -0.833333333333;
	sys.A( 2, vi + V2i(-1,0), 1 ) += 0.833333333333;
	sys.A( 2, vi + V2i(-1,1), 1 ) += 0.833333333333;
	sys.A( 2, vi + V2i(0,0), 1 ) += 0.833333333333;
	sys.A( 2, vi + V2i(0,1), 1 ) += 0.833333333333;
	sys.A( 2, vi + V2i(-1,0), 1 ) += 0.833333333333;
	sys.A( 2, vi + V2i(-1,1), 1 ) += 0.833333333333;
	sys.A( 2, vi + V2i(0,0), 1 ) += 0.833333333333;
	sys.A( 2, vi + V2i(0,1), 1 ) += 0.833333333333;
	sys.A( 2, vi + V2i(-1,1), 1 ) += -0.833333333333;
	sys.A( 2, vi + V2i(-1,2), 1 ) += -0.833333333333;
	sys.A( 2, vi + V2i(0,1), 1 ) += -0.833333333333;
	sys.A( 2, vi + V2i(0,2), 1 ) += -0.833333333333;

	// term 55 ----------------

	// term 56 ----------------
	sys.A( 0, vi + V2i(0,-1), 0 ) += -16.6666666667;
	sys.A( 0, vi + V2i(0,0), 0 ) += 16.6666666667;
	sys.A( 0, vi + V2i(0,0), 0 ) += 16.6666666667;
	sys.A( 0, vi + V2i(0,1), 0 ) += -16.6666666667;
	sys.A( 1, vi + V2i(0,-1), 1 ) += -20.0;
	sys.A( 1, vi + V2i(0,0), 1 ) += 20.0;
	sys.A( 1, vi + V2i(0,0), 1 ) += 20.0;
	sys.A( 1, vi + V2i(0,1), 1 ) += -20.0;
	sys.A( 2, vi + V2i(0,-1), 2 ) += -3.33333333333;
	sys.A( 2, vi + V2i(0,0), 2 ) += 3.33333333333;
	sys.A( 2, vi + V2i(0,0), 2 ) += 3.33333333333;
	sys.A( 2, vi + V2i(0,1), 2 ) += -3.33333333333;

	// term 57 ----------------

	// term 58 ----------------

	// term 59 ----------------

	// term 60 ----------------

	// term 61 ----------------
	sys.A( 1, vi + V2i(-1,-1), 1 ) += -std::complex<double>(0.0, 4.166666666666666);
	sys.A( 1, vi + V2i(-1,0), 1 ) += -std::complex<double>(0.0, 4.166666666666666);
	sys.A( 1, vi + V2i(0,-1), 1 ) += -std::complex<double>(0.0, 4.166666666666666);
	sys.A( 1, vi + V2i(0,0), 1 ) += -std::complex<double>(0.0, 4.166666666666666);
	sys.A( 1, vi + V2i(0,-1), 1 ) += std::complex<double>(0.0, 4.166666666666666);
	sys.A( 1, vi + V2i(0,0), 1 ) += std::complex<double>(0.0, 4.166666666666666);
	sys.A( 1, vi + V2i(1,-1), 1 ) += std::complex<double>(0.0, 4.166666666666666);
	sys.A( 1, vi + V2i(1,0), 1 ) += std::complex<double>(0.0, 4.166666666666666);
	sys.A( 1, vi + V2i(-1,0), 1 ) += std::complex<double>(0.0, 4.166666666666666);
	sys.A( 1, vi + V2i(-1,1), 1 ) += std::complex<double>(0.0, 4.166666666666666);
	sys.A( 1, vi + V2i(0,0), 1 ) += std::complex<double>(0.0, 4.166666666666666);
	sys.A( 1, vi + V2i(0,1), 1 ) += std::complex<double>(0.0, 4.166666666666666);
	sys.A( 1, vi + V2i(0,0), 1 ) += -std::complex<double>(0.0, 4.166666666666666);
	sys.A( 1, vi + V2i(0,1), 1 ) += -std::complex<double>(0.0, 4.166666666666666);
	sys.A( 1, vi + V2i(1,0), 1 ) += -std::complex<double>(0.0, 4.166666666666666);
	sys.A( 1, vi + V2i(1,1), 1 ) += -std::complex<double>(0.0, 4.166666666666666);

	// term 62 ----------------

	// term 63 ----------------
	sys.A( 1, vi + V2i(0,-1), 2 ) += std::complex<double>(0.0, 16.666666666666664);
	sys.A( 1, vi + V2i(1,-1), 2 ) += -std::complex<double>(0.0, 16.666666666666664);
	sys.A( 1, vi + V2i(0,0), 2 ) += -std::complex<double>(0.0, 16.666666666666664);
	sys.A( 1, vi + V2i(1,0), 2 ) += std::complex<double>(0.0, 16.666666666666664);

	// term 64 ----------------

	// term 65 ----------------
	sys.A( 1, vi + V2i(0,-1), 1 ) += -16.6666666667;
	sys.A( 1, vi + V2i(0,0), 1 ) += 16.6666666667;
	sys.A( 1, vi + V2i(0,0), 1 ) += 16.6666666667;
	sys.A( 1, vi + V2i(0,1), 1 ) += -16.6666666667;

	// term 66 ----------------

	// term 67 ----------------
	sys.A( 1, vi + V2i(0,-2), 2 ) += -4.16666666667;
	sys.A( 1, vi + V2i(0,-1), 2 ) += -4.16666666667;
	sys.A( 1, vi + V2i(1,-2), 2 ) += -4.16666666667;
	sys.A( 1, vi + V2i(1,-1), 2 ) += -4.16666666667;
	sys.A( 1, vi + V2i(0,-1), 2 ) += 4.16666666667;
	sys.A( 1, vi + V2i(0,0), 2 ) += 4.16666666667;
	sys.A( 1, vi + V2i(1,-1), 2 ) += 4.16666666667;
	sys.A( 1, vi + V2i(1,0), 2 ) += 4.16666666667;
	sys.A( 1, vi + V2i(0,-1), 2 ) += 4.16666666667;
	sys.A( 1, vi + V2i(0,0), 2 ) += 4.16666666667;
	sys.A( 1, vi + V2i(1,-1), 2 ) += 4.16666666667;
	sys.A( 1, vi + V2i(1,0), 2 ) += 4.16666666667;
	sys.A( 1, vi + V2i(0,0), 2 ) += -4.16666666667;
	sys.A( 1, vi + V2i(0,1), 2 ) += -4.16666666667;
	sys.A( 1, vi + V2i(1,0), 2 ) += -4.16666666667;
	sys.A( 1, vi + V2i(1,1), 2 ) += -4.16666666667;

	// term 68 ----------------

	// term 69 ----------------

	// term 70 ----------------
	sys.A( 0, vi + V2i(-1,-1), 0 ) += -std::complex<double>(0.0, 4.166666666666666);
	sys.A( 0, vi + V2i(-1,0), 0 ) += -std::complex<double>(0.0, 4.166666666666666);
	sys.A( 0, vi + V2i(0,-1), 0 ) += -std::complex<double>(0.0, 4.166666666666666);
	sys.A( 0, vi + V2i(0,0), 0 ) += -std::complex<double>(0.0, 4.166666666666666);
	sys.A( 0, vi + V2i(0,-1), 0 ) += std::complex<double>(0.0, 4.166666666666666);
	sys.A( 0, vi + V2i(0,0), 0 ) += std::complex<double>(0.0, 4.166666666666666);
	sys.A( 0, vi + V2i(1,-1), 0 ) += std::complex<double>(0.0, 4.166666666666666);
	sys.A( 0, vi + V2i(1,0), 0 ) += std::complex<double>(0.0, 4.166666666666666);
	sys.A( 0, vi + V2i(-1,0), 0 ) += std::complex<double>(0.0, 4.166666666666666);
	sys.A( 0, vi + V2i(-1,1), 0 ) += std::complex<double>(0.0, 4.166666666666666);
	sys.A( 0, vi + V2i(0,0), 0 ) += std::complex<double>(0.0, 4.166666666666666);
	sys.A( 0, vi + V2i(0,1), 0 ) += std::complex<double>(0.0, 4.166666666666666);
	sys.A( 0, vi + V2i(0,0), 0 ) += -std::complex<double>(0.0, 4.166666666666666);
	sys.A( 0, vi + V2i(0,1), 0 ) += -std::complex<double>(0.0, 4.166666666666666);
	sys.A( 0, vi + V2i(1,0), 0 ) += -std::complex<double>(0.0, 4.166666666666666);
	sys.A( 0, vi + V2i(1,1), 0 ) += -std::complex<double>(0.0, 4.166666666666666);
	sys.A( 1, vi + V2i(-1,-1), 1 ) += -std::complex<double>(0.0, 0.8333333333333333);
	sys.A( 1, vi + V2i(-1,0), 1 ) += -std::complex<double>(0.0, 0.8333333333333333);
	sys.A( 1, vi + V2i(0,-1), 1 ) += -std::complex<double>(0.0, 0.8333333333333333);
	sys.A( 1, vi + V2i(0,0), 1 ) += -std::complex<double>(0.0, 0.8333333333333333);
	sys.A( 1, vi + V2i(0,-1), 1 ) += std::complex<double>(0.0, 0.8333333333333333);
	sys.A( 1, vi + V2i(0,0), 1 ) += std::complex<double>(0.0, 0.8333333333333333);
	sys.A( 1, vi + V2i(1,-1), 1 ) += std::complex<double>(0.0, 0.8333333333333333);
	sys.A( 1, vi + V2i(1,0), 1 ) += std::complex<double>(0.0, 0.8333333333333333);
	sys.A( 1, vi + V2i(-1,0), 1 ) += std::complex<double>(0.0, 0.8333333333333333);
	sys.A( 1, vi + V2i(-1,1), 1 ) += std::complex<double>(0.0, 0.8333333333333333);
	sys.A( 1, vi + V2i(0,0), 1 ) += std::complex<double>(0.0, 0.8333333333333333);
	sys.A( 1, vi + V2i(0,1), 1 ) += std::complex<double>(0.0, 0.8333333333333333);
	sys.A( 1, vi + V2i(0,0), 1 ) += -std::complex<double>(0.0, 0.8333333333333333);
	sys.A( 1, vi + V2i(0,1), 1 ) += -std::complex<double>(0.0, 0.8333333333333333);
	sys.A( 1, vi + V2i(1,0), 1 ) += -std::complex<double>(0.0, 0.8333333333333333);
	sys.A( 1, vi + V2i(1,1), 1 ) += -std::complex<double>(0.0, 0.8333333333333333);
	sys.A( 2, vi + V2i(-1,-1), 2 ) += -std::complex<double>(0.0, 5.0);
	sys.A( 2, vi + V2i(-1,0), 2 ) += -std::complex<double>(0.0, 5.0);
	sys.A( 2, vi + V2i(0,-1), 2 ) += -std::complex<double>(0.0, 5.0);
	sys.A( 2, vi + V2i(0,0), 2 ) += -std::complex<double>(0.0, 5.0);
	sys.A( 2, vi + V2i(0,-1), 2 ) += std::complex<double>(0.0, 5.0);
	sys.A( 2, vi + V2i(0,0), 2 ) += std::complex<double>(0.0, 5.0);
	sys.A( 2, vi + V2i(1,-1), 2 ) += std::complex<double>(0.0, 5.0);
	sys.A( 2, vi + V2i(1,0), 2 ) += std::complex<double>(0.0, 5.0);
	sys.A( 2, vi + V2i(-1,0), 2 ) += std::complex<double>(0.0, 5.0);
	sys.A( 2, vi + V2i(-1,1), 2 ) += std::complex<double>(0.0, 5.0);
	sys.A( 2, vi + V2i(0,0), 2 ) += std::complex<double>(0.0, 5.0);
	sys.A( 2, vi + V2i(0,1), 2 ) += std::complex<double>(0.0, 5.0);
	sys.A( 2, vi + V2i(0,0), 2 ) += -std::complex<double>(0.0, 5.0);
	sys.A( 2, vi + V2i(0,1), 2 ) += -std::complex<double>(0.0, 5.0);
	sys.A( 2, vi + V2i(1,0), 2 ) += -std::complex<double>(0.0, 5.0);
	sys.A( 2, vi + V2i(1,1), 2 ) += -std::complex<double>(0.0, 5.0);

	// term 71 ----------------

	// term 72 ----------------
	sys.A( 1, vi + V2i(0,-1), 2 ) += std::complex<double>(0.0, 3.333333333333333);
	sys.A( 1, vi + V2i(1,-1), 2 ) += -std::complex<double>(0.0, 3.333333333333333);
	sys.A( 1, vi + V2i(0,0), 2 ) += -std::complex<double>(0.0, 3.333333333333333);
	sys.A( 1, vi + V2i(1,0), 2 ) += std::complex<double>(0.0, 3.333333333333333);

	// term 73 ----------------

	// term 74 ----------------
	sys.A( 0, vi + V2i(0,-1), 0 ) += -16.6666666667;
	sys.A( 0, vi + V2i(0,0), 0 ) += 16.6666666667;
	sys.A( 0, vi + V2i(0,0), 0 ) += 16.6666666667;
	sys.A( 0, vi + V2i(0,1), 0 ) += -16.6666666667;
	sys.A( 1, vi + V2i(0,-1), 1 ) += -3.33333333333;
	sys.A( 1, vi + V2i(0,0), 1 ) += 3.33333333333;
	sys.A( 1, vi + V2i(0,0), 1 ) += 3.33333333333;
	sys.A( 1, vi + V2i(0,1), 1 ) += -3.33333333333;
	sys.A( 2, vi + V2i(0,-1), 2 ) += -20.0;
	sys.A( 2, vi + V2i(0,0), 2 ) += 20.0;
	sys.A( 2, vi + V2i(0,0), 2 ) += 20.0;
	sys.A( 2, vi + V2i(0,1), 2 ) += -20.0;

	// term 75 ----------------

	// term 76 ----------------
	sys.A( 1, vi + V2i(0,-2), 2 ) += -0.833333333333;
	sys.A( 1, vi + V2i(0,-1), 2 ) += -0.833333333333;
	sys.A( 1, vi + V2i(1,-2), 2 ) += -0.833333333333;
	sys.A( 1, vi + V2i(1,-1), 2 ) += -0.833333333333;
	sys.A( 1, vi + V2i(0,-1), 2 ) += 0.833333333333;
	sys.A( 1, vi + V2i(0,0), 2 ) += 0.833333333333;
	sys.A( 1, vi + V2i(1,-1), 2 ) += 0.833333333333;
	sys.A( 1, vi + V2i(1,0), 2 ) += 0.833333333333;
	sys.A( 1, vi + V2i(0,-1), 2 ) += 0.833333333333;
	sys.A( 1, vi + V2i(0,0), 2 ) += 0.833333333333;
	sys.A( 1, vi + V2i(1,-1), 2 ) += 0.833333333333;
	sys.A( 1, vi + V2i(1,0), 2 ) += 0.833333333333;
	sys.A( 1, vi + V2i(0,0), 2 ) += -0.833333333333;
	sys.A( 1, vi + V2i(0,1), 2 ) += -0.833333333333;
	sys.A( 1, vi + V2i(1,0), 2 ) += -0.833333333333;
	sys.A( 1, vi + V2i(1,1), 2 ) += -0.833333333333;

	// term 77 ----------------

	// term 78 ----------------

	// term 79 ----------------

	// term 80 ----------------

	// term 81 ----------------

	// term 82 ----------------

	// term 83 ----------------

	// term 84 ----------------

	// term 85 ----------------

	// term 86 ----------------

	// term 87 ----------------

	// term 88 ----------------

	// term 89 ----------------

	// term 90 ----------------

	// term 91 ----------------

	// term 92 ----------------

	// term 93 ----------------

	// term 94 ----------------

	// term 95 ----------------

	// term 96 ----------------

	// term 97 ----------------

	// term 98 ----------------

	// term 99 ----------------

	// term 100 ----------------
	sys.A( 2, vi + V2i(0,0), 0 ) += (0.204124145232*(-(10.0*fields.sigma_t->eval(sys.voxelToWorld(vd+V2d(-0.5, 0.5))))+
			(10.0*fields.sigma_t->eval(sys.voxelToWorld(vd+V2d(0.5, 0.5))))));
	sys.A( 2, vi + V2i(-1,0), 0 ) += (0.204124145232*(-(10.0*fields.sigma_t->eval(sys.voxelToWorld(vd+V2d(-0.5, 0.5))))+
			(10.0*fields.sigma_t->eval(sys.voxelToWorld(vd+V2d(0.5, 0.5))))));

	// term 101 ----------------
	sys.A( 0, vi + V2i(0,1), 1 ) += -(0.204124145232*(-(10.0*fields.sigma_t->eval(sys.voxelToWorld(vd+V2d(0.0, 0.5))))+
			(10.0*fields.sigma_t->eval(sys.voxelToWorld(vd+V2d(1.0, 0.5))))));
	sys.A( 0, vi + V2i(0,0), 1 ) += -(0.204124145232*(-(10.0*fields.sigma_t->eval(sys.voxelToWorld(vd+V2d(0.0, 0.5))))+
			(10.0*fields.sigma_t->eval(sys.voxelToWorld(vd+V2d(1.0, 0.5))))));

	// term 102 ----------------
	sys.A( 1, vi + V2i(0,0), 0 ) += -(0.204124145232*(-(10.0*fields.sigma_t->eval(sys.voxelToWorld(vd+V2d(0.0, 0.0))))+
			(10.0*fields.sigma_t->eval(sys.voxelToWorld(vd+V2d(1.0, 0.0))))));
	sys.A( 1, vi + V2i(0,-1), 0 ) += -(0.204124145232*(-(10.0*fields.sigma_t->eval(sys.voxelToWorld(vd+V2d(0.0, 0.0))))+
			(10.0*fields.sigma_t->eval(sys.voxelToWorld(vd+V2d(1.0, 0.0))))));

	// term 103 ----------------
	sys.A( 0, vi + V2i(1,0), 2 ) += (0.204124145232*(-(10.0*fields.sigma_t->eval(sys.voxelToWorld(vd+V2d(0.0, 0.5))))+
			(10.0*fields.sigma_t->eval(sys.voxelToWorld(vd+V2d(1.0, 0.5))))));
	sys.A( 0, vi + V2i(0,0), 2 ) += (0.204124145232*(-(10.0*fields.sigma_t->eval(sys.voxelToWorld(vd+V2d(0.0, 0.5))))+
			(10.0*fields.sigma_t->eval(sys.voxelToWorld(vd+V2d(1.0, 0.5))))));

	// term 104 ----------------
	sys.A( 2, vi + V2i(0,0), 0 ) += -(std::complex<double>(0.0, 0.2041241452319315)*(-(10.0*fields.sigma_t->eval(sys.voxelToWorld(vd+V2d(0.0, 0.0))))+
			(10.0*fields.sigma_t->eval(sys.voxelToWorld(vd+V2d(0.0, 1.0))))));
	sys.A( 2, vi + V2i(-1,0), 0 ) += -(std::complex<double>(0.0, 0.2041241452319315)*(-(10.0*fields.sigma_t->eval(sys.voxelToWorld(vd+V2d(0.0, 0.0))))+
			(10.0*fields.sigma_t->eval(sys.voxelToWorld(vd+V2d(0.0, 1.0))))));

	// term 105 ----------------
	sys.A( 0, vi + V2i(0,1), 1 ) += (std::complex<double>(0.0, 0.2041241452319315)*(-(10.0*fields.sigma_t->eval(sys.voxelToWorld(vd+V2d(0.5, 0.0))))+
			(10.0*fields.sigma_t->eval(sys.voxelToWorld(vd+V2d(0.5, 1.0))))));
	sys.A( 0, vi + V2i(0,0), 1 ) += (std::complex<double>(0.0, 0.2041241452319315)*(-(10.0*fields.sigma_t->eval(sys.voxelToWorld(vd+V2d(0.5, 0.0))))+
			(10.0*fields.sigma_t->eval(sys.voxelToWorld(vd+V2d(0.5, 1.0))))));

	// term 106 ----------------
	sys.A( 1, vi + V2i(0,0), 0 ) += -(std::complex<double>(0.0, 0.2041241452319315)*(-(10.0*fields.sigma_t->eval(sys.voxelToWorld(vd+V2d(0.5, -0.5))))+
			(10.0*fields.sigma_t->eval(sys.voxelToWorld(vd+V2d(0.5, 0.5))))));
	sys.A( 1, vi + V2i(0,-1), 0 ) += -(std::complex<double>(0.0, 0.2041241452319315)*(-(10.0*fields.sigma_t->eval(sys.voxelToWorld(vd+V2d(0.5, -0.5))))+
			(10.0*fields.sigma_t->eval(sys.voxelToWorld(vd+V2d(0.5, 0.5))))));

	// term 107 ----------------
	sys.A( 0, vi + V2i(1,0), 2 ) += (std::complex<double>(0.0, 0.2041241452319315)*(-(10.0*fields.sigma_t->eval(sys.voxelToWorld(vd+V2d(0.5, 0.0))))+
			(10.0*fields.sigma_t->eval(sys.voxelToWorld(vd+V2d(0.5, 1.0))))));
	sys.A( 0, vi + V2i(0,0), 2 ) += (std::complex<double>(0.0, 0.2041241452319315)*(-(10.0*fields.sigma_t->eval(sys.voxelToWorld(vd+V2d(0.5, 0.0))))+
			(10.0*fields.sigma_t->eval(sys.voxelToWorld(vd+V2d(0.5, 1.0))))));

	// term 108 ----------------

	// term 109 ----------------

	// term 110 ----------------
	sys.A( 0, vi + V2i(0,0), 0 ) += std::pow(fields.sigma_t->eval(sys.voxelToWorld(vd+V2d(0.5, 0.5))), 2);
	sys.A( 1, vi + V2i(0,0), 1 ) += std::pow(fields.sigma_t->eval(sys.voxelToWorld(vd+V2d(0.5, 0.0))), 2);
	sys.A( 2, vi + V2i(0,0), 2 ) += std::pow(fields.sigma_t->eval(sys.voxelToWorld(vd+V2d(0.0, 0.5))), 2);

	// term 111 ----------------
	sys.A( 2, vi + V2i(0,0), 0 ) += -(0.204124145232*(-(10.0*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(-0.5, 0.5))))+
			(10.0*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, 0.5)))))*fields.f_p->eval(0, 0, sys.voxelToWorld(vd+V2d(0.0, 0.5))));
	sys.A( 2, vi + V2i(-1,0), 0 ) += -(0.204124145232*(-(10.0*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(-0.5, 0.5))))+
			(10.0*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, 0.5)))))*fields.f_p->eval(0, 0, sys.voxelToWorld(vd+V2d(0.0, 0.5))));

	// term 112 ----------------
	sys.A( 2, vi + V2i(0,0), 0 ) += -(0.204124145232*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.0, 0.5)))*(-(10.0*fields.f_p->eval(0, 0, sys.voxelToWorld(vd+V2d(-0.5, 0.5))))+
			(10.0*fields.f_p->eval(0, 0, sys.voxelToWorld(vd+V2d(0.5, 0.5))))));
	sys.A( 2, vi + V2i(-1,0), 0 ) += -(0.204124145232*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.0, 0.5)))*(-(10.0*fields.f_p->eval(0, 0, sys.voxelToWorld(vd+V2d(-0.5, 0.5))))+
			(10.0*fields.f_p->eval(0, 0, sys.voxelToWorld(vd+V2d(0.5, 0.5))))));

	// term 113 ----------------
	sys.A( 2, vi + V2i(-1,0), 0 ) += (4.08248290464*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.0, 0.5)))*fields.f_p->eval(0, 0, sys.voxelToWorld(vd+V2d(0.0, 0.5))));
	sys.A( 2, vi + V2i(0,0), 0 ) += -(4.08248290464*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.0, 0.5)))*fields.f_p->eval(0, 0, sys.voxelToWorld(vd+V2d(0.0, 0.5))));

	// term 114 ----------------
	sys.A( 0, vi + V2i(0,1), 1 ) += (0.204124145232*(-(10.0*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.0, 0.5))))+
			(10.0*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(1.0, 0.5)))))*fields.f_p->eval(1, 0, sys.voxelToWorld(vd+V2d(0.5, 0.5))));
	sys.A( 0, vi + V2i(0,0), 1 ) += (0.204124145232*(-(10.0*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.0, 0.5))))+
			(10.0*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(1.0, 0.5)))))*fields.f_p->eval(1, 0, sys.voxelToWorld(vd+V2d(0.5, 0.5))));

	// term 115 ----------------
	sys.A( 0, vi + V2i(0,1), 1 ) += (0.204124145232*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, 0.5)))*(-(10.0*fields.f_p->eval(1, 0, sys.voxelToWorld(vd+V2d(0.0, 0.5))))+
			(10.0*fields.f_p->eval(1, 0, sys.voxelToWorld(vd+V2d(1.0, 0.5))))));
	sys.A( 0, vi + V2i(0,0), 1 ) += (0.204124145232*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, 0.5)))*(-(10.0*fields.f_p->eval(1, 0, sys.voxelToWorld(vd+V2d(0.0, 0.5))))+
			(10.0*fields.f_p->eval(1, 0, sys.voxelToWorld(vd+V2d(1.0, 0.5))))));

	// term 116 ----------------
	sys.A( 0, vi + V2i(-1,0), 1 ) += -(1.02062072616*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, 0.5)))*fields.f_p->eval(1, 0, sys.voxelToWorld(vd+V2d(0.5, 0.5))));
	sys.A( 0, vi + V2i(-1,1), 1 ) += -(1.02062072616*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, 0.5)))*fields.f_p->eval(1, 0, sys.voxelToWorld(vd+V2d(0.5, 0.5))));
	sys.A( 0, vi + V2i(0,0), 1 ) += -(1.02062072616*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, 0.5)))*fields.f_p->eval(1, 0, sys.voxelToWorld(vd+V2d(0.5, 0.5))));
	sys.A( 0, vi + V2i(0,1), 1 ) += -(1.02062072616*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, 0.5)))*fields.f_p->eval(1, 0, sys.voxelToWorld(vd+V2d(0.5, 0.5))));
	sys.A( 0, vi + V2i(0,0), 1 ) += (1.02062072616*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, 0.5)))*fields.f_p->eval(1, 0, sys.voxelToWorld(vd+V2d(0.5, 0.5))));
	sys.A( 0, vi + V2i(0,1), 1 ) += (1.02062072616*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, 0.5)))*fields.f_p->eval(1, 0, sys.voxelToWorld(vd+V2d(0.5, 0.5))));
	sys.A( 0, vi + V2i(1,0), 1 ) += (1.02062072616*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, 0.5)))*fields.f_p->eval(1, 0, sys.voxelToWorld(vd+V2d(0.5, 0.5))));
	sys.A( 0, vi + V2i(1,1), 1 ) += (1.02062072616*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, 0.5)))*fields.f_p->eval(1, 0, sys.voxelToWorld(vd+V2d(0.5, 0.5))));

	// term 117 ----------------
	sys.A( 1, vi + V2i(0,0), 0 ) += (0.204124145232*(-(10.0*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.0, 0.0))))+
			(10.0*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(1.0, 0.0)))))*fields.f_p->eval(0, 0, sys.voxelToWorld(vd+V2d(0.5, 0.0))));
	sys.A( 1, vi + V2i(0,-1), 0 ) += (0.204124145232*(-(10.0*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.0, 0.0))))+
			(10.0*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(1.0, 0.0)))))*fields.f_p->eval(0, 0, sys.voxelToWorld(vd+V2d(0.5, 0.0))));

	// term 118 ----------------
	sys.A( 1, vi + V2i(0,0), 0 ) += (0.204124145232*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, 0.0)))*(-(10.0*fields.f_p->eval(0, 0, sys.voxelToWorld(vd+V2d(0.0, 0.0))))+
			(10.0*fields.f_p->eval(0, 0, sys.voxelToWorld(vd+V2d(1.0, 0.0))))));
	sys.A( 1, vi + V2i(0,-1), 0 ) += (0.204124145232*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, 0.0)))*(-(10.0*fields.f_p->eval(0, 0, sys.voxelToWorld(vd+V2d(0.0, 0.0))))+
			(10.0*fields.f_p->eval(0, 0, sys.voxelToWorld(vd+V2d(1.0, 0.0))))));

	// term 119 ----------------
	sys.A( 1, vi + V2i(-1,-1), 0 ) += -(1.02062072616*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, 0.0)))*fields.f_p->eval(0, 0, sys.voxelToWorld(vd+V2d(0.5, 0.0))));
	sys.A( 1, vi + V2i(-1,0), 0 ) += -(1.02062072616*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, 0.0)))*fields.f_p->eval(0, 0, sys.voxelToWorld(vd+V2d(0.5, 0.0))));
	sys.A( 1, vi + V2i(0,-1), 0 ) += -(1.02062072616*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, 0.0)))*fields.f_p->eval(0, 0, sys.voxelToWorld(vd+V2d(0.5, 0.0))));
	sys.A( 1, vi + V2i(0,0), 0 ) += -(1.02062072616*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, 0.0)))*fields.f_p->eval(0, 0, sys.voxelToWorld(vd+V2d(0.5, 0.0))));
	sys.A( 1, vi + V2i(0,-1), 0 ) += (1.02062072616*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, 0.0)))*fields.f_p->eval(0, 0, sys.voxelToWorld(vd+V2d(0.5, 0.0))));
	sys.A( 1, vi + V2i(0,0), 0 ) += (1.02062072616*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, 0.0)))*fields.f_p->eval(0, 0, sys.voxelToWorld(vd+V2d(0.5, 0.0))));
	sys.A( 1, vi + V2i(1,-1), 0 ) += (1.02062072616*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, 0.0)))*fields.f_p->eval(0, 0, sys.voxelToWorld(vd+V2d(0.5, 0.0))));
	sys.A( 1, vi + V2i(1,0), 0 ) += (1.02062072616*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, 0.0)))*fields.f_p->eval(0, 0, sys.voxelToWorld(vd+V2d(0.5, 0.0))));

	// term 120 ----------------
	sys.A( 0, vi + V2i(1,0), 2 ) += -(0.204124145232*(-(10.0*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.0, 0.5))))+
			(10.0*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(1.0, 0.5)))))*fields.f_p->eval(1, 0, sys.voxelToWorld(vd+V2d(0.5, 0.5))));
	sys.A( 0, vi + V2i(0,0), 2 ) += -(0.204124145232*(-(10.0*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.0, 0.5))))+
			(10.0*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(1.0, 0.5)))))*fields.f_p->eval(1, 0, sys.voxelToWorld(vd+V2d(0.5, 0.5))));

	// term 121 ----------------
	sys.A( 0, vi + V2i(1,0), 2 ) += -(0.204124145232*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, 0.5)))*(-(10.0*fields.f_p->eval(1, 0, sys.voxelToWorld(vd+V2d(0.0, 0.5))))+
			(10.0*fields.f_p->eval(1, 0, sys.voxelToWorld(vd+V2d(1.0, 0.5))))));
	sys.A( 0, vi + V2i(0,0), 2 ) += -(0.204124145232*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, 0.5)))*(-(10.0*fields.f_p->eval(1, 0, sys.voxelToWorld(vd+V2d(0.0, 0.5))))+
			(10.0*fields.f_p->eval(1, 0, sys.voxelToWorld(vd+V2d(1.0, 0.5))))));

	// term 122 ----------------
	sys.A( 0, vi + V2i(0,0), 2 ) += (4.08248290464*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, 0.5)))*fields.f_p->eval(1, 0, sys.voxelToWorld(vd+V2d(0.5, 0.5))));
	sys.A( 0, vi + V2i(1,0), 2 ) += -(4.08248290464*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, 0.5)))*fields.f_p->eval(1, 0, sys.voxelToWorld(vd+V2d(0.5, 0.5))));

	// term 123 ----------------
	sys.A( 2, vi + V2i(0,0), 0 ) += (std::complex<double>(0.0, 0.2041241452319315)*(-(10.0*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.0, 0.0))))+
			(10.0*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.0, 1.0)))))*fields.f_p->eval(0, 0, sys.voxelToWorld(vd+V2d(0.0, 0.5))));
	sys.A( 2, vi + V2i(-1,0), 0 ) += (std::complex<double>(0.0, 0.2041241452319315)*(-(10.0*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.0, 0.0))))+
			(10.0*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.0, 1.0)))))*fields.f_p->eval(0, 0, sys.voxelToWorld(vd+V2d(0.0, 0.5))));

	// term 124 ----------------
	sys.A( 2, vi + V2i(0,0), 0 ) += (std::complex<double>(0.0, 0.2041241452319315)*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.0, 0.5)))*(-(10.0*fields.f_p->eval(0, 0, sys.voxelToWorld(vd+V2d(0.0, 0.0))))+
			(10.0*fields.f_p->eval(0, 0, sys.voxelToWorld(vd+V2d(0.0, 1.0))))));
	sys.A( 2, vi + V2i(-1,0), 0 ) += (std::complex<double>(0.0, 0.2041241452319315)*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.0, 0.5)))*(-(10.0*fields.f_p->eval(0, 0, sys.voxelToWorld(vd+V2d(0.0, 0.0))))+
			(10.0*fields.f_p->eval(0, 0, sys.voxelToWorld(vd+V2d(0.0, 1.0))))));

	// term 125 ----------------
	sys.A( 2, vi + V2i(-1,-1), 0 ) += -(std::complex<double>(0.0, 1.0206207261596576)*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.0, 0.5)))*fields.f_p->eval(0, 0, sys.voxelToWorld(vd+V2d(0.0, 0.5))));
	sys.A( 2, vi + V2i(-1,0), 0 ) += -(std::complex<double>(0.0, 1.0206207261596576)*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.0, 0.5)))*fields.f_p->eval(0, 0, sys.voxelToWorld(vd+V2d(0.0, 0.5))));
	sys.A( 2, vi + V2i(0,-1), 0 ) += -(std::complex<double>(0.0, 1.0206207261596576)*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.0, 0.5)))*fields.f_p->eval(0, 0, sys.voxelToWorld(vd+V2d(0.0, 0.5))));
	sys.A( 2, vi + V2i(0,0), 0 ) += -(std::complex<double>(0.0, 1.0206207261596576)*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.0, 0.5)))*fields.f_p->eval(0, 0, sys.voxelToWorld(vd+V2d(0.0, 0.5))));
	sys.A( 2, vi + V2i(-1,0), 0 ) += (std::complex<double>(0.0, 1.0206207261596576)*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.0, 0.5)))*fields.f_p->eval(0, 0, sys.voxelToWorld(vd+V2d(0.0, 0.5))));
	sys.A( 2, vi + V2i(-1,1), 0 ) += (std::complex<double>(0.0, 1.0206207261596576)*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.0, 0.5)))*fields.f_p->eval(0, 0, sys.voxelToWorld(vd+V2d(0.0, 0.5))));
	sys.A( 2, vi + V2i(0,0), 0 ) += (std::complex<double>(0.0, 1.0206207261596576)*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.0, 0.5)))*fields.f_p->eval(0, 0, sys.voxelToWorld(vd+V2d(0.0, 0.5))));
	sys.A( 2, vi + V2i(0,1), 0 ) += (std::complex<double>(0.0, 1.0206207261596576)*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.0, 0.5)))*fields.f_p->eval(0, 0, sys.voxelToWorld(vd+V2d(0.0, 0.5))));

	// term 126 ----------------
	sys.A( 0, vi + V2i(0,1), 1 ) += -(std::complex<double>(0.0, 0.2041241452319315)*(-(10.0*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, 0.0))))+
			(10.0*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, 1.0)))))*fields.f_p->eval(1, 0, sys.voxelToWorld(vd+V2d(0.5, 0.5))));
	sys.A( 0, vi + V2i(0,0), 1 ) += -(std::complex<double>(0.0, 0.2041241452319315)*(-(10.0*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, 0.0))))+
			(10.0*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, 1.0)))))*fields.f_p->eval(1, 0, sys.voxelToWorld(vd+V2d(0.5, 0.5))));

	// term 127 ----------------
	sys.A( 0, vi + V2i(0,1), 1 ) += -(std::complex<double>(0.0, 0.2041241452319315)*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, 0.5)))*(-(10.0*fields.f_p->eval(1, 0, sys.voxelToWorld(vd+V2d(0.5, 0.0))))+
			(10.0*fields.f_p->eval(1, 0, sys.voxelToWorld(vd+V2d(0.5, 1.0))))));
	sys.A( 0, vi + V2i(0,0), 1 ) += -(std::complex<double>(0.0, 0.2041241452319315)*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, 0.5)))*(-(10.0*fields.f_p->eval(1, 0, sys.voxelToWorld(vd+V2d(0.5, 0.0))))+
			(10.0*fields.f_p->eval(1, 0, sys.voxelToWorld(vd+V2d(0.5, 1.0))))));

	// term 128 ----------------
	sys.A( 0, vi + V2i(0,0), 1 ) += (std::complex<double>(0.0, 4.08248290463863)*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, 0.5)))*fields.f_p->eval(1, 0, sys.voxelToWorld(vd+V2d(0.5, 0.5))));
	sys.A( 0, vi + V2i(0,1), 1 ) += -(std::complex<double>(0.0, 4.08248290463863)*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, 0.5)))*fields.f_p->eval(1, 0, sys.voxelToWorld(vd+V2d(0.5, 0.5))));

	// term 129 ----------------
	sys.A( 1, vi + V2i(0,0), 0 ) += (std::complex<double>(0.0, 0.2041241452319315)*(-(10.0*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, -0.5))))+
			(10.0*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, 0.5)))))*fields.f_p->eval(0, 0, sys.voxelToWorld(vd+V2d(0.5, 0.0))));
	sys.A( 1, vi + V2i(0,-1), 0 ) += (std::complex<double>(0.0, 0.2041241452319315)*(-(10.0*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, -0.5))))+
			(10.0*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, 0.5)))))*fields.f_p->eval(0, 0, sys.voxelToWorld(vd+V2d(0.5, 0.0))));

	// term 130 ----------------
	sys.A( 1, vi + V2i(0,0), 0 ) += (std::complex<double>(0.0, 0.2041241452319315)*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, 0.0)))*(-(10.0*fields.f_p->eval(0, 0, sys.voxelToWorld(vd+V2d(0.5, -0.5))))+
			(10.0*fields.f_p->eval(0, 0, sys.voxelToWorld(vd+V2d(0.5, 0.5))))));
	sys.A( 1, vi + V2i(0,-1), 0 ) += (std::complex<double>(0.0, 0.2041241452319315)*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, 0.0)))*(-(10.0*fields.f_p->eval(0, 0, sys.voxelToWorld(vd+V2d(0.5, -0.5))))+
			(10.0*fields.f_p->eval(0, 0, sys.voxelToWorld(vd+V2d(0.5, 0.5))))));

	// term 131 ----------------
	sys.A( 1, vi + V2i(0,-1), 0 ) += -(std::complex<double>(0.0, 4.08248290463863)*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, 0.0)))*fields.f_p->eval(0, 0, sys.voxelToWorld(vd+V2d(0.5, 0.0))));
	sys.A( 1, vi + V2i(0,0), 0 ) += (std::complex<double>(0.0, 4.08248290463863)*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, 0.0)))*fields.f_p->eval(0, 0, sys.voxelToWorld(vd+V2d(0.5, 0.0))));

	// term 132 ----------------
	sys.A( 0, vi + V2i(1,0), 2 ) += -(std::complex<double>(0.0, 0.2041241452319315)*(-(10.0*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, 0.0))))+
			(10.0*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, 1.0)))))*fields.f_p->eval(1, 0, sys.voxelToWorld(vd+V2d(0.5, 0.5))));
	sys.A( 0, vi + V2i(0,0), 2 ) += -(std::complex<double>(0.0, 0.2041241452319315)*(-(10.0*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, 0.0))))+
			(10.0*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, 1.0)))))*fields.f_p->eval(1, 0, sys.voxelToWorld(vd+V2d(0.5, 0.5))));

	// term 133 ----------------
	sys.A( 0, vi + V2i(1,0), 2 ) += -(std::complex<double>(0.0, 0.2041241452319315)*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, 0.5)))*(-(10.0*fields.f_p->eval(1, 0, sys.voxelToWorld(vd+V2d(0.5, 0.0))))+
			(10.0*fields.f_p->eval(1, 0, sys.voxelToWorld(vd+V2d(0.5, 1.0))))));
	sys.A( 0, vi + V2i(0,0), 2 ) += -(std::complex<double>(0.0, 0.2041241452319315)*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, 0.5)))*(-(10.0*fields.f_p->eval(1, 0, sys.voxelToWorld(vd+V2d(0.5, 0.0))))+
			(10.0*fields.f_p->eval(1, 0, sys.voxelToWorld(vd+V2d(0.5, 1.0))))));

	// term 134 ----------------
	sys.A( 0, vi + V2i(0,-1), 2 ) += (std::complex<double>(0.0, 1.0206207261596576)*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, 0.5)))*fields.f_p->eval(1, 0, sys.voxelToWorld(vd+V2d(0.5, 0.5))));
	sys.A( 0, vi + V2i(0,0), 2 ) += (std::complex<double>(0.0, 1.0206207261596576)*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, 0.5)))*fields.f_p->eval(1, 0, sys.voxelToWorld(vd+V2d(0.5, 0.5))));
	sys.A( 0, vi + V2i(1,-1), 2 ) += (std::complex<double>(0.0, 1.0206207261596576)*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, 0.5)))*fields.f_p->eval(1, 0, sys.voxelToWorld(vd+V2d(0.5, 0.5))));
	sys.A( 0, vi + V2i(1,0), 2 ) += (std::complex<double>(0.0, 1.0206207261596576)*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, 0.5)))*fields.f_p->eval(1, 0, sys.voxelToWorld(vd+V2d(0.5, 0.5))));
	sys.A( 0, vi + V2i(0,0), 2 ) += -(std::complex<double>(0.0, 1.0206207261596576)*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, 0.5)))*fields.f_p->eval(1, 0, sys.voxelToWorld(vd+V2d(0.5, 0.5))));
	sys.A( 0, vi + V2i(0,1), 2 ) += -(std::complex<double>(0.0, 1.0206207261596576)*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, 0.5)))*fields.f_p->eval(1, 0, sys.voxelToWorld(vd+V2d(0.5, 0.5))));
	sys.A( 0, vi + V2i(1,0), 2 ) += -(std::complex<double>(0.0, 1.0206207261596576)*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, 0.5)))*fields.f_p->eval(1, 0, sys.voxelToWorld(vd+V2d(0.5, 0.5))));
	sys.A( 0, vi + V2i(1,1), 2 ) += -(std::complex<double>(0.0, 1.0206207261596576)*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, 0.5)))*fields.f_p->eval(1, 0, sys.voxelToWorld(vd+V2d(0.5, 0.5))));

	// term 135 ----------------

	// term 136 ----------------

	// term 137 ----------------

	// term 138 ----------------

	// term 139 ----------------

	// term 140 ----------------

	// term 141 ----------------
	sys.A( 0, vi + V2i(0,0), 0 ) += -(fields.sigma_t->eval(sys.voxelToWorld(vd+V2d(0.5, 0.5)))*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, 0.5)))*fields.f_p->eval(0, 0, sys.voxelToWorld(vd+V2d(0.5, 0.5))));
	sys.A( 1, vi + V2i(0,0), 1 ) += -(fields.sigma_t->eval(sys.voxelToWorld(vd+V2d(0.5, 0.0)))*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, 0.0)))*fields.f_p->eval(1, 0, sys.voxelToWorld(vd+V2d(0.5, 0.0))));
	sys.A( 2, vi + V2i(0,0), 2 ) += -(fields.sigma_t->eval(sys.voxelToWorld(vd+V2d(0.0, 0.5)))*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.0, 0.5)))*fields.f_p->eval(1, 0, sys.voxelToWorld(vd+V2d(0.0, 0.5))));

	// term 142 ----------------
	sys.b(2) += (0.5*0.816496580928*((-10.0*fields.q->eval(0, 0, sys.voxelToWorld(vd+V2d(-0.5, 0.5))))+
			(10.0*fields.q->eval(0, 0, sys.voxelToWorld(vd+V2d(0.5, 0.5))))));

	// term 143 ----------------
	sys.b(0) += -(0.5*0.816496580928*((-10.0*fields.q->eval(1, -1, sys.voxelToWorld(vd+V2d(0.0, 0.5))))+
			(10.0*fields.q->eval(1, -1, sys.voxelToWorld(vd+V2d(1.0, 0.5))))));

	// term 144 ----------------
	sys.b(1) += -(0.5*0.816496580928*((-10.0*fields.q->eval(0, 0, sys.voxelToWorld(vd+V2d(0.0, 0.0))))+
			(10.0*fields.q->eval(0, 0, sys.voxelToWorld(vd+V2d(1.0, 0.0))))));

	// term 145 ----------------
	sys.b(0) += (0.5*0.816496580928*((-10.0*fields.q->eval(1, 1, sys.voxelToWorld(vd+V2d(0.0, 0.5))))+
			(10.0*fields.q->eval(1, 1, sys.voxelToWorld(vd+V2d(1.0, 0.5))))));

	// term 146 ----------------
	sys.b(2) += -(std::complex<double>(0.0, 0.5)*0.816496580928*((-10.0*fields.q->eval(0, 0, sys.voxelToWorld(vd+V2d(0.0, 0.0))))+
			(10.0*fields.q->eval(0, 0, sys.voxelToWorld(vd+V2d(0.0, 1.0))))));

	// term 147 ----------------
	sys.b(0) += (std::complex<double>(0.0, 0.5)*0.816496580928*((-10.0*fields.q->eval(1, -1, sys.voxelToWorld(vd+V2d(0.5, 0.0))))+
			(10.0*fields.q->eval(1, -1, sys.voxelToWorld(vd+V2d(0.5, 1.0))))));

	// term 148 ----------------
	sys.b(1) += -(std::complex<double>(0.0, 0.5)*0.816496580928*((-10.0*fields.q->eval(0, 0, sys.voxelToWorld(vd+V2d(0.5, -0.5))))+
			(10.0*fields.q->eval(0, 0, sys.voxelToWorld(vd+V2d(0.5, 0.5))))));

	// term 149 ----------------
	sys.b(0) += (std::complex<double>(0.0, 0.5)*0.816496580928*((-10.0*fields.q->eval(1, 1, sys.voxelToWorld(vd+V2d(0.5, 0.0))))+
			(10.0*fields.q->eval(1, 1, sys.voxelToWorld(vd+V2d(0.5, 1.0))))));

	// term 150 ----------------

	// term 151 ----------------

	// term 152 ----------------
	sys.b(0) += (fields.sigma_t->eval(sys.voxelToWorld(vd+V2d(0.5, 0.5)))*fields.q->eval(0, 0, sys.voxelToWorld(vd+V2d(0.5, 0.5))));
	sys.b(1) += (fields.sigma_t->eval(sys.voxelToWorld(vd+V2d(0.5, 0.0)))*fields.q->eval(1, -1, sys.voxelToWorld(vd+V2d(0.5, 0.0))));
	sys.b(2) += (fields.sigma_t->eval(sys.voxelToWorld(vd+V2d(0.0, 0.5)))*fields.q->eval(1, 1, sys.voxelToWorld(vd+V2d(0.0, 0.5))));
}