﻿#pragma once

#define _CRTDBG_MAP_ALLOC  
#include <stdlib.h>  
#include <crtdbg.h>  

/* constants */
const long long Bill = 1000000000;  // one billion = 10^9
const long long quintillion = Bill*Bill;  // one quintillion = 10^18

/* now using boost miltiprecision library as well as MPIR to handle bigints. 
see http://www.boost.org/doc/libs/1_65_1/libs/multiprecision/doc/html/boost_multiprecision */

/* putting the using namespace in the header file is considered bad practice,
but not putting it here would prevent the mpz_int type integers in the declarations
below from compiling.*/
using namespace boost::multiprecision;

/* functions that use bigints (mpz_int)*/

void ShowLargeXY(std::string x, std::string y, mpz_int Bi_H1, mpz_int Bi_K1,
	bool sol, const std::string eqX, const std::string eqY);
bool ShowHomoSols(int type, mpz_int Bi_H1, mpz_int Bi_K1, long long s, long long T,
	long long MagnifY, std::string eqX, std::string eqY);
void ShowLargeNumber(const mpz_int Bi_Nbr);
long long DivLargeNumber(const mpz_int Bi_Nbr, long long Coef, mpz_int *Bi_Dest);
long long tDivLargeNumber(const mpz_int n, long long d, mpz_int *q);
long long DoublePrecToLong(const mpz_int x);
long long DivDoublePrec(const mpz_int Dp_Dividend, const mpz_int Dp_Divisor);
void DivideDoublePrecLong(const mpz_int Dp_Dividend, const mpz_int Dp_Divisor, mpz_int *Dp_Quotient);
void gcd(const mpz_int Dp_Nbr1, const mpz_int Dp_Nbr2, mpz_int *Dp_Gcd);
long long gcd(long long M, long long N);    // same name, does same job, but for normal integers
std::string numToStr(const mpz_int Dp_Nbr);
std::string numToStr(long long num);        // same name, does same job, but for normal integers


/* other functions that need a separate declaration */

std::string par(long long num);
void ShowLin(long long D, long long E, long long F, std::string x, std::string y);
bool ContFrac(const mpz_int Dp_A, int type, int SqrtSign, long long s, long long T,
	long long MagnifY, long long A);
void ShowAllLargeSolutions();
void ShowEq(long long A, long long B, long long C, long long D, long long E, long long F,
	std::string x, std::string y);
void SolContFrac(long long H, long long T, long long A, long long B, long long C, std::string SCFstr);

long long MultMod(long long factor1, long long factor2, long long Mod);
long long ModInv(long long Val, long long Mod);
long long ModPow(long long Base, long long Exp, long long Mod);
void ShowBigEq(mpz_int Dp_A, mpz_int Dp_B, mpz_int Dp_C, std::string x, std::string y);
void GetRoot(mpz_int Dp_A, mpz_int Dp_B, mpz_int Dp_C);
//int Compare(const mpz_int Bi_array, const mpz_int Bi_K1);

/* enumerated variable used to classify type of equation */
enum equation_class {
	linear,                // A = B = C = 0.
	simple_hyperbolic,     // A = C = 0; B ≠ 0. (implies B ^ 2 - 4AC > 0)
	elliptical,            // B^2 - 4AC < 0.
	parabolic,             // B ^ 2 - 4AC = 0
	hyperbolic_homog,      // B ^ 2 - 4AC > 0,  D = 0, E = 0
	hyperbolic_gen,        // B ^ 2 - 4AC > 0, not in other hyperbolic classes above
	no_soln                // fails tests that check a solution exists.
						   // note, even if it passes the tests there may still be no solution
};