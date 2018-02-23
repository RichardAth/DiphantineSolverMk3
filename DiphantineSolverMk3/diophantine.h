#pragma once

#define _CRTDBG_MAP_ALLOC  
#include <crtdbg.h>  

/* constants */
const long long Bill = 1000000000;  // one billion = 10^9
const long long quintillion = Bill*Bill;  // one quintillion = 10^18

/* now using boost miltiprecision library as well as MPIR to handle bigints. 
see http://www.boost.org/doc/libs/1_65_1/libs/multiprecision/doc/html/boost_multiprecision */

typedef boost::multiprecision::mpz_int bigint;
/* access underlying mpz_t inside an bigint */
#define ZT(a) a.backend().data()

/* functions that use bigints (mpz_t AKA mpz_int)*/

void ShowLargeXY(const std::string &x, const std::string &y, const bigint &Bi_H1, const bigint &Bi_K1,
	bool sol, const std::string &eqX, const std::string &eqY);
bool ShowHomoSols(const int type, const bigint &H1, const bigint &K1, long long s, const bigint &T,
	const bigint &MagnifY, const std::string &eqX, const std::string &eqY, const bigint &F);
void ShowLargeNumber(const bigint &Bi_Nbr);
long long DivLargeNumberRem(const bigint &Bi_Nbr, long long Coef, bigint *Bi_Dest);
long long tDivLargeNumber(const bigint &n, const bigint &d, bigint &q);
long long MulPrToLong(const bigint &x);
long long DivLargeNumberLL(const bigint &Dp_Dividend, const bigint &Dp_Divisor);
void DivLargeNumber(const bigint &Dp_Dividend, const bigint &Dp_Divisor, bigint *Dp_Quotient);
void gcd(const bigint &Dp_Nbr1, const bigint &Dp_Nbr2, bigint &Dp_Gcd);
long long gcd(long long M, long long N);    // same name, does same job, but for normal integers
std::string numToStr(const bigint &Dp_Nbr);
std::string numToStr(long long num);        // same name, does same job, but for normal integers


/* other functions that need a separate declaration */

std::string par(long long num);
std::string par(const bigint &num);
void ShowLin(long long D, long long E, long long F, std::string x, std::string y);
void ShowLin(bigint D, bigint E, bigint F, std::string x, std::string y);
bool ContFrac(const bigint &Dp_A, int type, int SqrtSign, long long s, const bigint &T,
	bigint MagnifY, bigint A, const bigint &Disc, const bigint &F);
void ShowAllLargeSolutions();
void ShowEq(const bigint A, const bigint B, const bigint C, const bigint D, 
	const bigint E, const bigint F, std::string x, std::string y);
void SolContFrac(const bigint &H, long long T, bigint A, const long long B, bigint C, 
	std::string SCFstr);

long long MultMod(long long factor1, long long factor2, long long Mod);
long long MultMod(const bigint &factor1, const bigint &factor2, bigint Mod);
long long ModInv(const long long Val, const long long Mod);
void ModInv(bigint* op, const bigint &Val, const bigint &Mod);
long long ModPow(long long Base, long long Exp, long long Mod);
void ShowBigEq(const bigint &Dp_A, const bigint &Dp_B, const bigint &Dp_C, 
	const std::string &x, const std::string &y);
void GetRoot(const bigint &Dp_A, const bigint &Dp_B, const bigint &Dp_C, bigint *Disc);
//int Compare(const bigint Bi_array, const bigint Bi_K1);
void adjustGandK(long long &G, long long &K);

/* enumerated variable used to classify type of equation */
enum equation_class {
	linear,                // A = B = C = 0.
	simple_hyperbolic,     // A = C = 0; B ≠ 0. (implies B^2 - 4AC > 0)
	elliptical,            // B^2 - 4AC < 0.
	parabolic,             // B^2 - 4AC = 0
	hyperbolic_homog,      // B^2 - 4AC > 0,  D = 0, E = 0
	hyperbolic_disc_ps,	   // B^2 - 4AC > 0 and (B^2 - 4AC) is a perfect square
	hyperbolic_gen,        // B^2 - 4AC > 0, not in other hyperbolic classes above
	no_soln                // fails tests that check a solution exists.
						   // note, even if it passes the tests there may still be no solution
};