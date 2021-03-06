﻿#define _CRTDBG_MAP_ALLOC  
#include <stdlib.h>  
#include <crtdbg.h> 
#include <safeint.h>

#include "stdafx.h"
#include <windows.h> 
#include <set>
//#define temp 1

// Diophantine Quadratic Equation Solver
// ax^2 + bxy + cy^2 + dx + ey + f = 0 (unknowns x,y integer numbers)
//
// Written by Dario Alejandro Alpern (Buenos Aires - Argentina)
// (Last updated December 15th, 2003)
//
// No part of this code can be used for commercial purposes without
// the written consent from the author. Otherwise it can be used freely.
// 
// converted fron Java to C++ by R Atherton, Octobler 2017
// because of problems with java not being supported by web browsers.


/* following used to detect memory leakage */
_CrtMemState memState1;          // initial memory state
_CrtMemState memState2;          // current memory state
_CrtMemState memStatDiff;        // difference between initial and current memory state

std::string txtStr;

/* ² is U00b2 in unicode (decimal 178)*/
TCHAR sqchr = L'²';
const std::wstring sq2 =  L"²";        // would be nice to use, but causes problems 

const std::string sq = "^2";
const std::string msg = "There are no solutions !!!\n";
std::string g_s_number;
bool also = false, ExchXY = false, teach = false;
bool allSolsFound;
const std::string divgcd = "Dividing the equation by the GCD we obtain: \n";

long long g_A, g_B, g_C, g_D, g_E, g_F;  // coefficients for x², xy, y² x, y, and constant

// set by SolveEquation
bigint g_CY1, g_CY0;  // used by SolveDiscIsSq, ShowHomoSols & ContFrac

bigint g_A2, g_B2;   // note: g_A2, g_B2 set by ContFrac, used by SolveEquation

bigint Bi_Disc;     
unsigned long long SqrtDisc;  // square root of Disc, set by GetRoot

/* large numbers; shown by Bi_ prefix */
bigint Bi_H1, Bi_H2, Bi_K2;
bigint Bi_NUM, Bi_DEN;              // set by GetRoot

int NbrSols = 0, NbrCo;
int digitsInGroup = 6;     // used in ShowLargeNumber

/* typedef for vector-type arrays */
//typedef  std::vector <void *> ArrayLongs;

//ArrayLongs sortedSolsXv;
//ArrayLongs sortedSolsYv;

uint32_t *primeFlags = NULL; /* 1 bit for each odd number, 1 = not prime, 0 = prime */
unsigned long long *primeList = NULL;
unsigned int prime_list_count = 0;
unsigned long long int primeListMax = 0;
unsigned long long int max_prime = 100000000;  // arbitrary limit,

//void listLargeSolutions();

/* print a string variable. This is a relic from the original java program */
//void w(std::string texto) {
//	std::cout << texto;
//}
void CheckForLeaks(void) {
#ifdef _DEBUG
	_CrtMemCheckpoint(&memState2); /* check for memory leakage*/
	if (_CrtMemDifference(&memStatDiff, &memState1, &memState2)) {
		/* Heap memory state has changed. may be a memory leak */
		std::cerr << "** there seems to have been a memory leak";
		_CrtMemDumpStatistics(&memStatDiff);  // show changes
		//_CrtMemDumpAllObjectsSince(&memState1);
		//_CrtDumpMemoryLeaks();
		_CrtMemCheckpoint(&memState1);   // update checkpoint;
	}

#endif
}

/* calculate integer square root */
unsigned long long llSqrt(unsigned long long num) {
	unsigned long long num1 = 0;
	unsigned long long num2 = (long long)65536 * (long long)32768;  /* 2^31 */
	while (num2 != 0) {
		if ((num1 + num2)*(num1 + num2) <= num) {
			num1 += num2;
		}
		num2 /= 2;
	}
	return num1;
}

/* get sqare root of ext prec. num. If num > 128 bits will get overflow.*/
unsigned long long llSqrt(const bigint &num) {
	mpz_t rv;
	unsigned long long llrv;

	assert(num >= 0);  // god knows what would happen if num were -ve
	mpz_init(rv);
	mpz_sqrt(rv, ZT(num));  // rv = truncated integer part of the square root
	llrv = MulPrToLong(rv);  // convert to long long (check for overflow)
	mpz_clear(rv);
	return llrv;
}

/* calculate Greatest Common Divisor */
long long gcd(long long M, long long N) {
	long long P = M;
	long long Q = N;
	while (P != 0) {
		long long R = Q%P;
		Q = P;
		P = R;
	}
	return abs(Q);
}

/* calculate Greatest Common Divisor for Bigints*/
void gcd(const bigint &Dp_Nbr1, const bigint &Dp_Nbr2, bigint &Dp_Gcd) {
	mpz_t gcd;
	mpz_init(gcd);
	mpz_gcd(gcd, ZT(Dp_Nbr1), ZT(Dp_Nbr2));
	Dp_Gcd = gcd;
	mpz_clear(gcd);
}

/* convert double to a string */
//std::string doubToStr(double num) {
//	char numbr[21];
//	sprintf_s(numbr, sizeof(numbr), "%g", num);
//	number = numbr;   // copy to C++ string 
//	return number;
//}

void showAlso() {
	if (also) {
		printf("and also: ");
	}
	else {
		also = true;
	}
}

/* print values of X and Y */
void ShowXY(long long X, long long Y) {
	showAlso();
	if (!ExchXY) {
		//w("x = " + numToStr(X));
		printf("x = %lld  y = %lld \n", X, Y);
		//w("\ny = " + numToStr(Y) + "\n");
	}
	else {
		//w("x = " + numToStr(Y));
		printf("x = %lld  y = %lld \n", Y, X);
		//w("\ny = " + numToStr(X) + "\n");
	}
}
void ShowXY(bigint X, bigint Y) {
	showAlso();
	if (!ExchXY) {
		std::cout<< "x = " << X << "  y = " << Y << "\n";
	}
	else {
		std::cout << "x = " << Y << "  y = " << X << "\n";
	}
}

/* calculate and print values of X and Y. */
void ShowX1Y1(long long X1, long long Y1, long long A, long long B, long long D,
	long long D1, long long E1) {
	long long X, Y;

	if ((Y1 - E1) % D1 == 0) {
		Y = (Y1 - E1) / D1;
		if ((X1 - B*Y - D) % (2 * A) == 0) {
			X = (X1 - B*Y - D) / (2 * A);
			ShowXY(X, Y);          // display solution
		}
		if (X1 != 0 && (-X1 - B*Y - D) % (2 * A) == 0) {
			X = (-X1 - B*Y - D) / (2 * A);
			ShowXY(X, Y);           // display solution
		}
	}

	if (Y1 != 0 && (-Y1 - E1) % D1 == 0) {
		Y = (-Y1 - E1) / D1;
		if ((X1 - B*Y - D) % (2 * A) == 0) {
			X = (X1 - B*Y - D) / (2 * A);
			ShowXY(X, Y);           // display solution
		}
		if (X1 != 0 && (-X1 - B*Y - D) % (2 * A) == 0) {
			X = (-X1 - B*Y - D) / (2 * A);
			ShowXY(X, Y);           // display solution
		}
	}
}
void ShowX1Y1(bigint X1, bigint Y1, bigint A, bigint B, bigint D,
	bigint D1, bigint E1) {
	bigint X, Y;

	if ((Y1 - E1) % D1 == 0) {
		Y = (Y1 - E1) / D1;
		if ((X1 - B*Y - D) % (2 * A) == 0) {
			X = (X1 - B*Y - D) / (2 * A);
			ShowXY(X, Y);          // display solution
		}
		if (X1 != 0 && (-X1 - B*Y - D) % (2 * A) == 0) {
			X = (-X1 - B*Y - D) / (2 * A);
			ShowXY(X, Y);           // display solution
		}
	}

	if (Y1 != 0 && (-Y1 - E1) % D1 == 0) {
		Y = (-Y1 - E1) / D1;
		if ((X1 - B*Y - D) % (2 * A) == 0) {
			X = (X1 - B*Y - D) / (2 * A);
			ShowXY(X, Y);           // display solution
		}
		if (X1 != 0 && (-X1 - B*Y - D) % (2 * A) == 0) {
			X = (-X1 - B*Y - D) / (2 * A);
			ShowXY(X, Y);           // display solution
		}
	}
}

/* output value of num. Uses global variable txtStr.
iff t is odd, use '+' for +ve numbers.
if t =2 result is left in txtStr and not output. (this feature does not seem to be used anywhere)
if num is zero do nothing */
int Show(long long num, std::string str, int t) {
	if (t == 2) {
		txtStr = "";    // 1st output of sequence; clear contents of txtStr
	}
	if (num != 0) {
		std::string str1 = "";
		if ((t & 1) != 0 && num>0) {
			str1 += " +";
		}
		if (num<0) {
			str1 += " -";
		}
		if (abs(num) != 1) {
			str1 += numToStr(abs(num));
		}

		str1 += str;    // add suffix if any
		if ((t & 2) != 0) {
			txtStr += str1;
		}
		else {
			std::cout << str1;
		}
		return (t | 1);  // set bit 0, return type for next time
	}
	return t;
}
int Show(bigint num, std::string str, int t) {
	if (t == 2) {
		txtStr = "";    // 1st output of sequence; clear contents of txtStr
	}
	if (num != 0) {
		std::string str1 = "";
		if ((t & 1) != 0 && num>0) {
			str1 += " +";
		}
		if (num<0) {
			str1 += " -";
		}
		if (abs(num) != 1) {
			str1 += numToStr(abs(num));
		}

		str1 += str;    // add suffix if any
		if ((t & 2) != 0) {
			txtStr += str1;
		}
		else {
			std::cout << str1;
		}
		return (t | 1);  // set bit 0, return type for next time
	}
	return t;
}

/* Uses global variable txtStr. 
if t =2 result is left in txtStr and not output. (this feature does not seem to be used anywhere)*/
void Show1(long long num, int t) {
	char u = Show(num, "", t);
	if ((u & 1) == 0 || abs(num) == 1) {
		if ((u & 2) != 0) {
			txtStr += numToStr(abs(num));  // output left in txtStr, not printed
		}
		else {
			//w("" + numToStr(abs(num)));
			std::cout << abs(num);
		}
	}
}
void Show1(bigint num, int t) {
	char u = Show(num, "", t);
	if ((u & 1) == 0 || abs(num) == 1) {
		if ((u & 2) != 0) {
			txtStr += numToStr(abs(num));  // output left in txtStr, not printed
		}
		else {
			//w("" + numToStr(abs(num)));
			std::cout << abs(num);
		}
	}
}

/* print values of D, E and F */
void ShowLin(long long D, long long E, long long F, std::string x, std::string y) {
	int t = Show(D, x, 0);
	t = Show(E, y, t);
	Show1(F, t);
}
void ShowLin(bigint D, bigint E, bigint F, std::string x, std::string y) {
	int t = Show(D, x, 0);
	t = Show(E, y, t);
	Show1(F, t);
}

/* print "ax^2 + Bxy + Cy^2 +Dx + Ey + F" */
void ShowEq(const bigint A, const bigint B, const bigint C, const bigint D, 
	const bigint E, const bigint F, std::string x, std::string y) {
	int t = Show(A, x + sq, 0);
	t = Show(B, x + y, t);
	t = Show(C, y + sq, t);
	t = Show(D, x, t);
	t = Show(E, y, t);
	Show1(F, t);
}

/* print " so there are no integer solutions" and set also to 'true' */
void NoSol() {
	printf(" so there are no integer solutions.\n");
	also = true;
}

void NoGcd(bigint F) {
	if (teach) {
		std::cout << "This gcd is not a divisor of the constant term (" << F << ")," << "\n";
	}
	NoSol();
}

/* divide num by den, using floor division. 
(standard / operator uses truncation division.
DivLargeNumberLL does the same, but num and den are bigint numbers.*/
long long floordiv(long long num, long long den) {
	if ((num<0 && den>0 || num>0 && den<0) && num%den != 0) {
		return num / den - 1;
	}
	return num / den;
}

// commented out because it's not actually used
//long long ceildiv(long long num, long long den) {
//	if ((num>0 && den>0 || num<0 && den<0) && num%den != 0) {
//		return num / den + 1;
//	}
//	return num / den;
//}
/* q = n/d. return quotient as a bigint. Uses floor division.
i.e. rounds q down towards -infinity. The only difference from
< DivLargeNumberRem> is that this function does not return a value but
< DivLargeNumberRem> returns the remainder as a long long. */
void DivLargeNumber(const bigint &n, const bigint &d, bigint *q) {
	mpz_t qq;

	if (d == 0) {
		fprintf(stderr, "** divide by zero error\n");
		throw std::invalid_argument("divide by zero");
	}
	mpz_init(qq);
	mpz_fdiv_q(qq, ZT(n), ZT(d));          // note use of floor division
	*q = qq;
	mpz_clear(qq);
}

/* q = n/d    return value is remainder. Uses floor division
i.e. rounds q down towards -infinity, and r will have the same sign as d.
the 3 different types of division only give different results when
the remainder is non-zero and n or d (or both) are negative.
NB. The choice of floor or truncation type division is sometimes
critical, and the correct choice is not obvious. */
long long DivLargeNumberRem(const bigint &n, long long d, bigint *q) {
	long long remainder;
	mpz_t mpr, mpd, nsave, qq;

	if (d == 0) {
		fprintf(stderr, "** divide by zero error\n");
		throw std::invalid_argument("divide by zero");
	}

	/* because no mpz_fdiv_q_si exists in the library we need to convert the divisor to a bigint */
	mpz_init_set_si(mpd, d);       // mpd = d
	mpz_init_set(nsave, ZT(n));        // nsave = n
	mpz_inits(mpr, qq, NULL);
	mpz_fdiv_qr(qq, mpr, ZT(n), mpd);  // qq = n/d, mpr=remainder
	remainder = mpz_get_si(mpr);  // magnitude of remainder is < d, so it can't overflow
	if (remainder != 0) {
		//gmp_printf("**temp DivLargeNumberRem %Zd/%lld = %Zd  rem=%lld \n", nsave, d, qq, remainder);
	}
	*q = qq;
	mpz_clears(mpr, mpd, nsave, qq, NULL);
	return remainder;
}

/* q = n/d    return value is remainder. Uses 'truncate' division
i.e. rounds q towards zero, and remainder will have the same sign as n.
the 3 different types of division only give different results when
the remainder is non-zero and n or d (or both) are negative
NB. The choice of floor or truncation type division is sometimes
critical, and the correct choice is not obvious. */
long long tDivLargeNumber(const bigint &n, const bigint &d, bigint &q) {
	long long remainder;
	mpz_t mpr, nsave, qq;

	if (d == 0) {
		fprintf(stderr, "** divide by zero error\n");
		throw std::invalid_argument("divide by zero");
	}

	/* because no mpz_tdiv_q_si exists in the library we need to convert the divisor to a bigint */
	mpz_init_set(nsave, ZT(n));        // nsave = n
	mpz_inits(mpr, qq, NULL);
	mpz_tdiv_qr(qq, mpr, ZT(n), ZT(d));  // q = n/d, mpr=remainder
	remainder = mpz_get_si(mpr);  // magnitude of remainder is < d, so it can't overflow
	if (remainder != 0) {
		//gmp_printf("**temp tDivLargeNumber %Zd/%lld = %Zd  rem=%lld \n", nsave, d, qq, remainder);
	}
	q = qq;
	mpz_clears(mpr, nsave, qq, NULL);
	return remainder;
}

/* return quotient (= n/d) as a normal integer. uses floor division.
i.e. rounds q down towards -infinity, and r will have the same sign as d.
(MulPrToLong will throw an exception if the quotient > 64 bits) */
long long DivLargeNumberLL(const bigint &n, const bigint &d) {
	mpz_t q, r;
	long long llquot;

	if (d == 0) {
		fprintf(stderr, "** divide by zero error\n");
		throw std::invalid_argument("divide by zero");
	}

	mpz_inits(q, r, NULL);
	mpz_fdiv_qr(q, r, ZT(n), ZT(d));    // note use of floor division
	llquot = MulPrToLong(q);             // quotient convert to normal integer
	mpz_clears(q, r, NULL);		           // avoid memory leakage
	return llquot;
}

/*  Dest = CPrev*Bi_Prev + CAct*Bi_Act */
//void MultAddLargeNumbers(long long CPrev, const bigint Bi_Prev,
//	long long CAct, const bigint Bi_Act, bigint *Bi_Dest) {
//	
//	*Bi_Dest = CPrev*Bi_Prev + CAct*Bi_Act;
//}

/* calculate K and L. values are returned in  K and L
It returns true if K and L are valid.
called from ShowSols */
bool CalculateKandL(const long long A, const long long B, const long long C, 
	const long long D, const long long E, const bigint &Bi_R, const bigint &Bi_s, 
	bigint &L, bigint &K) {
	bigint Dp_H1, Dp_H2, Dp_K1;
#ifdef temp
	std::cout << "**temp CalculateKandL "
	<< " Bi_R=" << Bi_R << " Bi_s=" << Bi_s;
	std::cout << " A=" << A << " B=" << B << " C=" << C << " D=" << D
	<< " E=" << E << "\n";
#endif
	K = 2 * Bi_R + B*Bi_s;   // 
	Dp_H1 = K - 2;           /* Kd = 2r + B*s - 2 */
	Dp_K1 = Bi_R - 1;         /* r - 1 */
	Dp_K1 = -B*Dp_K1 - 2 * A*C* Bi_s;
	K = (bigint)C*D*Dp_H1 + E* Dp_K1;
	if (tDivLargeNumber(K, 4 * A*C - B*B, K) != 0) {  /* K */
#ifdef temp
		std::cout <<"**temp CalculateKandL: K not integer; return false\n";
#endif
		return false;               /* K not integer */
	}
	L = D*Dp_K1 + (bigint)A*E* Dp_H1;
#ifdef temp
	std::cout << "**temp CalculateKandL " << " L=" << L
	<< " Divisor = " <<   4 * A*C - B*B << "\n";
#endif
	if (tDivLargeNumber(L, 4 * A*C - B*B, L) != 0) {
#ifdef temp
		std::cout << "**temp CalculateKandL: L not integer; return false\n";
#endif
		return false;               /* L not integer */
	}
	L += D* Bi_s;
#ifdef temp
	std::cout << "**temp CalculateKandL: L=" << L << "\n";
#endif
	return true;
}

/* print large number, with - sign if negative. */
void ShowLargeNumber(const bigint &Bi_Nbr) {
	std::string nbrOutput = "";
	char* buffer = NULL;
	size_t msglen, index = 0;

	// convert to null-terminated ascii string, base 10, with sign if -ve
	buffer = mpz_get_str(NULL, 10, ZT(Bi_Nbr));
	msglen = strnlen(buffer, 500000);  // arbitrary limit of 500000
	if (buffer[0] == '-') {  // if number is -ve
		nbrOutput = "-";
		index = 1;
	}
	for (; index < msglen; index++) {
		if ((msglen - index) % digitsInGroup == 0) {
			// divide digits into groups, if there are more than 6 digits
			if ((index > 0 && buffer[0] != '-') || (index > 1)) {
				nbrOutput += " ";  // put space after group of digits 
				if ((msglen - index) > 100) {
					nbrOutput += "... ";   // if more than 100 digits left, don't print them all
					while (msglen > (digitsInGroup + index))
						index += digitsInGroup;
				}
			}
		}
		nbrOutput += buffer[index];
	}
	if (buffer[0] == '-')   // if number is -ve
		msglen--;
	free(buffer);		// avoid memory leakage
	std::cout << nbrOutput;
	if (msglen > 6)
		std::cout << " (" << msglen << " digits)";
}


/* called from ShowRecursionRoot. 
Uses global variables: A, B, C, D, E
type = hyperbolic_homog (homogenous) or hyperbolic_gen (general hyperbolic)
return 1 if K or L not integers, otherwise zero
If K and L are integers, print values for P, Q, K, R, S, L*/
int ShowSols(equation_class type, const bigint &m, const bigint &n) {
	bigint K, L, Mp_H1;
#ifdef temp
	std::cout << "**temp ShowSols  m=" << m << " n=" << n << " g_C=" << g_C << "\n";
#endif
	if (type == hyperbolic_gen) {
		if (!CalculateKandL(g_A, g_B, g_C, g_D, g_E, m, n, L, K)) {     /* if K or L not integers */
			return 1;                      /* bye */
		}
	}
	if (teach) {
		printf("m = ");
		ShowLargeNumber(m);
		printf("\nn = ");
		ShowLargeNumber(n);
		printf("\nUsing the formulas: P = m\n");
		printf("Q = -Cn \n");
		if (type == hyperbolic_gen) {
			printf("K =CD(P+S-2) + E(B-Bm-2ACn)4AC-B^2\n");
		}
		printf("R = An \nS = m + Bn\n");
		if (type == hyperbolic_gen) {
			printf("L =D(B-Bm-2ACn) + AE(P+S-2)4AC - B^2 + Dn\n");
		}
		printf("we obtain:\n");
	}
	printf("P = ");
	ShowLargeNumber(m);
	printf("\nQ = ");
	Mp_H1 = -g_C* n;   
	ShowLargeNumber(Mp_H1);
	if (type == hyperbolic_gen) {
		printf("\nK = "); 	ShowLargeNumber(K);
	}
	Mp_H1 = g_A* n;  
	printf("\nR = "); ShowLargeNumber(Mp_H1);

	Mp_H1 = m + g_B*n;
	printf("\nS = "); ShowLargeNumber(Mp_H1);

	if (type == hyperbolic_gen) {
		printf("\nL = ");  ShowLargeNumber(L);
	}
	putchar('\n');
	return 0;
}

/* called from ShowRecursion
type = hyperbolic_homog (homogenous) or hyperbolic_gen (general hyperbolic)
Uses global variables Bi_H2, Bi_K2, A, B, C, D, E */
void ShowRecursionRoot(equation_class type) {
	char t;   // 1 if K and L not integers, otherwise 0
	bigint Mp_tmp1, Mp_tmp2, Mp_H1, Mp_K1;
	t = ShowSols(type, Bi_H2, Bi_K2);
	if (type == hyperbolic_gen) {
		Bi_H2 = -Bi_H2; //ChangeSign(&Bi_R);
		Bi_K2 = -Bi_K2; //ChangeSign(&Bi_s);
		if (t == 0) putchar('\n');    /* separate 2 sets of values for  P, Q, K, R, S, L
				If there is only 1 set we just get an unneeded blank line.*/
				//if (t == 1 && ShowSols(type) == 1) {
				/* changed this so that ShowSols is called again even if t is 1.
				ShowSols has an important side-effect; it prints values for P, Q, K, R, S, L
				if it can find them. Sometimes each call produces a different but valid
				set of values */
		if (ShowSols(type, Bi_H2, Bi_K2) == 1 && t == 1) {
			/* if we get to here we haven't yet got any values for P, Q, K, R, S, L.
			The 3rd way below tends to produce very large numbers. */
			if (teach) {
				printf("m = ");
				ShowLargeNumber(Bi_H2);
				printf("\nn = ");
				ShowLargeNumber(Bi_K2);
				std::cout << "\nUsing the formulas: P = m" << sq << " - ACn" << sq <<
					" \nQ = -Cn(2m + Bn) \n K = -n(Em + CDn)\n" <<
					"R = An(2m + Bn)\nS = m" << sq << " + 2Bmn + (B" + sq << " - AC)n" << sq <<
					"\nL = n(Dm + (BD-AE)n) \nwe obtain:\n";
			}
			Mp_H1 = Bi_H2 * Bi_H2;    /* m^2 */
			Mp_K1 = Bi_H2 * Bi_K2;    /* mn */
			Mp_tmp1 = Bi_K2 * Bi_K2;    /* n^2 */

			Mp_tmp2 = Mp_H1 - g_A*g_C* Mp_tmp1;
			printf("P = ");  ShowLargeNumber(Mp_tmp2);

			Mp_tmp2 = (-2 * g_C* Mp_K1) + (-g_B*g_C* Mp_tmp1);
			printf("\nQ = ");   ShowLargeNumber(Mp_tmp2);

			Mp_tmp2 = (-g_E* Mp_K1) + (-g_C*g_D* Mp_tmp1);
			printf("\nK = "); 	ShowLargeNumber(Mp_tmp2);

			Mp_tmp2 = 2 * g_A* Mp_K1 + g_A*g_B* Mp_tmp1;
			printf("\nR = ");  ShowLargeNumber(Mp_tmp2);

			Mp_tmp2 = (g_B*g_B - g_A*g_C)* Mp_tmp1 + 2 * g_B* Mp_K1;
			Mp_tmp2 += Mp_H1; 
			printf("\nS = ");  ShowLargeNumber(Mp_tmp2);

			Mp_tmp2 = g_D* Mp_K1 + (g_B*g_D - g_A*g_E)* Mp_tmp1;
			printf("\nL = ");   ShowLargeNumber(Mp_tmp2);   putchar('\n');
		}
	}
}

/* called from classify. Returns true if equation has no solutions*/
bool CheckMod(long long R, long long S, long long X2, long long X1, long long X0) {
	long long Y2 = gcd(R, S);
	size_t i = 0;
	long long fact = primeList[i];
	int indH = 1;
	long long D1 = abs(Y2);
	long long factors[64];    // holds factors of Y2
	if (teach) {
		if (Y2 != 1) {
			printf("We try to check the equation modulo the prime divisors of ");
			if (R != 0 && S != 0) {
				//w("gcd(" + numToStr(abs(R)) + "," + numToStr(abs(S)) + ") = " + numToStr(Y2) + ".");
				printf("gcd(%lld, %lld) = %lld.", R, abs(S), Y2);
			}
			else {
				//w(numToStr(Y2) + ".\n");
				printf("%lld. \n", Y2);
			}
		}
	}
	/* store factors of D1 in array factors*/
	while (D1 >= fact*fact) {
		int T = 0;
		while (D1%fact == 0) {
			T++;
			if (T == 1) {
				factors[indH++] = fact;  // only save factor once
			}
			D1 /= fact;      // remove factor
		}
		i++;
		if (i >= prime_list_count) {
			std::cerr << "Cannot factorise " << D1 << ". Try increasing max_prime \n";
			throw std::invalid_argument("Number too large to factorise");
		}
		fact = primeList[i];
	}

	if (D1>1) {
		factors[indH++] = D1;
	}

	/* now have list of factors. indH is the number of factors*/
	for (int T = 1; T<indH; T++) {
		if (factors[T]>1) {
			bigint Z = (bigint)((X1*X1 - 4 * X0*X2) % factors[T]) + factors[T];
			long long N = (factors[T] - 1) / 2;
			long long Y = 1;
			while (N != 0) {
				if (N % 2 != 0) {
					Y = MulPrToLong((Y*Z) % factors[T]);
				}
				N /= 2;
				Z = (Z*Z) % factors[T];
			}
			if (Y>1) {
				if (teach) {
					printf("There are no solutions modulo %lld,", factors[T]);
					NoSol();
				}
				return true;   // equation has no solutions
			}
		}
	}
	return false;  // equation appears to have solution(s)
}

/* convert extended precision Nbr to string */
std::string numToStr(const bigint &Dp_Nbr) {
	std::string nbrOutput;
	char* buffer = NULL;
	// convert to null-terminated ascii string, base 10, with sign if -ve
	buffer = mpz_get_str(NULL, 10, ZT(Dp_Nbr));
	nbrOutput = buffer;  // copy from C-style string to STL-type string
	free(buffer);		// avoid memory leakage
	return nbrOutput;
}

/* convert integer to a string */
std::string numToStr(long long num) {
	char numbr[21];
	sprintf_s(numbr, sizeof(numbr), "%lld", num);
	g_s_number = numbr;  // copy to C++ string 
	return g_s_number;
}

/* print ax^2 +bxy +cy2  (with apppropriate signs)*/
void ShowBigEq(const bigint &Dp_A, const bigint &Dp_B, const bigint &Dp_C, 
	const std::string &x, const std::string &y) {
	if (Dp_A != 0) {
		if (Dp_A !=1)
			gmp_printf("%Zd%s%s ", ZT(Dp_A), x.c_str(), sq.c_str());
		else
			gmp_printf("%s%s ", x.c_str(), sq.c_str());
	}
	if (Dp_B !=0) {
		gmp_printf("%+Zd%s%s ", ZT(Dp_B), x.c_str(), y.c_str());
	}
	if (Dp_C != 0) {
		gmp_printf("%+Zd", ZT(Dp_C));
	}
	if (y.length() == 0)
		std::cout << " ";
	else
		std::cout << " " << y << sq << " ";
}

/* convert biginteger to normal. Checks for overflow */
long long MulPrToLong(const bigint &x) {
	long long rv;
	// note: do not use mpz_fits_slong_p because it checks whether x fits a 32 bit integer, rather than a 64 bit integer.
	if (x >= LLONG_MIN && x <= LLONG_MAX) { // is x value OK for normal integer?
		rv = mpz_get_si(ZT(x)); // convert to normal integer
		return rv;
	}
	else
		throw std::range_error ("big number cannot be converted to 64-bit integer");
	return 0;
}


/* prepare for solving by continued fractions.
return values in global variables pDisc, SqrtDisc, Bi_NUM, Bi_DEN */
void GetRoot(const bigint &BiA, const bigint &BiB, const bigint &BiC, bigint *pDisc) {
	long long A, B, T, Z;
	bigint  BiP, BiM, BiZ, BiG, BiK, BiDisc;
	bool NUMis0;

	BiDisc = BiB*BiB - 4 * BiA*BiC;
   	*pDisc = BiDisc;         // Discriminant = B^2 -4AC (this is used later by caller)
	Bi_NUM = -BiB;                    
	NUMis0 = (Bi_NUM == 0);					/* check whether NUM == 0 */
	Bi_DEN = BiA*2;  
	SqrtDisc = llSqrt(BiDisc);        // copy sqrt(Disc) to global var. Assume it won't overflow

									  /* temporary */
	/*std::cout << "**temp** getroot: BiA=" << BiA;
	std::cout << " BiB=" << BiB;
	std::cout << " BiC=" << BiC;
	std::cout << " Bi_NUM=" << Bi_NUM;
	std::cout << " Bi_DEN=" << Bi_DEN;
	std::cout << " SqrtDisc =" << SqrtDisc << "  BiDisc =" << BiDisc << "\n";*/
	/* end temporary */

	if (teach) {
		printf("We have to find the continued fraction expansion of the roots of \n");
		ShowBigEq(BiA, BiB, BiC, "t", "");
		printf("= 0, that is, ");
		if (!NUMis0) {
			std::cout << "(";
		}
		std::cout << "Sqrt(" << *pDisc << ") ";
		if (!NUMis0) {  // print NUM if it's not zero
			gmp_printf("%+Zd) ", ZT(Bi_NUM));  // print num with sign + or -
		}

		if (Bi_DEN< 0) {
			std::cout << " / (" << Bi_DEN << ")";  // print /(-den) if negative
		}
		else {
			std::cout << " / " << Bi_DEN;          // print /den if positive
		}
		putchar('\n');

		gcd(Bi_NUM, Bi_DEN, BiG);     // BiG = gcd(NUM,DEN) = gcd (B,2*A)
		BiZ = BiG * BiG;			   // BiZ = gcd(NUM,DEN)^2
		gcd(BiZ, BiDisc, BiG);	       // BiG = gcd(NUM^2, DEN^2, BiDisc)
		A = MulPrToLong(BiG);        // convert gcd. assume that gcd fit into 64 bits!! 
		B = 1;
		T = 3;

		// remove any repeated factors from A (=BiG), put sqrt(product of repeated factors) in B
		adjustGandK(A, B);
		/* B is the product of the factors removed from A (which is the GCD of NUM^2, DEN^2 , BiDisc)*/
		if (B != 1) {
			DivLargeNumberRem(Bi_NUM, B, &Bi_NUM);   // NUM /= B (floor division)
			DivLargeNumberRem(Bi_DEN, B, &Bi_DEN);   // DEN /= B (floor division)
			bool DENis1 = (Bi_DEN == 1);  /* check whether DEN == 1*/
			printf("Simplifying, ");
			if (!DENis1 && !NUMis0) {
				printf("(");
			}
			std::cout << "Sqrt(" << *pDisc / (B*B) << ")";
			if (!NUMis0) {
				gmp_printf("%+Zd", ZT(Bi_NUM));   // print NUM with sign + or - 
			}
			if (!DENis1 && !NUMis0) {
				printf(")");
			}
			/* print "/DEN"  unless DEN = 1*/
			if (!DENis1) {
				printf(" / ");
				if (Bi_DEN < 0) {
					std::cout << "(" << Bi_DEN << ")";  // print (-DEN)
				}
				else {
					std::cout << Bi_DEN;            // print DEN without brackets round it
				}
			}
			putchar('\n');

			Bi_NUM = -BiB;        // restore values of Bi_Num                 
			Bi_DEN = BiA * 2;       // and Bi_DEN;
		}

		printf("The continued fraction expansion is: \n");
		BiP = Bi_DEN;     

		/* if DEN >= 0 then K = SqrtDisc else K = SqrtDisc+1 */
		BiK = SqrtDisc + ((Bi_DEN < 0) ? 1 : 0); 
		BiK += Bi_NUM; 
		Z = DivLargeNumberLL(BiK, Bi_DEN);      // Z = K/DEN
		BiM = Z;             // M = K/DEN
		BiK = BiM * Bi_DEN;    // K = M*DEN
		BiM = BiK - Bi_NUM;    // M = K-NUM
												 
		printf("%lld", Z);
		std::string sep = "+ //";
		int cont = -1;
		{
			bigint B, P, C, M;
			B = P = C = M = -1;

			while (cont < 0 || B != P || C != M) {
				std::cout << sep;
				sep = ", ";

				if (cont < 0 && (BiP > 0) && (BiP < SqrtDisc + BiM)
					/*  is P > 0 and P < SqrtDisc+M ?*/
					&& (BiM > 0) && (BiM <= SqrtDisc))
					/*  is M > 0 and M <= SqrtDisc ?*/ {

					putchar('\n');
					B = P = BiP;
					C = M = BiM;
					cont = 0;
				}

				if (cont >= 0) {
					P = (*pDisc - M*M) / P;  /* both numerator and denominator are positive */
					Z = MulPrToLong((SqrtDisc + M) / P);
					M = Z*P - M;
					cont++;
				}
				else {
					BiG = BiDisc - BiM * BiM;        // G = Disc-M^2
					DivLargeNumber(BiG, BiP, &BiK);  // K=G/P
					BiP = BiK;
					BiZ = SqrtDisc + ((BiK < 0) ? 1 : 0);
					Z = DivLargeNumberLL(BiZ + BiM, BiK);
					BiM = Z * BiK - BiM;
				}
				printf("%lld", Z);
			}
		}
		printf("//\nwhere the periodic part is on 2nd line");
		if (cont>1) {
			printf(" (the period has %d coefficients)", cont);
		}
		putchar('\n');
	}
}

/* return gcd (P,Q,R), unless gcd is not a factor of S.
In that case return 0 */
bigint DivideGcd(bigint P, bigint Q, bigint R, bigint S,
	const std::string &x1, const std::string &y) {
	int t;
	bigint N;
	gcd(Q, R, N);
	gcd(N, P, N);
	if (N != 1) {
		if (teach) {
			printf("We must get the gcd of all terms except the constant: \n");
			//w("gcd(" + numToStr(P));
			std::cout << "gcd(" << P;

			if (Q != 0) {
				//w(", " + numToStr(Q));
				std::cout << ", " << Q;
			}
			if (R != 0) {
				//w(", " + numToStr(R));
				std::cout << ", " << R;
			}
			//w(") = " + numToStr(N) + "\n");
			std::cout << ") = " << N << "\n";
		}
		if (S%N != 0) {
			NoGcd(S);
			return 0;
		}
		if (teach) {
			std::cout << divgcd;    // show new values after division by gcd
			P /= N;
			Q /= N;
			R /= N;
			S /= N;
			ShowEq(P, 0, 0, Q, R, S, x1, y);
			printf(" = 0 \n");
		}
	}
	if (teach) {
		if (abs(R) != 1) {
			printf("This means that ");
			ShowEq(P, 0, 0, Q, 0, S, x1, y);
			//w(" should be a multiple of " + numToStr(abs(R)) + "\n");
			std::cout << " should be a multiple of " << abs(R) << "\n";
			std::cout << "To determine this, we should try all values of " << x1 << " from 0 to " <<
				(abs(R) - 1) << " to check if the condition holds.\n";
			t = 0;
			for (long long u = 0; u<abs(R); u++) {
				if ((P*u*u + Q*u + S) % R == 0) {
					if (t != 0) {
						std::cout << ", " << u;
					}
					else {
						std::cout << "The values of " << x1 << " (mod " <<
							abs(R) << ") are: " << u;
						t = 1;
					}
				}
			}
			if (t == 0) {
				std::cout << "The modular equation is not satisfied by any " << x1;
				NoSol();
				return 0;
			}
			putchar('\n');
		}
		else {
			t = 1;
		}
	}
	return N;
}

/* compare 2 Bi numbers*/
//int Compare(const bigint Bi_array, const bigint Bi_K1) {
//	if (Bi_array == Bi_K1) return 0;
//	if (Bi_array > Bi_K1) return 1;
//	else return -1;
//}

//int Comparev(const void * Bi_array, const mpz_t &Bi_K1) {
//	mpz_t tempZ;
//	memcpy(tempZ, Bi_array, sizeof(mpz_t));
//	return mpz_cmp(tempZ, Bi_K1);
//}

/* class used to store solutions to equation. Operator < is overloaded so that 
STL functions can handle this class without defining a special compare function. */
class solution {
public:
	bigint x, y;
	bool operator < (const solution &p) const {
		if (x == p.x)
			return y < p.y;
		else
			return x < p.x;
	}
	bool operator > (const solution &p) const {
		if (x == p.x)
			return y > p.y;
		else
			return x > p.x;
	}
	bool operator == (const solution &p) const {
		return (x == p.x) && (y == p.y);
	}
};

/* set used to store solutions. Arguably overkill for this use, but it has the
 convenient property that duplicate solutions are detected automatically */
std::set <solution> SortedSols;

/* Insert new number into sorted solutions set.  Duplicate solutions are not *
*  saved. Called from ShowLargeXY                                            */
void InsertNewSolution2(const bigint &Bi_x, const bigint &Bi_y) {
	solution newval;

	newval.x = Bi_x;
	newval.y = Bi_y;
	/* note that insert keeps the solutions in ascending order and does not 
	store duplicates. */
	if (SortedSols.insert(newval).second)
		NbrSols++;  // increment counter unless the solution is a duplicate.
}


/* convert num to digits, in brackets if -ve */
std::string par(long long num) {
	if (num<0) {
		return "(" + numToStr(num) + ")";
	}
	return "" + numToStr(num);
}
std::string par(const bigint &num) {
	if (num<0) {
		return "(" + numToStr(num) + ")";
	}
	return "" + numToStr(num);
}

/* convert num to digits, in brackets if -ve, but return empty string if num == 1 */
std::string par1(long long num) {
	return (num == 1) ? "" : par(num);
}
std::string par1(const bigint &num) {
	return (num == 1) ? "" : par(num);
}

/* this function either calls InsertNewSolution to save the solution
or outputs the the solution: "x=<val>  y=<val>"
In some cases the -ve value of x and/or y is also stored/printed */
void ShowLargeXY(const std::string &x, const std::string &y, const bigint &Bi_X, const bigint &Bi_Y,
	bool sol, const std::string &eqX, const std::string &eqY) {
	/* if sol == true use both +ve and -ve of solution */

	if ((x == "X") && (y == "Y")) {

		/*std::cout << "**temp ShowLargeXY X="; 	ShowLargeNumber(Bi_X);
		std::cout << " Y=";  ShowLargeNumber(Bi_Y);
		std::cout << "  sol=" << sol << "\n";*/

		if (!teach && !allSolsFound) {
			also = true;
			/* store solutions in sorted solutions vector and exit*/
			if (sol && (Bi_Y == 0)) {
				InsertNewSolution2(Bi_X, Bi_Y);
				InsertNewSolution2(-Bi_X, Bi_Y);
			}
			else {
				if (sol && (Bi_Y < 0)) {
					InsertNewSolution2(-Bi_X, -Bi_Y);
				}
				else {
					/* sol is false and/or Bi_Y > 0 */
					InsertNewSolution2(Bi_X, Bi_Y);
				}
			}
			//if (!allSolsFound)
			std::cout << NbrSols << " solutions \n";
			return;
		}
	}

	if (teach && y == "Y") {
		putchar('\n');
	}
	if (y == "Y") {
		showAlso();
	}
	std::cout << x << " = ";
	if (teach) {
		std::cout << eqX;
	}
	ShowLargeNumber(Bi_X);
	std::cout << "    " << y << " = ";
	if (teach) {
		std::cout << eqY;
	}
	ShowLargeNumber(Bi_Y);
	putchar('\n');

	if (y == "Y" && sol) {
		/* show -ve of solution as well */
		std::cout << "\n  and also:  " << x << "0 = ";
		if (teach && eqX.length() >0) {
			std::cout << "-" << eqX;
		}
		ShowLargeNumber(-Bi_X);
		std::cout << "\n" << y << "0 = ";
		if (teach && eqY.length() >0) {
			std::cout << "-" << eqY;
		}
		ShowLargeNumber(-Bi_Y);
		putchar('\n');
	}
}


/* show all solutions stored in sorted solutions set, stored by
InsertNewSolution function.
Note that solutions are deleted after they are printed, making the set
ready for reuse, but this function can only be called once. */
void  ShowAllLargeSolutions2() {
	for (auto sol : SortedSols) {
		std::cout << "X= ";     ShowLargeNumber(sol.x);
		std::cout << "  \tY= "; ShowLargeNumber(sol.y);  // tab makes output a bit more tidy
		std::cout << "\n";
	}
	SortedSols.clear();
}


/* quick print of sortedSols - for testing */
void listLargeSolutions() {
	for (auto sol : SortedSols) {
		std::cout << "X= ";     ShowLargeNumber(sol.x);
		std::cout << "  \tY= "; ShowLargeNumber(sol.y);  // tab makes output a bit more tidy
		std::cout << "\n";
	}
}

/* print x = Xi +X1*t  y = Yi +Y1*t 
output is 'prettied up' by:
1.	if Xi =0 & X1  = 0 just print 'x=0'
2.	if Xi =0 & X1 NE 0 don't print 'xi +' i.e just print x = X1*t 
3. if value of X1 is +/- 1 just print x = Xi +/-t without the digit 1.
4. same rules apply for y */
void PrintLinear(long long Xi, long long Xl, long long Yi, long long Yl, const std::string &va) {
	if (va == "t") {         // actually, va always = t
		showAlso();
	}

	if (ExchXY) {
		long long T = Xi; Xi = Yi; Yi = T;  /* swap Xi and Yi */
		long long T2 = Xl; Xl = Yl; Yl = T2;            /* swap X1 and Y1 */
	}

	printf("x = ");
	if (Xi == 0 && Xl == 0) {
		printf("0");  // print x = 0
	}
	else {
		if (Xi != 0) {
			std::cout << Xi;
		}
		if (Xl < 0) {
			printf(" - ");
		}
		if (Xl > 0 && Xi != 0) {
			printf(" + ");
		}
		if (Xl != 0) {
			if (abs(Xl) == 1) {
				std::cout << " " << va;
			}
			else {
				std::cout << abs(Xl) << va;
			}
		}
	}

	printf("   y = ");
	if (Yi == 0 && Yl == 0) {
		printf("0");
	}
	else {
		if (Yi != 0) {
			std::cout << Yi;
		}
		if (Yl < 0) {
			printf(" - ");
		}
		if (Yl > 0 && Yi != 0) {
			printf(" + ");
		}
		if (Yl != 0) {
			if (abs(Yl) == 1) {
				std::cout << " " << va;
			}
			else {
				std::cout << abs(Yl) << va;
			}
		}
	}
	putchar('\n');
	return;
}
void PrintLinear(bigint Xi, bigint Xl, bigint Yi, bigint Yl, const std::string &va) {
	if (va == "t") {         // actually, va always = t
		showAlso();
	}

	if (ExchXY) {
		bigint T = Xi; Xi = Yi; Yi = T;  /* swap Xi and Yi */
		bigint T2 = Xl; Xl = Yl; Yl = T2;            /* swap X1 and Y1 */
	}

	printf("x = ");
	if (Xi == 0 && Xl == 0) {
		printf("0");  // print x = 0
	}
	else {
		if (Xi != 0) {
			std::cout << Xi;
		}
		if (Xl < 0) {
			printf(" - ");
		}
		if (Xl > 0 && Xi != 0) {
			printf(" + ");
		}
		if (Xl != 0) {
			if (abs(Xl) == 1) {
				std::cout << " " << va;
			}
			else {
				std::cout << abs(Xl) << va;
			}
		}
	}

	printf("   y = ");
	if (Yi == 0 && Yl == 0) {
		printf("0");
	}
	else {
		if (Yi != 0) {
			std::cout << Yi;
		}
		if (Yl < 0) {
			printf(" - ");
		}
		if (Yl > 0 && Yi != 0) {
			printf(" + ");
		}
		if (Yl != 0) {
			if (abs(Yl) == 1) {
				std::cout << " " << va;
			}
			else {
				std::cout << abs(Yl) << va;
			}
		}
	}
	putchar('\n');
	return;
}

/* called from solveEquation, SolveParabolic & SolveDiscIsSq.
return 0 if solution found, 1 if no solutions exist, 2 if there are an infinite
number of solutions.
The equation to be solved is of the form Dx + Ey + F =0 */
int Linear(long long D, bigint E, bigint F) {
	long long Tx;
	long long Yl;
	bigint Xi, Xl, Yi;
	int t;

	if (teach) {
		printf("This is a linear equation ");
		ShowLin(D, E, F, "x", "y");
		printf(" = 0 \n");
	}
	if (D == 0) {
		if (E == 0) {
			if (F != 0) {
				return 1;             // No solutions
			}
			else {
				printf("x, y: any integer\n");  /* 0X + 0Y +0 = 0 */
				return 2;             // Infinite number of solutions
			}
		}
		if (F%E != 0) {
			return 1;               // No solutions
		}
		else {
			/*Xi = 0;
			Xl = 1;
			Yi = -F / E;
			Yl = 0;*/
			PrintLinear(0, 1, -F/E, 0, "t");
			return 0;               // Solution found
		}
	}
	if (E == 0) {
		if (F%D != 0) {
			return 1;               // No solutions
		}
		else {
			/*Xi = -F / D;
			Xl = 0;
			Yi = 0;
			Yl = 1;*/
			PrintLinear(-F / D, 0, 0, 1, "t");
			return 0;               // Solution found
		}
	}

	bigint Q;  // gcd of D & E
	gcd(D, E, Q);   // since D is < 64 bits, Q also < 64 bits
	if (Q != 1 && Q != -1) {
		if (teach) {
			std::cout << "To solve it, we first find the gcd of the linear coefficients, that is: gcd(" <<
				D << ", " << E << ") = " << Q << ".\n";
		}
		if (F%Q != 0) {
			NoGcd(F);
			return 1;               // No solutions
		}
		D = MulPrToLong(D / Q);   // no overflow problem here; magnitude of quotient < D
		E = E / Q;
		F = F / Q;   // divide by gcd
	}

	if (teach) {
		if (Q != 1) {
			std::cout << divgcd;    // show new values after division by gcd
			ShowLin(D, E, F, "x", "y");
			printf(" = 0 \n");
		}
		printf("Now we must apply the Generalized Euclidean algorithm: \n");
	}

	bigint U1 = 1;
	long long U2 = 0;
	bigint U3 = D;
	bigint V1 = 0;
	long long V2 = 1;
	bigint V3 = E;
	t = 1;
	while (V3 != 0) {
		if (teach) {
			std::cout << "Step " << t << ": " << par(U1) << "*" << par(D) << " + " << par(U2) << "*"
				<< par(E) << " = " << U3 << "\n";
		}
		long long q = DivLargeNumberLL(U3, V3);
		bigint T1 = U1 - q*V1;
		long long T2 = U2 - q*V2;
		bigint T3 = U3 - q*V3;
		U1 = V1;  
		V1 = T1; 
		U2 = V2;
		V2 = T2; 
		U3 = V3;
		V3 = T3;
		t++;
	}
	Xi = -U1*F / U3; 
	Xl = E; 
	Yi = -U2*F / U3; 
	Yl = -D;
	if (teach) {
		std::cout << "Step " << t << ": " << par(U1) << "*" << par(D) << " + " << par(U2) << "*" <<
			par(E) << " = " << U3 << "\n";
		std::cout << "\nMultiplying the last equation by " << par(-F / U3) << " we obtain:\n";
		std::cout << par(Xi) << "*" << par(D) << " + " << par(Yi) << "*" << par(E) << " = " << (-F) << "\n";
		std::cout << "Adding and subtracting " << par(D) << "*" << par(E) << "t' we obtain:\n";
		std::cout << "(" << Xi << "+" << par(E) << "t')*" << par(D) << " + (" << Yi << "-" << par(D)
			<< "t')*" << par(E) << " = " << (-F) << "\n";
		printf("So, the solution is given by the set:\n");
		PrintLinear(Xi, Xl, Yi, Yl, "t'");
	}
	V1 = D*D + E*E;
	Tx = DivLargeNumberLL((D*Yi - E*Xi) + V1 / 2, V1);
	if (teach) {
		std::cout << "By making the substitution t = " << Tx << " + t' we finally obtain:\n";
	}
	Xi += E*Tx;
	Yi += -D*Tx;
	if (Xl<0 && Yl<0) {
		Xl = -Xl;
		Yl = -Yl;
	}
	PrintLinear(Xi, Xl, Yi, Yl, "t");
	return 0;   // 0 indicates a solution was found
}

/* uses global variables A, B, C, D, E, F.
returns true if there are no solutions. */
bool Mod(int mod) {
	/* make safeint copies of a, b, c, d, e, f */
	msl::utilities::SafeInt <long long> sA(g_A), sB(g_B), sC(g_C), sD(g_D), sE(g_E), sF(g_F);

	for (msl::utilities::SafeInt <long long> x = 0; x<mod; x++) {
		/* calculate z & t, but check for overflow */
		msl::utilities::SafeInt <long long> z = sA*x*x + sD*x + sF;
		msl::utilities::SafeInt <long long> t = sB*x + sE;
		for (msl::utilities::SafeInt <long long> y = 0; y<mod; y++) {
			if ((z + y*(t + sC*y)) % mod == 0) {
				long long xx = x;
				long long yy = y;
				if (teach) {
					std::cout << "solution found using x= " << xx << ", y=" << yy;
					std::cout << " mod = " << mod << "\n";
				}
				return false;  // there is a solution
			}
		}
	}

	/* no solution found */
	if (teach) {
		std::cout << "No solutions found using mod " << mod << ",";
		NoSol();
	}
	return true;
}

/* called from SolveElliptical */
void ShowElipSol(long long A, long long B, long long D, long long u, const std::string &x, 
	const std::string &x1, const std::string &y,const bigint &w2) {
	if (teach) {
		putchar('\n');
		ShowLin(2 * A, B, D, x, y);
		std::cout << " = " << w2 << "\n";
		ShowLin(2 * A, B, D, x, par(u));
		std::cout << " = " << w2 << "\n";
		std::cout << par1(2 * A) << x << " = " << (w2 - D - B*u) << "\n";
		if ((w2 - D - B*u) % (2 * A) == 0) {

			std::cout << x << " = " << ((w2 - D - B*u) / (2 * A));
		}
		else {
			std::cout << "There is no integer solution for " << x << " \n";
		}
	}
	if ((w2 - D - B*u) % (2 * A) == 0) {
		if (teach) {
			std::cout << "\n-----------------------------------\n";
		}
		ShowXY((w2 - D - B*u) / (2 * A), u);        // display solution
		if (teach) {
			std::cout << "-----------------------------------";
		}
		putchar('\n');
		also = true;
	}
}

/* called from solveEquation
type = hyperbolic_homog (homogenous) or hyperbolic_gen (general hyperbolic)
uses global variables g_A, g_B, g_C, Bi_H2, Bi_K2*/
void ShowRecursion(equation_class type) {
	const std::string t1 = " integer solution of the equation \nm" + sq + " + bmn + acn" + sq + " = ";
	bigint Disc;

	/*std::cout << "**temp ShowRecursion A=" << A << " B=" << B << " C=" << g_C << " Bi_H2=";
	ShowLargeNumber(Bi_H2);
	std::cout << " type=" << type << "\n";*/

	if (type == hyperbolic_gen) {
		std::cout << "X(n+1) = P X(n) + Q Y(n) + K \n";
		std::cout << "Y(n+1) = R X(n) + S Y(n) + L \n";
	}
	else {
		std::cout << "X(n+1) = P X(n) + Q Y(n) \n";
		std::cout << "Y(n+1) = R X(n) + S Y(n) \n";
	}

	if (teach) {
		std::cout << "In order to find the values of P, Q, R, S we have to find first an" << t1;
		ShowEq(1, g_B, g_A*g_C, 0, 0, 0, "m", "n");
		printf(" = 1. \n");
	}

	
	GetRoot(1, g_B, g_A * g_C, &Disc);      /* return values in Disc, SqrtDisc, Bi_NUM, Bi_DEN, 
											used by ContFrac */
	ContFrac(1, 2, 1, 0, 0, 1, g_A, Disc, g_F);

	if (teach) {
		std::cout << "An" << t1;
		ShowEq(1, g_B, g_A*g_C, 0, 0, 0, "m", "n");
		printf(" = 1 is: \n");
	}

	ShowRecursionRoot(type);  // note: sign of Bi_H2, Bi_K2 may be reversed on return

	if (g_B != 0) {
		if (teach) {
			std::cout << "\nAnother" << t1;
			ShowEq(1, g_B, g_A*g_C, 0, 0, 0, "m", "n");
			printf(" = 1 is: \n");
		}
		else {
			printf("\nas well as\n");
		}
		//MultAddLargeNumbers(1, Bi_H2, g_B, Bi_K2, Bi_H2); /* r <- r + Bs */
		Bi_H2 += g_B*Bi_K2;        /* r <- r + Bs */
		Bi_K2 = -Bi_K2;            /* s <- -s */
		ShowRecursionRoot(type);   
	}

}

/* called from SolveSimpleHyperbolic */
void SolByFact(const bigint &R, long long T, long long B, long long D, long long E) {
	if (teach) {
		std::cout << T << " is a factor of " << R << ", so we can set:\n";
		Show1(E, Show(B, "x", 0));
		printf(" =  %lld", T);
		if (E != 0) {
			printf("\nThis means that ");
			Show(B, "x = ", 0);
			printf("%lld", T - E);
		}
		putchar('\n');
		Show1(D, Show(B, "y", 0));
		std::cout << " = " << R << "/" << par(T) << " = " << (R / T);
		if (D != 0) {
			printf("\nThis means that ");
			Show(B, "y = ", 0);
			std::cout << R / T - D;
		}
		putchar('\n');
	}

	if ((T - E) % B == 0 && (R / T - D) % B == 0) {
		if (teach) {
			std::cout << "-----------------------------------\n";
		}
		ShowXY((T - E) / B, (R / T - D) / B);
		if (teach) {
			std::cout << "-----------------------------------\n";
		}
	}
	else {
		if (teach) {
			printf("These equations are not valid in integers.\n");
		}
	}
}

/* The discriminant is a perfect square;
the equation can be represented as the product of 2 linear expressions.
uses global variables A, B, g_D, g_CY0, g_CY1 */
void SolveDiscIsSq(const bigint &BiN0, const std::string &x, const std::string &y, 
	long long SqrtD, const bigint &Bi_g, const bigint &Bi_h, const bigint &BiN1, 
	const std::string &y1, const std::string &x1) {

	long long Yc = llSqrt(Bi_g / Bi_h);
	long long Xc = llSqrt(abs(g_CY1 / Bi_h));
	long long tempL, Fact1;
	bigint Fact2, X1, Y1;

	if (teach) {
		printf("(");
		ShowLin(Yc, Xc, 0, y1, x1);
		printf(") (");
		ShowLin(Yc, -Xc, 0, y1, x1);
		//w(") = " + numToStr(N0 / h) + "\n");
		std::cout << ") = " << (BiN0 / Bi_h) << "\n";
	}

	if (BiN0 == 0) {
		if (teach) {
			printf("\nOne of the parentheses must be zero, so: \n");
		}

		/* If SqrtD is not zero perform loop twice, 2nd time with sign of Xc reversed
		if SqrtD is zero,  for loop is only executed once. */
		for (tempL = (SqrtD == 0 ? 1 : 0); tempL<2; tempL++) {
			if (teach) {
				//w("<LI>");
				ShowLin(Yc, Xc, 0, y1, x1);
				std::cout << " = 0 \n" << par(Yc) << " (";
				ShowLin(0, g_CY1, g_CY0, x, y);
				std::cout << ") + " << par1(Xc) << " (";
				ShowLin(2 * g_A, g_B, g_D, x, y);
				printf(") = 0 \n");
				ShowLin(2 * g_A*Xc, Xc*g_B + g_CY1*Yc, g_D*Xc + g_CY0*Yc, x, y);
				printf(" = 0 \n");
			}
			Linear(2 * g_A*Xc, Xc*g_B + g_CY1*Yc, g_D*Xc + g_CY0*Yc);
			Xc = -Xc;          // reverse sign of Xc
		}

		if (teach) {
			//w("</UL>");
		}
		return;
	}

	/* N0 is not zero, discriminant is a perfect square */
	if (teach) {
		//w("Now we have to find all factors of " + numToStr(N1) + ". \n");
		std::cout << "Now we have to find all factors of " << BiN1 << "\n";
	}
	for (unsigned long long t1 = 1; t1 <= llSqrt(BiN1); t1++) {
		if (BiN1%t1 == 0) {
			Fact1 = t1;       // t1 is a factor of N1, so copy it to Fact1
			Fact2 = BiN0 / Bi_h / t1;
			if (teach) {
				std::cout << "Since " << (Fact1*Fact2) << " is equal to " << Fact1 <<
					" times " << Fact2 << ", we can set:\n";
				ShowLin(Yc, Xc, 0, y1, x1);
				//w(" = " + numToStr(Fact1) + "\n");
				printf(" = %lld \n", Fact1);
				ShowLin(Yc, -Xc, 0, y1, x1);
				//w(" = " + numToStr(Fact2) + "\n");
				std::cout << " = " << Fact2 << "\n";
				if ((Fact1 - Fact2) % (2 * Xc) == 0 && (Fact1 + Fact2) % (2 * Yc) == 0) {
					X1 = (Fact1 - Fact2) / (2 * Xc);
					Y1 = (Fact1 + Fact2) / (2 * Yc);
					std::cout << x1 << " = " << X1 << "\n" << y1
						<< " = " << Y1 << "\n";
					std::cout <<"-----------------------------------------\n";
					ShowX1Y1(X1, Y1, g_A, g_B, g_D, g_CY1, g_CY0);  // display solution
					std::cout << "-----------------------------------------\n";

				}
				else {
					//w("Solving this system we do not obtain integer values for " + x1 + " and " + y1 + ". \n");
					std::cout << "Solving this system we do not obtain integer values for "
						<< x1 << " and " << y1 << ". \n";
				}
			}
			else {
				if ((Fact1 - Fact2) % (2 * Xc) == 0 && (Fact1 + Fact2) % (2 * Yc) == 0) {
					X1 = (Fact1 - Fact2) / (2 * Xc);
					Y1 = (Fact1 + Fact2) / (2 * Yc);
					ShowX1Y1(X1, Y1, g_A, g_B, g_D, g_CY1, g_CY0);  // display solution
				}
			}
		}
	}
	if (teach) {
		putchar('\n');
	}
	return;
}

/* called from solveEquation
uses global variables g_B, g_D, g_E and g_F*/
void SolveSimpleHyperbolic(void) {
	/* simple hyperbolic; A = C = 0; B ≠ 0 */
	long long S, T;
	bigint R;

	R = g_D*g_E - g_B*g_F;
	if (teach) {
		printf("Multiplying by %lld we obtain; \n", g_B);
		ShowEq(0, g_B*g_B, 0, g_D*g_B, g_E*g_B, 0, "x", "y");
		printf(" = %lld \n", (-g_F*g_B));
		printf("Adding %lld to both sides of the equal sign: \n", (g_D*g_E));
		ShowEq(0, g_B*g_B, 0, g_D*g_B, g_E*g_B, g_D*g_E, "x", "y");
		//printf(" = %lld \n", R);
		std::cout << " = " << R << "\n";
		printf("Now the left side can be factored as follows:\n");
		printf("(");
		Show1(g_E, Show(g_B, "x", 0));
		printf(") (");
		Show1(g_D, Show(g_B, "y", 0));
		//printf(") = %lld \n", R);
		std::cout << ") = " << R << "\n";
		if (R != 0) {
			printf("Then ");
			Show1(g_E, Show(g_B, "x", 0));
			//printf(" must be a factor of %lld, so we must find the factors of %lld: \n", R, R);
			std::cout << " must be a factor of " << R << ", so we must find the factors of " << R << ":\n";
		}
		else {
			printf("One of the parentheses must be zero, so: \n");
			Show1(g_E, Show(g_B, " x", 0));
			printf(" = 0");
			if (g_E%g_B == 0) {
				std::cout << " means that x = " << (-g_E / g_B) << " and y could be any integer.\n";
				also = true;
			}
			else {
				printf("This equation cannot be solved in integers.\n");
			}
			Show1(g_D, Show(g_B, " y", 0));
			printf(" = 0");
			if (g_D%g_B == 0) {
				printf(" means that y = %lld and x could be any integer. \n", (-g_D / g_B));
				also = true;
			}
			else {
				printf("This equation cannot be solved in integers. \n");
			}
			return;
		}
	}
	if (R != 0) {
		S = llSqrt(abs(R));  // because R NE 0, S > 0
		for (T = 1; T <= S; T++) {
			if (R%T == 0) {
				SolByFact(R, T, g_B, g_D, g_E);
				SolByFact(R, -T, g_B, g_D, g_E);
				if ((bigint)T*T != abs(R)) {
					SolByFact(R, MulPrToLong(R / T), g_B, g_D, g_E);
					SolByFact(R, MulPrToLong(-R / T), g_B, g_D, g_E);
				}
			}
		} /* end for */
		if (teach) {
			//w("</UL>");
			putchar('\n');
		}
		return;
	}

	if (g_E%g_B == 0) {
		PrintLinear(-g_E/g_B, 0, 0, 1, "t");
	}
	if (g_D%g_B == 0) {
		PrintLinear(0, 1, -g_D/g_B, 0, "t");
	}
	return;
}

/* uses global variables g_A, g_B, g_C, g_D, g_E, g_F*/
void SolveParabolic(std::string x, std::string y, std::string x1) {
	
	const bigint E1 = 4 * g_A*g_E - 2 * g_B*g_D;
	const bigint F1 = 4 * g_A*g_F - g_D*g_D;
	const long long r = gcd(2 * g_A, g_B);
	const long long s = 2 * g_A / r;
	int t;
	std::string t1;
	long long r1 = 1, r2=1, u;
	bigint P, P1, P2, Q, Q1, Q2, R, R1, S, S1, S2,T;
	P = r / 2;
	Q = g_D;
	R = (bigint)(2 * g_A*g_E - g_B*g_D) / r;
	S = (bigint)2 * g_A*g_F / r;

	if (teach) {
		if (s != 1) {
			std::cout << "Multiplying the equation by " << par(s) << ":\n";
			ShowEq(g_A*s, g_B*s, g_C*s, g_D*s, g_E*s, g_F*s, x, y);
			printf(" = 0 \n");
		}
		if (r != 2) {
			std::cout << "Extracting the factor " << (r / 2) << " in the quadratic terms:\n";
			std::cout << par1(r / 2) << "(";
			ShowEq(s*s, 2 * s*g_B / r, 2 * s*g_C / r, 0, 0, 0, x, y);
			printf(")");
			t = Show(g_D*s, " " + x, 1);
			t = Show(g_E*s, " " + y, t);
			Show1(g_F*s, t);
			printf(" = 0 \n");
		}
		if (g_B != 0) {
			std::cout << par1(r / 2) << "(";
			ShowLin(s, g_B / r, 0, x, y);
			std::cout << ")" << sq;
			if (g_D != 0 || g_E != 0 || g_F != 0) {
				printf(" + (");
				ShowLin(g_D*s, g_E*s, g_F*s, x, y);
				printf(")");
			}
			printf(" = 0 \n");
			if (g_D != 0) {
				std::cout << "Adding and subtracting " << par1(g_B*g_D / r) << y << ":\n";
				std::cout << par1(r / 2) << "(";
				ShowLin(s, g_B / r, 0, x, y);
				std::cout << ")" << sq;
				std::cout << " + " << par1(g_D) << " (";
				ShowLin(s, g_B / r, 0, x, y);
				printf(")");
			}
			if (E1 != 0 || g_F != 0) {
				printf(" + (");
				ShowLin(0, R, g_F*s, x, y);
				printf(")");
			}
		}
		else {
			std::cout << par1(r / 2) << "(";
			ShowLin(s, g_B / r, 0, x, y);
			std::cout << ")" << sq;
			if (g_D != 0) {
				std::cout << " + " << par1(g_D) << " (";
				ShowLin(s, g_B / r, 0, x, y);
				printf(")");
			}
			if (E1 != 0 || g_F != 0) {
				printf(" + (");
				ShowLin(0, R, g_F*s, x, y);
				printf(")");
			}
		}
		printf(" = 0 \nNow we perform the substitution:\n");
		std::cout << x1 << " = ";
		ShowLin(s, g_B / r, 0, x, y);
		printf("\nThis gives:\n");
		ShowEq(r / 2, 0, 0, g_D, R, g_F*s, x1, y);
		printf(" = 0 \n");
	}

	if (E1 == 0) {
		if (teach) {
			printf("This can be solved by the standard quadratic equation formula: \n");
			std::cout << "The roots are: " << x1 << " = ";
			if (g_D != 0) {
				std::cout << "(" << par(-g_D) << " +/- sqrt(" << (-F1) << "))";
			}
			else {
				std::cout << " +/- Sqrt(" << (-F1) << ")";
			}
			std::cout << " / " << par(r) << "\n";
		}

		if (F1>0) {
			printf("This quadratic equation has no solution in reals, so it has no solution in integers. \n");
			also = true;
			return;
		}
		/* get 1st solution */
		T = (long long)floor((-g_D + sqrt((double)-F1)) / r + 0.5);
		/* check whether 1st root is an integer */
		if (r*T*T / 2 + g_D*T + 2 * g_A*g_F / r == 0) {
			if (teach) {
				std::cout << "The first root is: " << x1 << " = " << T << "\n";
				ShowLin(s, g_B / r, 0, x, y);
				std::cout << " = " << T << "\n";
			}
			r1 = Linear(s, g_B / r, -T);
		}
		else {
			printf("The first root is not an integer. \n");
		}

		if (F1 == 0) {
			return;   // both roots are equal, skip checking 2nd root
		}
		/* get 2nd solution */
		T = (long long)floor((-g_D - sqrt((double)-F1)) / r + 0.5);

		if (r*T*T / 2 + g_D*T + 2 * g_A*g_F / r == 0) {
			if (teach) {
				std::cout << "The second root is: " << x1 << " = " << T << "\n";
				ShowLin(s, g_B / r, 0, x, y);
				std::cout << " = " << T << "\n";
			}
			r2 = Linear(s, g_B / r, -T);
		}
		else {
			printf("The second root is not an integer. \n");
		}

		if (r1 == 1 && r2 == 1)
			std::cout << msg;     // there are no solutions

		return;
	}

	bigint N = DivideGcd(P, Q, R, S, x1, y);
	if (N == 0) {
		return;   // cannot divide S by gcd(P,Q,R)
	}
	P /= N; Q /= N; R /= N; S /= N;   // divide by their gcd

	int numEq = 1;
	for (u = 0; u<abs(R); u++) {
		T = P*u*u + Q*u + S;
		if (T%R == 0) {
			int numEq2 = 3;
			P1 = g_B*P*R / r; 
			Q1 = abs(R) + g_B*(2 * P*u + Q)*abs(R) / R / r;
			R1 = -s; 
			S1 = g_B*T / R / r + u;
			P2 = -P*R; 
			Q2 = (2 * P*u + Q)*abs(R) / -R; 
			S2 = -T / R;
			if (teach) {
				t1 = ((u == 0) ? "" : numToStr(u) + " + ") + numToStr(abs(R)) + "t";
				std::cout << x1 << " = " << t1 << "   (" << numEq << ")";
				printf("\nReplacing this in the equation shown above:\n");
				std::cout << par1(-R) << y << " = ";
				ShowEq(P, 0, 0, Q, 0, S, "(" + t1 + ")", y);
				std::cout << "\n" << par1(-R) << y << " = ";
				ShowEq(P*R*R, 0, 0, (2 * P*u + Q)*abs(R), 0, T, "t", y);
				std::cout << "\n" << y << " = ";
				ShowEq(P2, 0, 0, Q2, 0, S2, "t", y);
				std::cout << "   (" << (numEq + 1) << ")\nFrom (" << numEq << "): ";
				ShowLin(s, g_B / r, 0, x, y);
				std::cout << " = " << t1 << "\nReplacing (" << (numEq + 1) << ") here:\n";
				std::cout << par1(s) << x << " = ";
				ShowEq(P1, 0, 0, Q1, 0, S1, "t", y);
				std::cout << "   (" << (numEq + 2) << ") \n";
			}
			bigint N1 = DivideGcd(P1, Q1, R1, S1, "t", x);
			if (N1 == 0) {
				continue;   // cannot divide S1 by gcd(P1,Q1,R1)
			}
			P1 /= N1; Q1 /= N1; R1 /= N1; S1 /= N1;
			for (long long u1 = 0; u1<abs(R1); u1++) {
				bigint T1 = P1*u1*u1 + Q1*u1 + S1;
				if (T1%R1 == 0) {
					if (teach && abs(R1) != 1) {
						t1 = ((u1 == 0) ? "" : u1 + " + ") + numToStr(abs(R1)) + "u";
						std::cout << "t = " << t1 << "   (" << (numEq + numEq2) << ")";
						std::cout << "\nReplacing this in (" << (numEq + 2) << "):\n";
						std::cout << par1(-R1) << x << " = ";
						ShowEq(P1, 0, 0, Q1, 0, S1, "(" + t1 + ")", y);
						std::cout << "\n" << par1(-R1) << x << " = ";
						ShowEq(P1*R1*R1, 0, 0, (2 * P1*u1 + Q1)*abs(R1), 0, T1, "u", "");
						std::cout << "\n" << x << " = ";
						ShowEq(-P1*R1, 0, 0, (2 * P1*u1 + Q1)*abs(R1) / -R1, 0, -T1 / R1, "u", y);
						std::cout << "\nFrom (" << (numEq + 1) << ") and (" <<
							(numEq + numEq2) << "):\n" << y << " = ";
						ShowEq(P2, 0, 0, Q2, 0, S2, "(" + t1 + ")", y);
						std::cout << "\n" << y << " = ";
						ShowEq(P2*R1*R1, 0, 0, (2 * P2*u1 + Q2)*abs(R1), 0, P2*u1*u1 + Q2*u1 + S2, "u", "");
						putchar('\n');
						numEq2++;
					}
					if (teach) {
						std::cout << "-----------------------------------------------------\n";
					}
					showAlso();
					printf("x = ");
					if (!ExchXY) {
						ShowEq(-P1*R1, 0, 0, (2 * P1*u1 + Q1)*abs(R1) / -R1, 0, -T1 / R1, "u", "");
						printf("\ny = ");
						ShowEq(P2*R1*R1, 0, 0, (2 * P2*u1 + Q2)*abs(R1), 0, P2*u1*u1 + Q2*u1 + S2, "u", "");
					}
					else {
						ShowEq(P2*R1*R1, 0, 0, (2 * P2*u1 + Q2)*abs(R1), 0, P2*u1*u1 + Q2*u1 + S2, "u", "");
						printf("\ny = ");
						ShowEq(-P1*R1, 0, 0, (2 * P1*u1 + Q1)*abs(R1) / -R1, 0, -T1 / R1, "u", "");
					}
					if (teach) {
						std::cout << "\n-----------------------------------------------------";
					}
					putchar('\n');
				}
			}
			if (teach && abs(R1) != 1) {
				//w("</UL>");
			}
			numEq += numEq2;
		}
	}

	if (teach) {
		//w("</UL>");
	}
	return;   /* end of specific processing for parabolic case*/
}

/* uses global variables g_A, g_B, g_C, g_D, g_E, g_F */
void SolveElliptical(const std::string &x, const std::string &x1, const std::string &y) {
	
	const bigint NegDisc = -Bi_Disc;   // get -ve of discriminant
	const bigint E1 = 4 * (bigint)g_A*g_E - 2 * (bigint)g_B*g_D;
	const bigint F1 = 4 * (bigint)g_A*g_F - (bigint)g_D*g_D;
	const bigint g = gcd(NegDisc, E1/2);
	const bigint CY1 = NegDisc / g;
	const bigint CY0 = E1/2/ g;
	const bigint N0 = CY0*CY0*g - CY1*F1;
	// If N0 is -ve, square root is not a real number; test later whether N0  is negative
	const double sqrtgN0 = sqrt((double)g*(double)N0);  

	int b;
	long long u;
	bigint w1, w2;

	const double R3 = (double(-E1/2) - sqrtgN0) / (double)NegDisc;    // get minimum for x as a real number
	const double R4 = (double(-E1/2) + sqrtgN0) / (double)NegDisc;    // get maximum for x as a real number

	if (R4 > LLONG_MAX || R3 < LLONG_MIN) {
		throw std::range_error("floating-point number outside range of 64-bit integer");
	}
	long long R2 = (long long)floor(R4);          // get maximum for x as an integer
	long long R1 = (long long)ceil(R3);           // get minimum for x as a real number

	if (teach) {
		std::cout << "Since " << x1 << sq << " is always greater than, or equal to zero,\n";
		ShowEq(0, 0, NegDisc, 0, E1, F1, x1, y);
		printf(" must be less than, or equal to zero. \n"
			"This is verified in the segment limited by the roots. \n");
		if (N0 < 0) {
			std::cout << "The polynomial in " << y << " is always positive,";
			NoSol();
			return;
		}
		std::cout << "The roots are: (-" << par(E1) << " - sqrt(" << E1 << sq <<
			" - 4*" << par(NegDisc) << "*" << par(F1) << ")) / (2*" << par(NegDisc) << ") = " <<
			R3 << "\n";
		std::cout << "and: (-" << par(E1) << " + sqrt(" << E1 << sq << " - 4*" << par(NegDisc)
			<< "*" << par(F1) << ")) / (2*" << par(NegDisc) + ") = " << R4 << " \n";
		if (R2<R1) {
			printf("There are no integers in this range,");
			NoSol();
			return;
		}

		std::cout << "All values of " << y << " from " << R1 << " to " << R2 << " should be replaced in \n";
		ShowEq(0, 0, NegDisc, 0, E1, F1, x1, y);
		printf(". The result should be the negative of a perfect square. \n");
		b = 0;  // flag used to control output formatting 
		for (u = R1; u <= R2; u++) {
			w1 = -NegDisc*u*u - E1*u - F1;
			w2 = llSqrt(w1);
			if (w2*w2 == w1) {   // w1 is a perfect square 
				if (b != 0) {
					/* not 1st solution */
					printf(", %lld", u);
				}
				else {
					/* first solution found*/
					std::cout << "The values of " << y << " are: " << u;
					b = 1;  /* set flag to show a solution was found */
				}
			}
		}
		if (b == 0) {
			std::cout << "This is not satisfied by any value of " << y;
			NoSol();
			return;
		}
		putchar('\n');
	}

	if (N0 < 0) {
		std::cout << "equation has no real solutions \n";
		NoSol();
		return;
	}

	/* try all integer values between R1 and R2*/
	bool found = false;
	for (u = R1; u <= R2; u++) {
		w1 = -NegDisc*u*u - E1*u - F1;
		w2 = llSqrt(w1);
		if (w2*w2 == w1) {
			/* we have found a solution */
			found = true;
			if (teach) {
				std::cout << y << " = " << u << "\n";
				std::cout << x1 << " = ";
				ShowLin(2 * g_A, g_B, g_D, x, y);
				std::cout << " = +/-Sqrt(" << (w2*w2) << ") = +/-" << w2 << "\n";
			}
			ShowElipSol(g_A, g_B, g_D, u, x, x1, y, w2);
			if (w2 != 0) {
				ShowElipSol(g_A, g_B, g_D, u, x, x1, y, -w2);
			}
		}
	}
	if (!found) {
		NoSol();
		return;
	}
}


// return value of bit in array corresponding to x. 0 indicates a prime number
static int getBit(const unsigned long long int x, uint32_t * array)
{
	unsigned long long int index;
	int bit;

	//index = x / 64;
	index = x >> 6;
	//bit = (x % 64) / 2;
	bit = (x & 0x3f) >> 1;   // get last 6 bits, then divide by 2
							 /* printf ("x = %d, index = %d, bit = %d, mask = %d\n", x, index, bit, mask) ; */
	return ((array[index] >> bit) & 1);

}

/* sets bit in 'array' corresponding to 'x' to 1 */
/* assume 'num' is odd. no bits are stored for even numbers */
/* assume an element of 'array' is 32 bits */
static int setBit(const unsigned long long int x, uint32_t * array) {
	unsigned long long int index;
	int bit, mask;

	//index = x / 64;
	index = x >> 6;
	//bit = (x % 64) / 2;
	bit = (x & 0x3f) >> 1;   // get last 6 bits, then divide by 2
	mask = (1 << bit);

	/* printf ("x = %d, index = %d, bit = %d, mask = %d\n", x, index, bit, mask) ; */
	array[index] |= mask; /* OR mask bit into array word */
	return (0);
}

void generatePrimes(unsigned long long int max_val) {
	unsigned long long int numsave, count = 1, num = 3;
	size_t plist_size = 0;
	unsigned long long int sqrt_max_val;
	int cutoff = 0;

	if (primeListMax >= max_val)
		return;  // nothing to do if prime list already generated

	max_val += 63 - (max_val + 63) % 64;  // round up to next multiple of 64
	sqrt_max_val = llSqrt(max_val);

	/* initialise flags */
	/* allocate 1 byte if max_val = 1-15, 2 bytes if 16-31, etc */
	if (primeFlags != NULL) free(primeFlags);	// if generatePrimes has been called before
												// clear out old prime list
	primeFlags = (uint32_t*)calloc((max_val / 16) + 1, 1);
	assert(primeFlags != NULL);

	// allocate storage for primeList if required
	{
		fprintf(stdout, "Expected no of primes is %.0f\n",
			(double)max_val / (log((double)max_val) - 1));
		if (primeList != NULL) free(primeList);
		plist_size = (size_t)((double)max_val / (log((double)max_val) - 1)) * 102 / 100;
		// add 2% for safety
		primeList = (unsigned long long *)malloc(plist_size * sizeof(long long));
		assert(primeList != NULL);
		prime_list_count = 1;
		primeList[0] = 2;  // put 1st prime into list
	}

	while (num <= max_val)
	{
		if (getBit(num, primeFlags) == 0)
		{   /* we have found a prime number */
			primeList[count] = num;
			count++;
			numsave = num;

			if (num <= sqrt_max_val)  /* check whether starting value for i is already
									  more than max_val, while avoiding an integer overflow*/
				for (unsigned long long int i = num*num; i <= max_val; i += (num * 2))
				{ /* mark all odd multiples of num in range as not prime */
				  /* starting value is num*num rather than num*3 because the flag bits for
				  multiples of num less than num*num are already set and using num*num
				  makes the program run significantly faster */
					setBit(i, primeFlags);
				}
		}
		num += 2;	// advance to next possible prime
	}

	// after completing the for loop we have found all the primes < max_val
	printf("  prime %9lld is %11lld\n", count, numsave);
	primeList[count] = ULLONG_MAX;		// set end marker
	prime_list_count = (unsigned int)count;
	primeListMax = primeList[count - 1];
	return;  // return value is used for problem 10
}
gmp_randstate_t state;
/* remove any double factors from G, add same factor to to K. The two versions 
are functionally the same but G has different types. This function, using a prime 
number list, was added to speed things up. If G is too large, we still have a 
problem. This method using a list of prime numbers is simple, but gets too slow 
for very large numbers. */
void adjustGandK(bigint &G, long long &K) {
	long long T;
	size_t i = 0;
	K = 1;
  	/* remove any double factors from G, add same factor to to K*/
	T = primeList[i];
	while (abs(G) >= T*T) {     // remove factor T
		while (G % (T*T) == 0) {
			G /= T*T; K *= T;
		}
		i++;
		if (i >= prime_list_count) {
			if (!mpz_probable_prime_p(ZT(G), state, 10, primeListMax)) {
				std::cerr << "Cannot factorise " << G << ". Try increasing max_prime \n";
				throw std::invalid_argument("Number too large to factorise");
			}
			else
				break;
		}
		T = primeList[i];     // advance to next prime
	}
}
void adjustGandK(long long &G, long long &K) {
	long long T;
	size_t i = 0;
	K = 1;
	/* remove any double factors from G, add same factor to to K*/
	T = primeList[i];
	while (abs(G) >= T*T) {     // remove factor T
		while (G % (T*T) == 0) {
			G /= T*T; K *= T;
		}
		i++;
		if (i >= prime_list_count) {
			std::cerr << "Cannot factorise " << G << ". Try increasing max_prime \n";
			throw std::invalid_argument("Number too large to factorise");
		}
		T = primeList[i];     // advance to next odd number
	}
}

/* determine type of equation. Also some basic checks for cases where there is no solution,
divide all coefficients by their gcd if gcd > 1, and calculate the discriminant */
equation_class	classify(const long long a, const long long b, const long long c,
	const long long d, const long long e, const long long f) {
	long long gcdA_E;
	bool DiscIsPerfSqare;
	bigint Bi_a, Bi_b, Bi_c;

	/* get gcd of A, B, C, D, E. gcd = zero only if they are all zero*/
	gcdA_E = gcd(a, gcd(b, gcd(c, gcd(d, e))));
	if (teach) {
		std::cout << "First of all we must determine the gcd of all coefficients but the constant term." << "\n";
		std::cout << "that is : gcd (" << a << ", " << b << ", " << c << ", " << d << ", " << e << ") = "
			<< gcdA_E << "\n";
	}
	if (gcdA_E != 0) {
		/* protect against divide-by-zero error */
		if (f%gcdA_E != 0) {
			NoGcd(f);    // output message - no solution 
			return no_soln;
		}
		else {
			/* divide all coefficients by gcd, set global values */
			g_A = a/gcdA_E;
			g_B = b/gcdA_E;
			g_C = c/gcdA_E;
			g_D = d/gcdA_E;
			g_E = e/gcdA_E;
			g_F = f/gcdA_E;
			if (teach && (gcdA_E != 1)) {
				std::cout << divgcd;    // show new values after division by gcd
				ShowEq(g_A, g_B, g_C, g_D, g_E, g_F, "x", "y");
				printf(" = 0\n");
			}
		}
	}
	else {  /* gcd = 0, presumably a-f are all 0 */
		g_A = a;           // x² coefficient
		g_B = b;           // xy coefficient 
		g_C = c;           // y² coefficient
		g_D = d;           // x coefficient
		g_E = e;           // y coefficient
		g_F = f;           // constant
	}

	if (g_D == 0 && g_A != 0 && g_C != 0)
		if (CheckMod(g_A, g_B, g_C, g_E, g_F))
			return no_soln;

	if (g_E == 0 && g_A != 0 && g_C != 0)
		if (CheckMod(g_B, g_C, g_A, g_D, g_F))
			return no_soln;

	Bi_a = g_A; Bi_b = g_B; Bi_c = g_C;   // copy A, B and C to extended-precision ints
	Bi_Disc = Bi_b*Bi_b - 4 * Bi_a*Bi_c; // get discriminant (now eliminated possibility of overflow)
	
	DiscIsPerfSqare = mpz_perfect_square_p(ZT(Bi_Disc));  // true if Disc is a perfect square

	if (Bi_Disc > 0 && !DiscIsPerfSqare &&
		g_D == 0 && g_E == 0 && g_F != 0) {
		if (teach)
			std::cout << "Equation is Hyperbolic Homogeneous\n";
		return hyperbolic_homog;  /* Bi_Disc is not a perfect square  and D, E are 0,
							and F is non-zero. This is a type of homogeneous equation*/
	}

	if (g_A == 0 && g_C == 0) {
		if (g_B == 0) {
			/* linear equation: A=B=C=0 */
			if (teach)
				std::cout << "Equation is Linear\n";
			return linear;
		}
		else {
			/* simple hyperbolic; A = C = 0; B ≠ 0. It follows that Disc > 0 */
			if (teach)
				std::cout << "Equation is Simple Hyperbolic\n";
			return simple_hyperbolic;
		}
	}

	/* not homogeneous, linear or simple hyperbolic equation */
	if (teach) {
		printf("We try now to solve this equation modulo 9, 16 and 25.\n");
	}
	if (Mod(9) || Mod(16) || Mod(25)) {
		return no_soln;
	}
	if (teach) {
		printf("There are solutions, so we must continue.\n");
	}

	if (Bi_Disc == 0) {
		if (teach)
			std::cout << "Equation is Parabolic\n";
		return parabolic;
	}

	if (Bi_Disc < 0) {
		if (teach)
			std::cout << "Equation is Elliptical\n";
		return elliptical;
	}

	if (!DiscIsPerfSqare) {
		if (teach)
			std::cout << "Equation is General Hyperbolic\n";
		return hyperbolic_gen;  // hyperbolic, not in other hyperbolic classes above
	}
	else {
		if (teach)
			std::cout << "Equation is Hyperbolic and the Discriminant is a perfect square\n";
		return hyperbolic_disc_ps;  // hyperbolic, disc is a perfect square
	}
}

/* Solve Diophantine equations of the form: Ax^2 + Bxy + Cy^2 + Dx + Ey + F = 0
This is the heart of the whole program.
there are six basic types which are solved using different methods:
•	Linear case:            A = B = C = 0.
•	Simple hyperbolic case: A = C = 0; B ≠ 0.
•	Elliptical case:        B^2 - 4AC < 0.
•	Parabolic case:			B^2 - 4AC = 0, A, B, C non-zero.
•	Hyperbolic case:        B^2 - 4AC > 0.
	 (1) Hyp. homogeneous equation Ax^2 + Bxy + Cy^2 + F = 0   i.e. D = 0, E = 0
	 (2) Hyp. general   equation:  Ax^2 + Bxy + Cy^2 + Dx + Ey + F = 0
Note: B^2 -4AC is known as the discriminant of the equation
*/
void solveEquation(const long long ax2, const long long bxy, const long long cy2,
	const long long dx, const long long ey, const long long f) {

	bool teachaux;
	long long G, H, K, T;
	std::string t1, x, y, x1, y1;

	/* code to check for heap corruption and memory leakage
	No code is generated unless compiled in debug mode */
	assert(_CrtCheckMemory());  //see https://msdn.microsoft.com/en-us/library/e73x0s4b.aspx
	_CrtMemCheckpoint(&memState1);  // see https://docs.microsoft.com/en-gb/cpp/c-runtime-library/reference/crtmemcheckpoint

									/* initialise some global variables */
	allSolsFound = false;
	NbrSols = 0;
	NbrCo = -1;
	also = false;
	//gmp_randinit_default(state);  // seed random number generator

	std::cout << "\nSolve Diophantine equation:   ";
	ShowEq(ax2, bxy, cy2, dx, ey, f, "x", "y");
	printf("  =  0\n by Dario Alejandro Alpern\n");

	equation_class eqnType = classify(ax2, bxy, cy2, dx, ey, f);   /* get equation type */

	if (g_A == 0 && g_C != 0) {
		T   = g_A;   /* swap A and C */
		g_A = g_C;
		g_C = T;
		T   = g_D;   /* swap D and E */
		g_D = g_E;
		g_E = T;
		ExchXY = true;
	}
	x = (ExchXY ? "y" : "x");
	y = (ExchXY ? "x" : "y");
	x1 = x + "'";
	y1 = y + "'";

	switch (eqnType) {
	case no_soln:
		std::cout << msg;   // no solutions
		break;

	case linear: { /* linear equation: A=B=C=0 */
		int rv = Linear(g_D, g_E, g_F);
		if (rv == 1) {
			std::cout << msg;   // no solutions
		}
		break;
	}

	case parabolic:   /* Parabolic case B^2 - 4AC = 0  A, B, C non-zero*/
		SolveParabolic(x, y, x1);
		break;

	case elliptical: {   /* B^2 - 4AC < 0 */
		SolveElliptical(x, x1, y);
		break;
	}

	case simple_hyperbolic:   // A = C = 0; B ≠ 0; therefore B^2 - 4AC > 0
		SolveSimpleHyperbolic();
		break;

		/* B^2 - 4AC > 0, D = 0, E = 0 */
	case hyperbolic_homog: {
		/* solve using continued fractions */
		bigint Disc2;
		teachaux = teach;  
		if (abs(g_F) != 1) {
			teach = false;   // save value of teach
		}
		GetRoot(g_A, g_B, g_C, &Disc2);       //  sneaky!! values returned in Disc, 
										 // SqrtDisc, Bi_NUM, Bi_DEN, used by SolContFrac
		teach = teachaux;    // restore saved value of teach

		G = H = g_F;
		adjustGandK(G, K);   /* remove any double factors from G, add same factor to to K*/
		/* K is the product of all the prime factors removed from G*/
		/* find all divisors of K */
		for (T = 1; T*T <= K; T++) {
			if (K%T == 0) {
				/* call SolContFrac for each factor of K  <= sqrt(K)*/
				SolContFrac(H, T, g_A, g_B, g_C, "");
			}
		}

		for (T = T - 1; T > 0; T--) {
			if (K%T == 0 && T*T < K) {
				/* call SolContFrac for each factor of K  > sqrt(K)*/
				SolContFrac(H, K / T, g_A, g_B, g_C, "");
			}
		}
		putchar('\n');

		if (also) {
			/* solutions found, so print them  */
			if (teach) {
				putchar('\n');
			}
			else {
				ShowAllLargeSolutions2();
				printf("If (x,y) is a solution, (-x,-y) is also a solution.\n");
			}
			ShowRecursion(eqnType);
		}
		else
			std::cout << msg;   // no solutions

		break;             // finished processing homogeneous equation
	}

	case hyperbolic_disc_ps:
	case hyperbolic_gen: {  /* equation fits none of the other categories */
		/* B^2 - 4AC > 0. There are a couple of special cases:
		(1)   B^2 - 4AC is a perfect square
		(2) N0 (calculated below) is zero
		otherwise the equation is solved using continued fractions */

		/* to avoid any risk of overflow the calculations below use Extended Precision! */
		bigint Disc2;
		bigint	 g, h, G, H;  
		const bigint NegDisc = -Bi_Disc;        // here, discriminant is +ve, so NegDisc is always -ve
		const bigint E1 = 4 * (bigint)g_A*g_E - (bigint)2 * g_B*g_D;
		const bigint F1 = 4 * (bigint)g_A*g_F - (bigint)g_D*g_D;
		gcd(NegDisc, E1/2, g);  // g = gcd of NegDisc and E1/2
		//const bigint Bi_CY1 = NegDisc/g;              // always a -ve number
		g_CY1 = NegDisc / g;              // always a -ve number
		//const bigint g_CY0 = E1/2/g;
		g_CY0 = E1 / 2 / g;
		const bigint BiN0 = g_CY0*g_CY0*g - g_CY1*F1;  

		gcd(g, BiN0, h);
		gcd(h, g_CY1, h);  // h = gcd of NegDisc, E1/2, BiN0, g_CY1
		//h = gcd(g_CY1, gcd(g, N0));       // h is gcd(Cy1, NegDisc, E1/2, N0)
		if (teach) {
			printf("We want to convert this equation to one of the form:\n");
			std::cout << x1 << sq << " + B " << y << sq << " + C " << y << " + D = 0 \n";
			std::cout << "Multiplying the equation by " << par(4 * g_A) << ":\n";
			ShowEq(4 * g_A*g_A, 4 * g_A*g_B, 4 * g_A*g_C, 4 * g_A*g_D, 4 * g_A*g_E, 4 * g_A*g_F, x, y);
			printf(" = 0\n");
			ShowLin(4 * g_A*g_A, 0, 0, x + sq, y);
			if (g_B != 0 || g_D != 0) {
				printf(" + (");
				ShowLin(0, 4 * g_A*g_B, 4 * g_A*g_D, x, y);
				std::cout << ")" << x;
			}
			if (g_C != 0 || g_E != 0 || g_F != 0) {
				printf(" + (");
				ShowEq(0, 0, 4 * g_A*g_C, 0, 4 * g_A*g_E, 4 * g_A*g_F, x, y);
				printf(") = 0\n");
			}
			if (g_B != 0 || g_D != 0) {
				printf("To complete the square we should add and subtract:\n(");
				ShowLin(0, g_B, g_D, x, y);
				std::cout << ")" << sq << "\nThen the equation converts to:\n(";
				ShowLin(2 * g_A, g_B, g_D, x, y);
				std::cout << ")" << sq << " + (";
				ShowEq(0, 0, 4 * g_A*g_C, 0, 4 * g_A*g_E, 4 * g_A*g_F, x, y);
				printf(") - (");
				ShowEq(0, 0, g_B*g_B, 0, 2 * g_B*g_D, g_D*g_D, x, y);
				printf(") = 0 \n");
			}

			printf("(");
			ShowLin(2 * g_A, g_B, g_D, x, y);
			std::cout << ")" << sq << " + (";
			ShowEq(0, 0, NegDisc, 0, E1, F1, x, y);
			printf(") = 0\nNow we perform the substitution:\n");
			std::cout << x1 << " = ";
			ShowLin(2 * g_A, g_B, g_D, x, y);
			printf("\nThis gives:\n");
			ShowEq(1, 0, NegDisc, 0, E1, F1, x1, y);
			printf(" = 0 \n");
		}
		if (teach) {
			if (NegDisc != g*h) {
				std::cout << "Multiplying the equation by " << (g_CY1 / h) << ":\n";
				ShowEq(g_CY1 / h, 0, NegDisc*g_CY1 / h, 0, NegDisc*E1 / g / h, NegDisc*F1 / g / h, x1, y);
			}
			printf(" = 0 \n");
			if (E1 != 0) {
				Show(g / h, "(", Show(g_CY1 / h, x1 + sq, 0));
				ShowEq(g_CY1*g_CY1, 0, 0, 2 * g_CY0*g_CY1, 0, 0, y, "");
				printf(")");
				Show1(NegDisc*F1 / g / h, 1);
				printf(" = 0 \n");
				Show(g / h, "(", Show(g_CY1 / h, x1 + sq, 0));
				if (g_CY1 != 1) {
					std::cout << par(g_CY1) << sq << " ";
				}
				std::cout << y << sq << " + 2*";
				if (g_CY1 != 1) {
					std::cout << par(g_CY1) << "*";
				}
				std::cout << par(g_CY0) << " " << y << ")";
				Show1(NegDisc*F1 / g / h, 1);
				printf(" = 0 \n");
				std::cout << "Adding and subtracting " << (g == h ? "" : numToStr(g / h) + " * ")
					<< par(E1 / 2 / g) + sq + ":\n";
			}
			Show(g / h, "(", Show(g_CY1 / h, x1 + sq, 0));
			if (g_CY1 != 1) {
				std::cout << par(g_CY1) << sq << " ";
			}
			std::cout << y << sq << " + 2*";
			if (g_CY1 != 1) {
				std::cout << par(g_CY1) << "*";
			}
			std::cout << par(g_CY0) << " " << y << " + " << par(g_CY0) << sq << ")";
			Show1(NegDisc*F1 / g / h, 1);
			std::cout << " - " << (g == h ? "" : numToStr(g / h) + " * ")
				<< par(E1 / 2 / g) << sq << " = 0 \n";
			Show(g / h, "(", Show(g_CY1 / h, x1 + sq, 0));
			ShowLin(0, g_CY1, g_CY0, x, y);
			std::cout << ")" << sq;
			Show1(-BiN0 / h, 1);
			std::cout << " = 0 \nMaking the substitution " << y1 << " = ";
			ShowLin(0, g_CY1, g_CY0, x, y);
			printf(":\n");
			ShowLin(g_CY1 / h, g / h, -BiN0 / h, x1 + sq, y1 + sq);
			printf(" = 0 \n");
		}

		const long long SqrtD = llSqrt(-NegDisc);
		const bigint N1 = abs(BiN0 / h);

		if (eqnType == hyperbolic_disc_ps) {    // is discriminant a perfect square?
			SolveDiscIsSq(BiN0, x, y, SqrtD, g, h, N1, y1, x1);
			break;
		}
		/* discriminant is not a perfect square */

		if (BiN0 == 0) {
			ShowX1Y1(0, 0, g_A, g_B, g_D, NegDisc, E1 / 2);   // display solution
			break;
		}

		/* Solve by continued fractions.
		Firstly, test if we need two cycles or four cycles */
		GetRoot(1, g_B, g_A * g_C, &Disc2);          //  sneaky!! values returned in Disc, 
								// SqrtDisc, Bi_NUM, Bi_DEN, used later by ContFrac
		ContFrac(g_A, 5, 1, 0, Bi_Disc, 1, g_A, Disc2, g_F); /* A2, B2 solutions */

		G = (2 * g_A2 + g_B*g_B2) % Disc2;    // note: g_A2, g_B2 set by ContFrac
		H = (g_B*g_A2 + 2 * g_A*g_C*g_B2) % Disc2;

		if (((g_C*g_D*(G - 2) + g_E*(g_B - H)) % Bi_Disc != 0 ||
			(g_D*(g_B - H) + g_A*g_E*(G - 2)) % Bi_Disc != 0) && ((g_C*g_D*(-G - 2) + g_E*(g_B + H)) % Bi_Disc != 0 ||
			(g_D*(g_B + H) + g_A*g_E*(-G - 2)) % Bi_Disc != 0)) {
			NbrCo *= 2;    
			//std::cout << "**temp SolveEquation  NbrCo=" << NbrCo << "\n";
		}

		//A = NegDisc/g/h;     
		//B = 0;                
		//C = g/h; 
		GetRoot(NegDisc/g/h, 0, g/h, &Disc2);  // values returned in Disc, SqrtDisc, Bi_NUM, Bi_DEN, used by SolContFrac

		G = H = -BiN0 / h;
		adjustGandK(G, K);  /* remove any double factors from G, add same factor to to K*/
		/* K is now the product of all the prime factors removed from G*/
		/* find all divisors of K */
		for (T = 1; T*T <= K; T++) {
			if (K%T == 0) {
				/* call SolContFrac for each divisor of K  <= sqrt(K) */
				SolContFrac(H, T, NegDisc/g/h, 0, (g/h), "'");
			}
		}
		for (T = T - 1; T > 0; T--) {
			if (K%T == 0 && T*T < K) {
				/* call SolContFrac for each divisor of K  > sqrt(K) */
				SolContFrac(H, K/T, NegDisc/g/h, 0, (g/h), "'");
			}
		}
		putchar('\n');
		if (also) {
			/* solutions found */
			if (!teach) {
				ShowAllLargeSolutions2();
			}
			ShowRecursion(eqnType);
		}
		else
			std::cout << msg;		// no solutions

		break;
	}

	default:
		assert(1 == 0);			// cause breakpoint if compiled in debug mode
		abort();				// must not happen
	}

	assert(_CrtCheckMemory());        /* check for heap corruption &  */
#ifdef _DEBUG
	mpz_realloc2(ZT(g_CY1), 32);   // avoid spurious memory leak reports, as far as possible
	mpz_realloc2(ZT(g_CY0), 32);
	mpz_realloc2(ZT(g_A2), 32);
	mpz_realloc2(ZT(g_B2), 32);
	mpz_realloc2(ZT(Bi_Disc), 32);
	mpz_realloc2(ZT(Bi_H1), 32);
	mpz_realloc2(ZT(Bi_H2), 32);
	mpz_realloc2(ZT(Bi_K2), 32);
	mpz_realloc2(ZT(Bi_NUM), 32);
	mpz_realloc2(ZT(Bi_DEN), 32);
#endif
	CheckForLeaks(); // check for memory leaks;
}


/* get a number from the input stream */
long long getnumber(const std::string &msg) {
	long long result;
	std::cout << msg;
	while (1) {
		std::cin >> result;     // as result is an integer only numeric input is allowed
		if (std::cin.good()) {
			break;			// valid number
		}
		else {
			// not a valid number
			std::cout << "Invalid Input! Please input a numerical value.  ";
			std::cin.clear();               // clear error flags
			std::cin.ignore(100000, '\n');  // discard invalid input up to newline/enter
		}
	}
	return result;
}

// generate a pseudo-random 64 bit integer.
static long long random64seed = 6497134987403476168;
void srandom64(unsigned long long seed) {
	random64seed = seed;
}

// generate a pseudo-random 64 bit integer.
/* The state must be seeded so that it is not everywhere zero. */
uint64_t random64(void) {
	static	uint64_t s[16] = { 464,463,822369, 78528, 6874, 16971, 27436, 13934645, 38572,
		2374348, 371941, 69429462649, 32716, 36586, 3258852, 2193 };
	static int p = 0;
	const uint64_t s0 = s[p];
	uint64_t s1 = s[p = (p + 1) & 15];

	s1 ^= s1 << 31; // a
	s[p] = s1 ^ s0 ^ (s1 >> 11) ^ (s0 >> 30); // b, c
	return s[p] * UINT64_C(1181783497276652981);
}


/* Entry point for program. Get parameters and call solveEquation */
int main(int argc, char* argv[]) {
	/* putting try here means that any exception thrown anywhere will be caught */
	try {
		long long int a, b, c, d, e, f;
		char yn = '\0';
		bool test = false;
		std::cout << "Compiled on " << __DATE__ " at " __TIME__ "\n";

		/* now generate the prime list used to factorise numbers. If max_prime  
		is too small the program will fail, but there has to be a limit. */
		//max_prime = 100000000;  // increase if necessary
		generatePrimes(max_prime);
		while (toupper(yn) != 'Y' && toupper(yn) != 'N') {
			std::cout << "Enter own data (Y/N)? (N = run standard tests) : ";
			std::cin >> yn;
			if (std::cin.good()) {
				test = (toupper(yn) == 'N');
			}
			else 
				return EXIT_FAILURE;   //  error reading input
		}

		if (!test) {
			std::cout << "Solve Diophantine equations of the form: Ax^2 + Bxy + Cy^2 + Dx + Ey + F = 0" << "\n";
			a = getnumber("Enter value for A ");
			b = getnumber("Enter value for B ");
			c = getnumber("Enter value for C ");
			d = getnumber("Enter value for D ");
			e = getnumber("Enter value for E ");
			f = getnumber("Enter value for F ");
			//std::cout << "solve " << a << "x² + " << t << "xy + " << cy2 << "y² + " << dx << "x + " << ey << "y + " << f << " = 0" << endl;
			wprintf(L"solve %lldx^2 + %lldxy + %lldy^2 + %lldx + %lldy + %lld = 0\n", a, b, c, d, e, f);

			yn = '\0';
			while (toupper(yn) != 'Y' && toupper(yn) != 'N') {
				std::cout << "Detailed explanation required? (Y/N): ";
				std::cin >> yn;
				if (std::cin.good()) {
					teach = (toupper(yn) == 'Y');
				}
				else
					return EXIT_FAILURE;   //  error reading input
			}

			solveEquation(a, b, c, d, e, f);
			system("PAUSE");   // press any key to continue
			return EXIT_SUCCESS;
		}

		else {
			yn = '\0';
			while (toupper(yn) != 'Y' && toupper(yn) != 'N') {
				std::cout << "Detailed explanation required? (Y/N): ";
				std::cin >> yn;
				if (std::cin.good()) {
					teach = (toupper(yn) == 'Y');
				}
				else
					return EXIT_FAILURE;   //  error reading input
			}

			printf("Run standard tests\n");

			printf("\n1st test. No solution exists \n");
			a = 0; b = 0; c = 0; d = 10; e = 84; f = 15;
			solveEquation(a, b, c, d, e, f);

			/* example numbers below refer to Dario Alpert's document. */
			printf("\nExample 1 (linear)\n");
			/* x = -136 + 42t,  y = 16 - 5t,  where t is any integer number.
			equivalent to: x = -10 +42t     y = 1-5t      (substitute t+3 for t) */
			a = 0; b = 0; c = 0; d = 10; e = 84; f = 16;
			solveEquation(a, b, c, d, e, f);

			printf("\nExample 2 (simple hyperbolic)\n");
			a = 0; b = 2; c = 0; d = 5; e = 56; f = 7;
			/*
			x = (266-56)/2  =  105, y = (266/266-5)/2    =  -2
			x = (-266-56)/2 = -161, y = [266/(-266)-5]/2 =  -3
			x = (2-56)/2    =  -27, y = (266/2-5)/2      =  64
			x = (-2-56)/2   =  -29, y = [266/(-2)-5]/2   = -69
			x = (38-56)/2   =   -9, y = (266/38-5)/2     =   1
			x = (-38-56)/2  =  -47, y = [266/(-38)-5]/2  =  -6
			x = (14-56)/2   =  -21, y = (266/14-5)/2     =   7
			x = (-14-56)/2  =  -35, y = [266/(-14)-5]/2  = -12
			The only 8 solutions to the equation are above
			*/
			solveEquation(a, b, c, d, e, f);

			printf("\nExample 3 (elliptical)\n");
			a = 42; b = 8; c = 15; d = 23; e = 17; f = -4915;
			/* x = -11, y=-1 is the only solution */
			solveEquation(a, b, c, d, e, f);

			printf("\nExample 4 (parabolic)\n");
			a = 8; b = -24; c = 18; d = 5; e = 7; f = 16;
			/*
			set 1: x = -174t^2 - 17t - 2, y = -116t^2 - 21t - 2
			set 2: x = -174t^2 - 41t - 4, y = -116t^2 - 37t - 4
			*/
			solveEquation(a, b, c, d, e, f);

			printf("\nProject Euler problem 140 (general hyperbolic)\n");
			a = -1; b = 0; c = 5; d = 0; e = 14; f = 1;
			solveEquation(a, b, c, d, e, f);
			/*
			basic solutions are:
			x=±1, 	y=0
			x=±2, 	y=-3
			x=±5, 	y=-4
			x=±7, 	y=2

			also, for each of the solutions above we can generate an infinite number of other
			solutions using the formulas:
			Xn+1 = PXn + QYn + K
			Yn+1 = RXn + SYn + L
			where:
			P = -9
			Q = 20		or -20
			K = 28		or -28
			R = 4		or -4
			S = -9
			L = -14
			*/

			printf("\nExample 5 (hyperbolic homogeneous)\n");
			//teach = true;
			a = 18; b = 41; c = 19; d = 0; e = 0; f = -24;
			/*
			basic solutions are:
			x = -10,							y =  6
			X = -7								Y =  11
			X = -202							Y =  312
			X = -14267							Y =  8751
			x = -10130							Y =  15646
			X = -284123							Y =  438834
			X =  -4 680127 (7 digits)			Y =  2 870666 (7 digits
			X =  -14 247838 (8 digits)			Y =  22 006088 (8 digits)
			If (x,y) is a solution, (-x,-y) is also a solution.
			*/
			solveEquation(a, b, c, d, e, f);

			printf("\nExample 6\n");
			a = 3; b = 13; c = 5; d = -11; e = -7; f = -92;
			/* x=-4, y=0, or  x=2, y=3 etc
			P=8351, Q=32625, R=-19575, S=-76474, K=-28775, L=67450*/
			solveEquation(a, b, c, d, e, f);

			printf("\nExample 6 - modified to produce 2 sets of values for P, Q, R, ...\n");
			/* x=-4, y=0, x=12, y=-2 etc
			P= 8351, Q=32625, R=-19575, S=-76474, K=-775, L= 1825
			or P=-8351, Q=32625, R= 19575, S= 76474, K= 783, L=-1827*/
			a = 3; b = 13; c = 5; d = -11; e = -42; f = -92;
			solveEquation(a, b, c, d, e, f);

			printf("\nExample 7\n");
			/* x=4, y=7 etc
			P = -1188641 Q = -4979520 K = 5146869 R = 2489760 S = 10430239 L = -10780770 */
			a = 3; b = 14; c = 6; d = -17; e = -23; f = -505;
			solveEquation(a, b, c, d, e, f);

			a = 2; b = 5; c = 2; d = 6; e = 6; f = 4;
			std::cout << "\nproduct of two linear expressions\n";
			std::cout << "(general hyperbolic with discriminant a perfect square)";
			solveEquation(a, b, c, d, e, f);

			/* D == 0 && A != 0 && C != 0*/
			std::cout << "Hyperbolic homogeneous\n";
			a = 18; b = 41; c = 19; d = 0; e = 0; f = 13;
			solveEquation(a, b, c, d, e, f);

			std::cout << "\nElliptical, no solutions\n";
			a = -1; b = -1; c = -2; d = -3; e = -4; f = -5;
			solveEquation(a, b, c, d, e, f);

			std::cout << "\nElliptical, no solutions\n";
			a = 1; b = 0; c = 1; d = -0; e = 0; f = -6;
			solveEquation(a, b, c, d, e, f);

			std::cout << "\nHyperbolic homogeneous - no solutions\n";
			a = 1; b = 0; c = -34; d = -0; e = 0; f = 1;
			solveEquation(a, b, c, d, e, f);

			std::cout << "\nHyperbolic, not solved using continued fractions\n";
			a = -1000; b = -1000; c = -249; d = -1000; e = -570; f = 975;
			solveEquation(a, b, c, d, e, f);

			printf("\nExample 4a (parabolic)\n");
			a = 8; b = -24; c = 18; d = -2; e = 3; f = 0;
			/*
			set 1: x = -174t^2 - 17t - 2, y = -116t^2 - 21t - 2
			set 2: x = -174t^2 - 41t - 4, y = -116t^2 - 37t - 4
			*/
			solveEquation(a, b, c, d, e, f);

			teach = false;
			printf("\n test using large pseudo-random numbers as coefficients\n");
			a = random64();    a >>= 48;    // note that in Visual Studio an arithmetic right shift is used
			b = random64();    b >>= 48; 
			c = random64();    c >>= 48;
			d = random64();    d >>= 48;
			e = random64();    e >>= 48;
			f = 0;              // this  ensures that there is at least 1 solution; x=y=0
			solveEquation(a, b, c, d, e, f);

			printf("\n 2nd test using large pseudo-random numbers as coefficients\n");

			a = random64();    a >>= 48;
			b = random64();    b >>= 48; 
			c = random64();    c >>= 48;  c = -c;
			d = random64();    d >>= 48;
			e = random64();    e >>= 48;
			f = 0;               // this  ensures that there is at least 1 solution; x=y=0
			solveEquation(a, b, c, d, e, f);


			printf("\n 3rd test using large pseudo-random numbers as coefficients\n");
			teach = true;
			a = random64();    a >>= 32; a = -a;
			b = random64();    b >>= 48;
			c = random64();    c >>= 32;  
			d = random64();    d >>= 48;
			e = random64();    e >>= 16;
			f = 0;               // this  ensures that there is at least 1 solution; x=y=0
			solveEquation(a, b, c, d, e, f);


			printf("\n 4th test using large pseudo-random numbers as coefficients\n");
			teach = false;
			a = random64();    a >>= 47; 
			b = random64();    b >>= 47;
			c = random64();    c >>= 47;
			d = random64();    d >>= 47;
			e = random64();    e >>= 47;
			f = 0;               // this  ensures that there is at least 1 solution; x=y=0

			solveEquation(a, b, c, d, e, f);


			// code to find values for hyperbolic eqn where N0=0
	/*		for (a = -1000; a <=1000; a++)
				for (b = -1000; b <= 1000; b++)
					for (c = -1000; c <= 1000; c++) {
						g_Disc = b*b - 4 * a*c;
						if (g_Disc <= 0)
							continue;
						for (d = -1000; d <= 1000; d++)
							for (e = -1000; e <= 1000; e++) {
								long long e1 = 4 * a*e - 2 * b*d;
								long long g = gcd(g_Disc, e1 / 2);
								g_CY0 = e1 / 2 / g;
								for (f = -1000; f <= 1000; f++) {
									long long f1 = 4 * a*f - d*d;
									g_CY1 = -g_Disc / g;
									long long N0 = g_CY0*g_CY0*g - g_CY1*f1;
									if (N0 == 0)
										std::cout << "a=" << a << " b=" << b << " c=" << c
										<< "  d=" << d << " e=" << e << " f=" << f << "\n";
								}
							}
					}*/
		}

		system("PAUSE");   // press any key to continue
		return EXIT_SUCCESS;
	}

	/* code below is executed only if an exception occurs*/
	catch (const std::exception &exc)
	{
		// catch anything thrown within try block that derives from std::exception
		std::cerr << "\n** Exception:" << exc.what() << std::endl;
		Beep(1200, 1000);              // sound at 1200 Hz for 1 second
		system("PAUSE");               // press any key to continue
		exit(EXIT_FAILURE);
	}

	catch (const char *str)
	{
		std::cerr << "Caught exception: " << str << std::endl;
		Beep(1200, 1000);              // sound at 1200 Hz for 1 second
		system("PAUSE");               // press any key to continue
		exit(EXIT_FAILURE);
	}

	catch (...) {   // catch block probably only would be executed under /EHa 
					/* most likely to be a SEH-type exception */
					//Eval_Exception(GetExceptionCode());
		std::cerr << "Caught unknown exception in catch(...)." << std::endl;
		Beep(1200, 1000);              // sound at 1200 Hz for 1 second
		system("PAUSE");   // press any key to continue
		exit(EXIT_FAILURE);
	}
}

/* SEH exception handler */
//int Eval_Exception(int n_except) {
//	if (n_except != STATUS_INTEGER_OVERFLOW &&
//		n_except != STATUS_INTEGER_DIVIDE_BY_ZERO )   // Pass on most exceptions  
//		return EXCEPTION_CONTINUE_SEARCH;
//
//	// Execute some code to handle problems  
//	if (n_except == STATUS_INTEGER_OVERFLOW) {
//		std::cout << "\n** INTEGER_OVERFLOW - program continues ** \n";
//		Beep(1000, 1000);
//		return EXCEPTION_CONTINUE_EXECUTION;
//	}
//
//	if (n_except == STATUS_INTEGER_DIVIDE_BY_ZERO) {
//		std::cout << "\n** INTEGER DIVIDE BY ZERO - program terminates ** \n";
//		Beep(1200, 1000);              // sound at 1200 Hz for 1 second
//		exit(EXIT_FAILURE);
//	}
//
//
//	return EXCEPTION_CONTINUE_EXECUTION;
//}


