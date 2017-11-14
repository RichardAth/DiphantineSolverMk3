#define _CRTDBG_MAP_ALLOC  
#include <stdlib.h>  
#include <crtdbg.h> 
#include "stdafx.h"

// get access to the mpz_t inside  mpz_int a
#define ZT(a) a.backend().data()

const std::string sq = "^2";
extern mpz_int Bi_L1, Bi_L2, Bi_H1, Bi_H2, Bi_K1, Bi_K2;
extern long long g_Disc;
extern long long g_A, g_B, g_D, g_F;
extern mpz_int g_NUM, g_DEN;
extern unsigned long long SqrtDisc;
extern long long g_CY1, g_CY0;
extern long long g_A2, g_B2;
extern int NbrSols, NbrCo;
extern bool teach;

long long g_A1, g_B1;
int NbrEqs, EqNbr;
std::string UU = "";
std::string VU = "";
std::string UL = "";
std::string VL = "";
std::string UL1 = "";
std::string VL1 = "";
std::string FP = "";

/*******************************************/
/* NextConv:                               */
/*  BigInteger tmp = Prev * A1 + Act * B1; */
/*  Act = Prev * A2 + Act * B2;            */
/*  Prev = Tmp;                            */
/*******************************************/
void NextConv(mpz_int *Bi_Prev, mpz_int *Bi_Act, const long long A1,
	const long long A2, const long long B1, const long long B2) {
	mpz_int t1, t2, tmp;

	/*std::cout << "**temp NextConv: Prev = " << *Bi_Prev;
	std::cout << " Act ="<< *Bi_Act;
	std::cout << " A1 =" << A1 << " A2=" << A2 << " B1=" << B1 << " B2=" << B2 << "\n";	*/

	tmp = *Bi_Prev * A1 + *Bi_Act * B1;
	*Bi_Act = *Bi_Prev * A2 + *Bi_Act * B2;
	*Bi_Prev = tmp;

	/* temporary */
	/*std::cout << "**exit NextConv: Prev = " << *Bi_Prev;
	std::cout << " Act ="<< *Bi_Act << "\n\n";*/
	/* end temporary */
}

/* uses global variables Bi_L1, Bi_L2, g_A, g_B, g_D, g_F, g_CY0, g_CY1*/
bool ShowHomoSols(int type, mpz_int Bi_SHH, mpz_int Bi_SHK, long long s, long long T,
	const long long MagnifY, const std::string eqX, const std::string eqY) {

	assert(_CrtCheckMemory());

	/*std::cout << "**temp ShowHomoSols: s=" << s << "  T=" << T << "  MagnifY=" << MagnifY << "\n";
	std::cout << "Bi_SHH=" << Bi_SHH;
	std::cout << "  Bi_SHK=" << Bi_SHK << "\n";
*/
	int i;
	std::string U = (type == 4 ? "'" : "");
	std::string X1 = "X";
	std::string Y1 = (MagnifY == 1 ? "Y" : "Y'");

	if (teach) {
		ShowLargeXY(Y1, "Z", Bi_SHH, Bi_SHK, true, eqX, eqY);
		std::cout << "Since " << X1 << U << "0 = ";
		if (T>1) {
			std::cout << T << " " << UU << U << "0 = " << T << " (";
		}
		ShowLin(s, -g_F, 0, VU + "0", "Z0");
		if (T>1) {
			std::cout << ")\nand " << Y1 << "0 = " << T << " V0";
		}
		printf(":\n");
	}
	//MultAddLargeNumbers(s*T, Bi_SHH, -g_F*T, Bi_SHK, &Bi_L1);  // result stored in Bi_L1
	Bi_L1 = s*T* Bi_SHH - g_F*T* Bi_SHK;
	Bi_L2 = T* Bi_SHH; 

	/*std::cout << "**temp ShowHomoSols(1) Bi_L1=" << Bi_L1;
	std::cout << "  Bi_L2=" << Bi_L2 << "\n";*/

	if (type == 4) {
		if (teach) {
			ShowLargeXY(X1 + U, Y1 + U, Bi_L1, Bi_L2, false, "", "");
		}
		for (i = ((Bi_L1 == 0) && (Bi_L2 == 0) ? 1 : 0); i<2; i++) {
			Bi_L2 -= g_CY0; //AddLarge(Bi_L2, -g_CY0, &Bi_L2);
			//std::cout << "**temp ShowHomoSols(2) Bi_L2=" << Bi_L2 << "\n";

			if (teach) {
				std::cout << Y1 << "0 = (";
				ShowLin(0, i == 0 ? 1 : -1, -g_CY0, "", Y1 + "0");
				std::cout << ")/" << par(g_CY1) << "\n";
			}
			if (tDivLargeNumber(Bi_L2, g_CY1, &Bi_L2) != 0) {
				//std::cout << "**temp ShowHomoSols(3) Bi_L2=" << Bi_L2 << "\n";
				if (teach) {
					printf("It is not an integer number. \n");
				}
			}
			else {
				//std::cout << "**temp ShowHomoSols(4) Bi_L2=" << Bi_L2 << "\n";
				if (teach) {
					std::cout << X1 << "0 = (";
					ShowLin(i == 0 ? 1 : -1, -g_B, -g_D, X1 + "0", Y1 + "0");
					std::cout << ") / " << par(2 * g_A) << "\n";
				}
				Bi_L1 -= g_D; 
				//MultAddLargeNumbers(1, Bi_L1, -g_B, Bi_L2, &Bi_L1);  // store result in L1
				Bi_L1 -=  g_B* Bi_L2;
				//std::cout << "**temp ShowHomoSols(5) Bi_L1=" << Bi_L1 << "\n";
				if (tDivLargeNumber(Bi_L1, 2 * g_A, &Bi_L1) != 0) {
					//std::cout << "**temp ShowHomoSols(6) Bi_L1=" << Bi_L1 << "\n";
					if (teach) {
						printf("It is not an integer number. \n");
					}
				}
				else {
					//std::cout << "**temp ShowHomoSols(7) Bi_L1=" << Bi_L1 << "\n";
					if (teach) {
						if (MagnifY != 1) {
							ShowLargeXY(X1, Y1, Bi_L1, Bi_L2, false, "", "");
							std::cout << "Since Y = " << MagnifY << Y1;
						}
						std::cout << "---------------------------------------";
					}
					Bi_L2 *= MagnifY; 
					ShowLargeXY("X", "Y", Bi_L1, Bi_L2, false, "", "");
					if (teach)
						std::cout << "---------------------------------------\n";
					/*std::cout << "**temp ShowHomoSols(8) returns true: Bi_L1=" << Bi_L1;
					std::cout << "  Bi_L2=" << Bi_L2 << "\n";*/

					return true;
				}
			}
			/*MultAddLargeNumbers*/(-s*T, Bi_SHH, g_F*T, Bi_SHK, &Bi_L1);
			Bi_L1 = -s*T* Bi_SHH + g_F*T* Bi_SHK;
			/*gmp_printf("**temp ShowHomoSols(9) Bi_L1= %lld*%Zd + %lld*%Zd =", -s*T, ZT(Bi_SHH), g_F*T, ZT(Bi_SHK));
			std::cout << Bi_L1;*/

			Bi_L2 = -T*Bi_SHH; 
			//std::cout << "  Bi_L2=" << Bi_L2 << "\n";
		}
	}
	else {
		if (teach) {
			if (MagnifY != 1) {
				ShowLargeXY(X1, Y1, Bi_L1, Bi_L2, false, "", "");
				std::cout << "Since " << Y1 << " = " << MagnifY << "Y";
			}
			putchar('\n');
		}
		Bi_L2 *= MagnifY; // MultLargeNumber(MagnifY, Bi_L2, &Bi_L2);

		if (teach) {
			ShowLargeXY("X", "Y", Bi_L1, Bi_L2, false, "", "");
			ShowLargeXY("X", "Y", -Bi_L1, -Bi_L2, false, "", "");
		}

	/*	std::cout << "**temp ShowHomoSols(10) returns true: Bi_L1=" << Bi_L1;
		std::cout << "  Bi_L2=" << Bi_L2 << "\n";*/

		return true;
	}
	//std::cout << "**temp ShowHomoSols(11) returns false \n";
	return false;   // solution not found
}

/***************************************************************************
* type = 1: Find convergents                                               *
* type = 2: Find convergents for x^2 + Bxy + ACy^2 = 1 (recursion)         *
* type = 3: Find convergents for modified equation in homogeneous equation *
* type = 4: Find convergents for modified equation in complete solution    *
* type = 5: Find convergents for x^2 + Bxy + ACy^2 = 1 (mod B^2-4AC)       *
* returns true if there are solutions, otherwise false                     *
* uses global variables DP_NUM, DP_DEN, Bi_H1, Bi_H2, Bi_K1, Bi_K2, NbrCo  *
*     NbrEqs, Eqnbr, NbrSols, g_F, g_A1, g_A2, g_B1, g_B2                  *
****************************************************************************/
bool ContFrac(const mpz_int Dp_A, int type, const int SqrtSign, long long s, long long T,
	long long MagnifY, long long A) {

	/*std::cout << "**temp ContFrac: type=" << type << "  SqrtSign=" << SqrtSign << "  s=" << s;
	std::cout << "  T=" << T << "  MagnifY=" << MagnifY << "  A=" << A << "\n";
	std::cout << "g_NUM=" << g_NUM;
	std::cout << "  g_DEN=" << g_DEN << "\n";*/

	long long P, Z, M, P1, M1, Tmp, K, L, Mu;
	mpz_int Dp_P, Dp_M, Dp_Z, Dp_G, Dp_Mu, Dp_K, Dp_L, Dp_M1, Dp_P1, Dp_zz;
	long long H1ModCY1 = 1, H2ModCY1 = 0, K1ModCY1 = 0, K2ModCY1 = 1;
	bool Sols = true, secondDo = true;
	int Conv;
	int Co = -1;
	std::string U = (type == 4 ? "'" : "");
	std::string X1 = "X";
	std::string Y1 = (MagnifY == 1 ? "Y" : "Y'");

	assert(_CrtCheckMemory());

	if (Dp_A == 1) {
		/* Dp_A = 1 */
		Bi_H1 = SqrtSign; 
		Bi_K1 = 0;       
		//std::cout << "**temp ContFrac Bi_H1=" << Bi_H1 << "\n";

		if (type == 1) {
			ShowLargeXY(X1, Y1, Bi_H1, Bi_K1, true, "", "");
			//std::cout << "** temp ContFrac returns true (1) - solution(s) found\n";
			assert(_CrtCheckMemory());
			return true;            /* Indicate there are solutions */
		}

		if ((type == 3 || type == 4) && (g_Disc != 5 || A*g_F<0)) {
			if (ShowHomoSols(type, Bi_H1, Bi_K1, s, T, MagnifY, "", "")) {
				//std::cout << "**temp ContFrac(2) - solution found\n";
				assert(_CrtCheckMemory());
				return true;          /* Indicate there are solutions */
			}
		}
	}

	/* Paso = 1: Quick scan for solutions */
	/* Paso = 2: Show actual solutions */
	for (int Paso = (type == 2 || g_Disc == 5 && A*g_F>0 && (type == 3 || type == 4) ? 2 : 1);
		Sols && Paso <= 2; Paso++) {
		Conv = 0;
		Sols = false;

		Dp_P = g_DEN; 
		if (SqrtSign < 0) {
			Dp_P = -Dp_P;
		}
		Dp_K = SqrtDisc + (Dp_P < 0 ? 1 : 0); 
		//std::cout << "  Dp_K=" << Dp_K << " (1)\n";
		if (SqrtSign < 0) {
			Dp_K -= g_NUM; 
			//std::cout << "  Dp_K=" << Dp_K << " (2)\n";
		}
		else {
			Dp_K += g_NUM;       
				//std::cout << "  Dp_K=" << Dp_K << " (3)\n";
		}

		//std::cout << "  g_DEN=" << g_DEN;
		Z = DivDoublePrec(Dp_K, Dp_P);           // Z = K/P
		//std::cout << "  Dp_K=" << Dp_K << "  Dp_P=" << Dp_P << "\n";
		Dp_M = Z;       // M = Z (=K/P)
		Dp_K = Dp_M * g_DEN;   // K = M*DEN
		//std::cout << "  Dp_K=" << Dp_K << " Z=" << Z;

		Dp_M = Dp_K - g_NUM;       
		//std::cout << "  Dp_M=" << Dp_M << " (4)\n";

		if (SqrtSign < 0) {
			Dp_M = -Dp_M; 
		}

		/* type = 4: Find convergents for modified equation in complete solution */
		if (type == 4) {
			H2ModCY1 = Z%g_CY1;
		}

		/* type = 5: Find convergents for x^2 + Bxy + ACy^2 = 1 (mod B^2-4AC) */
		if (type == 5) {
			g_A1 = g_B2 = 1;
			g_A2 = Z%T;
			g_B1 = 0;
		}
		else {
			Bi_H1 = SqrtSign;		
			Bi_H2 = Z*SqrtSign;		
			Bi_K1 = 0;				
			Bi_K2 = Bi_H1;			
			/*std::cout << "**temp ContFrac Bi_H1=" << Bi_H1;
			std::cout << "  Bi_H2=" << Bi_H2;
			std::cout << "  Bi_K1=" << Bi_K1;
			std::cout << "  Bi_K2=" << Bi_K2 << "\n";*/

			g_A1 = g_B2 = 1;
			g_A2 = g_B1 = 0;
		}

		Co = -1;

		Dp_K = -1;   
		Dp_L = -1;  
		//std::cout << "  Dp_K=" << Dp_K << " (6)\n";
		/* set Mu */
		switch (type) {
		case 1:   // find convergents 
			Dp_Mu = -2 * g_F*SqrtSign; 
			break;

		case 3:     // find convergents for modified equation in homogeneous equation
		case 4:     // find convergents for modified equation in complete solution
			Dp_Mu = -2 * SqrtSign; 
			break;

		default:  // type = 2 or 5
			Dp_Mu = g_DEN;    
			if (SqrtSign > 0) {
				Dp_Mu = -Dp_Mu; 
				//std::cout << "**temp ContFrac(3A) Dp_Mu= " << Dp_Mu << "  (revsign)\n";
			}
		}

		/*std::cout << "**temp ContFrac(3B) Dp_Mu= " << Dp_Mu;
		std::cout << "  Dp_K=" << Dp_K;
		std::cout << "  Dp_P1=" << Dp_P1 << "\n";*/

		do {
			Dp_Z = g_Disc;  
			//std::cout << "**temp ContFrac  Dp_Z=" << Dp_Z << " (3D)\n";

			Dp_G = Dp_M * Dp_M; 
			//std::cout << "**temp ContFrac  Dp_G=" << Dp_G << "(3D)\n";

			Dp_G = Dp_Z - Dp_G;     // G = Z-G  = Disc -M*M
			//std::cout << "**temp ContFrac  Dp_G=" << Dp_G << " (3E)\n";

			DivideDoublePrecLong(Dp_G, Dp_P, &Dp_P1);  // P1 = (Disc-M*M)/P
			 //std::cout << "**temp ContFrac Dp_P1=" << Dp_P1 << "\n";
			/* Z = SqrtDisc +(1 or 0, depending on sign of P1) */
			Dp_Z = SqrtDisc + ((Dp_P1 < 0) ? 1 : 0); 
			Dp_K = Dp_M + Dp_Z;        
			//std::cout << "**temp ContFrac  Dp_K=" << Dp_K << "\n";

			/* round Z to a multiple of P1 */
			Z = DivDoublePrec(Dp_K, Dp_P1);        // Z = K/P1
			Dp_G = Z;         // G = Z = K/P1
			Dp_Z = Dp_G * Dp_P1;   
			//std::cout << "**temp ContFrac  Dp_Z=" << Dp_Z << "\n";

			Dp_M1 = Dp_Z - Dp_M; 
			//std::cout << " **temp ContFrac Dp_M1=" << Dp_M1 << "\n";

			Dp_zz = SqrtDisc + Dp_M;    
			if (Co<0 && (Dp_P > 0) && (Dp_P < Dp_zz) &&  (Dp_M > 0) && (Dp_M <= SqrtDisc)) {  
				Co = 0;
				Dp_K = Dp_P;    
				Dp_L = Dp_M;  
				/*std::cout << "**temp ContFrac(4) Dp_K=" << Dp_K;
				std::cout << "  Dp_L=" << numToStr(Dp_L) << "\n";*/
			}

			/*std::cout << "**temp ContFrac(5)  type=" << type << "  Co=" << Co;
			std::cout << "  Dp_P="  << Dp_P;
			std::cout << "  Dp_Mu=" << Dp_Mu;
			std::cout << "  Dp_K="  << Dp_K;
			std::cout << "  Dp_P1=" << Dp_P1 << "\n";*/

			if (type == 1 && Dp_P == Dp_Mu) {
				// Solution found
				if (Co % 2 == 0 || (Dp_K != Dp_P1)) {
					if (Paso == 2) {
						if (g_A2 != 0) {
							NextConv(&Bi_H1, &Bi_H2, g_A1, g_A2, g_B1, g_B2);
							NextConv(&Bi_K1, &Bi_K2, g_A1, g_A2, g_B1, g_B2);
							g_A1 = g_B2 = 1;
							g_A2 = g_B1 = 0;
						}
						ShowLargeXY(X1, Y1, Bi_H2, Bi_K2, true, "NUM(" + numToStr(Conv) + ") = ",
							"DEN(" + numToStr(Conv) + ") = ");
					}
					Sols = true;
				}
				secondDo = false;
				//std::cout << "**temp ContFrac(5A)\n";
				break;
			}

			if (type == 3 || type == 4) {
				if (Co == 0 && A*g_F>0 && g_Disc == 5) {  /* Solution found */
					if (Paso == 1) {
						secondDo = false;
						Sols = true;
						//std::cout << "**temp ContFrac(6)\n";
						break;
					}
					else {
						NextConv(&Bi_H1, &Bi_H2, g_A1, g_A2, g_B1, g_B2);
						NextConv(&Bi_K1, &Bi_K2, g_A1, g_A2, g_B1, g_B2);
						g_A1 = g_B2 = 1;
						g_A2 = g_B1 = 0;
						Bi_H2 -= Bi_H1; 
						Bi_K2 -= Bi_K1; 
						if (ShowHomoSols(type, Bi_H2, Bi_K2, s, T, MagnifY,
							"NUM(" + numToStr(Conv) + ") - NUM(" + numToStr(Conv - 1) + ") = ",
							"DEN(" + numToStr(Conv) + ") - DEN(" + numToStr(Conv - 1) + ") = ")) {
							secondDo = false;
							Sols = true;
							//std::cout << "**temp ContFrac(7) - solution found\n";
							break;
						}
						Bi_H2 += Bi_H1;   
						Bi_K2 += Bi_K1;   
					}
				}
				if (Dp_P1 == Dp_Mu) {
					// Solution found
					if (Co % 2 == 0 || (Dp_K != Dp_P1) || (Dp_L != Dp_M1)) {
						if (Paso == 2) {
							if (g_A2 != 0) {
								NextConv(&Bi_H1, &Bi_H2, g_A1, g_A2, g_B1, g_B2);
								NextConv(&Bi_K1, &Bi_K2, g_A1, g_A2, g_B1, g_B2);
								g_A1 = g_B2 = 1;
								g_A2 = g_B1 = 0;
							}
							if (ShowHomoSols(type, Bi_H2, Bi_K2, s, T, MagnifY, "NUM(" + numToStr(Conv) + ") = ",
								"DEN(" + numToStr(Conv) + ") = ")) {
								secondDo = false;
								Sols = true;
								//std::cout << "**temp ContFrac(8) - solution found\n";
								break;
							}
						}
						else {
							if (type == 4) {
								Tmp = H2ModCY1*T;
								if ((Tmp - g_CY0) % g_CY1 == 0 || (Tmp + g_CY0) % g_CY1 == 0) {
									secondDo = false;
									Sols = true;
									//std::cout << "**temp ContFrac(9)\n";
									break;
								}
							}
							else {
								secondDo = false;
								Sols = true;
								//std::cout << "**temp ContFrac(10)\n";
								break;
							}
						}
					}
				}

				if (Paso == 1 && type == 4) {
					Tmp = (H1ModCY1 + Z*H2ModCY1) % g_CY1;
					H1ModCY1 = H2ModCY1;
					H2ModCY1 = Tmp;
					Tmp = (K1ModCY1 + Z*K2ModCY1) % g_CY1;
					K1ModCY1 = K2ModCY1;
					K2ModCY1 = Tmp;
				}
			}
			Dp_M = Dp_M1;   
			Dp_P = Dp_P1;  
			if (Co == 0) {
				Co = 1;
			}
			if (type == 5) {
				Tmp = (g_A1 + Z*g_A2) % T;
				g_A1 = g_A2; 
				g_A2 = Tmp;
				Tmp = (g_B1 + Z*g_B2) % T;
				g_B1 = g_B2; 
				g_B2 = Tmp;
			}
			Dp_Mu = -Dp_Mu; //ChangeSign
			if (Paso == 2) {
				if (g_A2 != 0 && Z>(quintillion / 10 - g_A1) / g_A2 ||
					g_B2 != 0 && Z>(quintillion / 10 - g_B1) / g_B2) {
					NextConv(&Bi_H1, &Bi_H2, g_A1, g_A2, g_B1, g_B2);
					NextConv(&Bi_K1, &Bi_K2, g_A1, g_A2, g_B1, g_B2);
					//std::cout << "**temp ContFrac (10A)\n";
					g_A1 = g_B2 = 1;
					g_A2 = g_B1 = 0;
				}
				g_A1 += Z*g_A2;
				g_B1 += Z*g_B2;
				Tmp = g_A1; g_A1 = g_A2; g_A2 = Tmp;   // swap A1 and A2
				Tmp = g_B1; g_B1 = g_B2; g_B2 = Tmp;   // swap B1 and B2
			}
			Conv++;
		} while (Co<0);

		if (!secondDo) {
			continue;    // go to next step (paso = 2)
		}

		Mu = DoublePrecToLong(Dp_Mu);
		L = DoublePrecToLong(Dp_L);
		K = DoublePrecToLong(Dp_K);
		M = DoublePrecToLong(Dp_M);
		P = DoublePrecToLong(Dp_P);

		do {
			P1 = (g_Disc - M*M) / P;    /* P & Q should be > 0 (See Knuth Ex 4.5.3-12) */
			Z = (SqrtDisc + M) / P1;
			M1 = Z*P1 - M;
			if (type == 1 && P == Mu) {    /* Solution found */
				if (Co % 2 == 0 || K != P1 || L != M1) {
					if (Paso == 2) {
						if (g_A2 != 0) {
							NextConv(&Bi_H1, &Bi_H2, g_A1, g_A2, g_B1, g_B2);
							NextConv(&Bi_K1, &Bi_K2, g_A1, g_A2, g_B1, g_B2);
							g_A1 = g_B2 = 1;
							g_A2 = g_B1 = 0;
						}
						ShowLargeXY(X1, Y1, Bi_H1, Bi_K1, true, "NUM(" + numToStr(Conv) + ") = ",
							"DEN(" + numToStr(Conv) + ") = ");
					}
					Sols = true;
				}
				//std::cout << "**temp ContFrac(11)\n";
				break;
			}
			if (type == 3 || type == 4) {
				if ((Co & 1) == 0 && A*g_F>0 && g_Disc == 5) {   /* Solution found */
					if (Paso == 1) {
						Sols = true;
						//std::cout << "**temp ContFrac(12)\n";
						break;
					}
					else {
						//std::cout << "**temp ContFrac(12B)\n";
						NextConv(&Bi_H1, &Bi_H2, g_A1, g_A2, g_B1, g_B2);
						NextConv(&Bi_K1, &Bi_K2, g_A1, g_A2, g_B1, g_B2);
						g_A1 = g_B2 = 1;
						g_A2 = g_B1 = 0;
						Bi_H2 -= Bi_H1;	
						Bi_K2 -= Bi_K1;
						if (ShowHomoSols(type, Bi_H2, Bi_K2, s, T, MagnifY,
							"NUM(" + numToStr(Conv) + ") - NUM(" + numToStr(Conv - 1) + ") = ",
							"DEN(" + numToStr(Conv) + ") - DEN(" + numToStr(Conv - 1) + ") = ")) {
							Sols = true;
							//std::cout << "**temp ContFrac(13) - solution found\n";
							break;
						}
						Bi_H2 += Bi_H1;   
						Bi_K2 += Bi_K1;   
					}
				}
				if (P1 == Mu) {   /* Solution found */
					if (Co % 2 == 0 || K != P1 || L != M1) {
						if (Paso == 2) {
							if (g_A2 != 0) {
								//std::cout << "**temp ContFrac(13B) - solution found\n";
								NextConv(&Bi_H1, &Bi_H2, g_A1, g_A2, g_B1, g_B2);
								NextConv(&Bi_K1, &Bi_K2, g_A1, g_A2, g_B1, g_B2);
								g_A1 = g_B2 = 1;
								g_A2 = g_B1 = 0;
							}
							if (ShowHomoSols(type, Bi_H2, Bi_K2, s, T, MagnifY,
								"NUM(" + numToStr(Conv) + ") = ",
								"DEN(" + numToStr(Conv) + ") = ")) {
								Sols = true;
								//std::cout << "**temp ContFrac(14) - solution found\n";
								break;
							}
						}
						else {
							if (type == 4) {
								Tmp = H2ModCY1*T;
								if ((Tmp - g_CY0) % g_CY1 == 0 || (Tmp + g_CY0) % g_CY1 == 0) {
									Sols = true;
									break;
								}
							}
							else {
								Sols = true;
								//std::cout << "**temp ContFrac(15)\n";
								break;
							}
						}
					}
				}

				if (Paso == 1 && type == 4) {
					Tmp = (H1ModCY1 + Z*H2ModCY1) % g_CY1;
					H1ModCY1 = H2ModCY1;
					H2ModCY1 = Tmp;
					Tmp = (K1ModCY1 + Z*K2ModCY1) % g_CY1;
					K1ModCY1 = K2ModCY1;
					K2ModCY1 = Tmp;
				}
			}

			Co++;
			if (Co % 5000 == 0) {
				std::cout << "Conv: " << Co << " (Eq " << EqNbr << " of " << NbrEqs << ") \n";
				std::cout << NbrSols << " solution" << (NbrSols == 1 ? "\n" : "s\n");
			}
			if (Co % 5000 == 2500) {
				std::cout << NbrSols << " solution" << (NbrSols == 1 ? "\n" : "s\n");
			}

			M = M1;
			P = P1;
			if (type == 2 && P1 == Mu) {
				NextConv(&Bi_H1, &Bi_H2, g_A1, g_A2, g_B1, g_B2);
				NextConv(&Bi_K1, &Bi_K2, g_A1, g_A2, g_B1, g_B2);
				Sols = true;
				//std::cout << "**temp ContFrac(17)\n";
				break;
			}
			if (type == 5) {
				if (P1 == Mu) {
					NbrCo = Co;
					Sols = true;
					//std::cout << "**temp ContFrac(18) NbrCo=" << NbrCo << "\n";
					break;
				}
				else {
					Tmp = (g_A1 + Z*g_A2) % T;
					g_A1 = g_A2; g_A2 = Tmp;
					Tmp = (g_B1 + Z*g_B2) % T;
					g_B1 = g_B2; 
					g_B2 = Tmp;
				}
			}
			Mu = -Mu;
			if (Paso == 2) {
				if (g_A2 != 0 && Z>(quintillion / 10 - g_A1) / g_A2 || g_B2 != 0 && Z>(quintillion / 10 - g_B1) / g_B2) {
					NextConv(&Bi_H1, &Bi_H2, g_A1, g_A2, g_B1, g_B2);
					NextConv(&Bi_K1, &Bi_K2, g_A1, g_A2, g_B1, g_B2);
					//std::cout << "**temp ContFrac(18B)\n";
					g_A1 = g_B2 = 1;
					g_A2 = g_B1 = 0;
				}
				g_A1 += Z*g_A2; 
				g_B1 += Z*g_B2;
				Tmp = g_A1; g_A1 = g_A2; g_A2 = Tmp;   // swap A1 and A2
				Tmp = g_B1; g_B1 = g_B2; g_B2 = Tmp;   // swap B1 and B2
			}
			Conv++;
			/*std::cout << "**temp ContFrac(18A) NbrCo=" << NbrCo << "  Co=" << Co;
			std::cout << "  K=" << K << "  P=" << P << "  L=" << L << "  M=" << M << "\n";*/
		} while (NbrCo>0 ? Co != NbrCo : Co % 2 != 0 || K != P || L != M);

		//std::cout << "**temp ContFrac(18C)\n";
		/* type = 5: Find convergents for x^2 + Bxy + ACy^2 = 1 (mod B^2-4AC) */
		if (type == 5) {
			//std::cout << "**temp ContFrac(19)\n";
			break;   // break out of for loop
		}
	}                        /* end for */

	//std::cout << "**temp ContFrac returns " << Sols << "\n";
	assert(_CrtCheckMemory());
	return Sols;
}

/******************************************************************************
* H = Constant term                                                           *
* T = Divisor of the square part of the constant term                         *
* A = X^2 coefficient                                                         *
* B = XY coefficient                                                          *
* C = Y^2 coefficient                                                         *
* SCFstr = Nothing or apostrophe (complete quad equation)                     *
* called from solveEquation                                                   *
* overwrites global variable g_F, Bi_L1, Bi_L2, NbrSols                       *
******************************************************************************/
void SolContFrac(long long H, long long T, long long A, long long B, long long C,
	std::string SCFstr) {
	long long factor[64] = { 0 };
	long long P[64] = { 0 };
	long long Q[64] = { 0 };
	long long Dif[64] = { 0 };   /* Holds difference */
	long long mod[64] = { 0 };
	long long pos[64] = { 0 };
	long long Tmp, q, s, t, v, Pp, dif, Sol1, Sol2, Modulo;
	long long SqrtDiscCopy = SqrtDisc;
	mpz_int NUMcopy, DENcopy;
	long long DiscCopy = g_Disc;
	long long ValA, ValB, ValC, ValF, ValAM, ValBM, ValCM;
	long long VarD, VarK, VarQ, VarR, VarV, VarW, VarX, VarY, VarY1;
	mpz_int Dp_A, Dp_B, Dp_C, Dp_R, Dp_S, Dp_T;
	mpz_int Bi_Xcopy, Bi_Ycopy;
	int index, index2, cont;
	int NbrFactors;
	long long gcdAF, MagnifY;
	int cuenta = 0;
	long long OrigA, OrigC;
	bool ShowHR = false;

	//std::cout << "**temp SolContFrac:  H=" << H << "  T=" << T << "  A=" << A << "  B=" << B << "  C=" << C << "\n";;

	NUMcopy = g_NUM;    // save global NUM to local 
	DENcopy = g_DEN;    // save global DEN to local 

	g_F = H / T / T;
	if (teach && T>1) {
		std::cout << "Since " << T << " * " << T << " is a divisor of the constant term ("
			<< H << "), the solutions should be " << T << " times the solutions of ";
		ShowEq(A, B, C, 0, 0, g_F, "u", "v");
		printf(" = 0.\n");
		if (abs(g_F) != 1) {
			printf(" Let F be the constant term.");
		}
		putchar('\n');
		UU = "U";
		VU = "V";
		UL = "u";
		VL = "v";
		FP = "F";
	}
	if (teach && T == 1) {
		UU = "X" + SCFstr;
		VU = "Y" + SCFstr;
		UL = "x" + SCFstr;
		VL = "y" + SCFstr;
		FP = "f" + SCFstr;
	}
	gcdAF = gcd(A, g_F);
	OrigA = A;
	OrigC = C;
	if (teach && gcdAF > 1) {
		std::cout << "Since gcd(A,F) = gcd(" << A << "," << g_F << ") = " << gcdAF
			<< " > 1, we have to replace y = ny' where n is a divisor of gcd(A,F).\n";
	}
	for (MagnifY = 1; MagnifY*MagnifY <= gcdAF; MagnifY++) {
		do {
			if (gcdAF / MagnifY*MagnifY != gcdAF) {
				continue;   // if MagnifY^2 is a not factor of gcdAF skip to next value 
			}
			if (teach) {
				if (ShowHR) {
					//w("<HR>"); 
				}
				else {
					ShowHR = true;
				}
			}
			MagnifY = gcdAF / MagnifY;
			g_F = H / T / T / MagnifY;
			ValF = abs(g_F);
			A = OrigA / MagnifY;
			C = OrigC*MagnifY;
			ValA = (A + ValF) % ValF;
			ValB = (B + ValF) % ValF;
			ValC = (C + ValF) % ValF;
			if (teach) {
				if (MagnifY != 1) {
					std::cout << "Let y = " << MagnifY << "y' We obtain: ";
				}
				putchar('\n');
				ShowEq(A, B, C, 0, 0, g_F, "x", "y'");
				printf(" = 0\n");
			}
			/* Find factors of F, store in array factors */
			NbrFactors = 0;
			Tmp = ValF;
			if (Tmp == 1) {
				factor[NbrFactors++] = 1;
			}
			else {
				while ((Tmp % 2) == 0) {
					factor[NbrFactors++] = 2;
					Tmp /= 2;
				}
				while ((Tmp % 3) == 0) {
					factor[NbrFactors++] = 3;
					Tmp /= 3;
				}
				s = 5;        /* Sequence of divisors 5, 7, 11, 13, 17, 19,... */
				do {
					while ((Tmp%s) == 0) {
						factor[NbrFactors++] = s;
						Tmp /= s;
					}
					s += 2;
					while ((Tmp%s) == 0) {
						factor[NbrFactors++] = s;
						Tmp /= s;
					}
					s += 4;
				} while (s*s <= Tmp);
				if (Tmp != 1) {
					factor[NbrFactors++] = Tmp;
				}
			}
			/* complete list of prime factors of F now in array F */

			mod[NbrFactors] = Tmp = 1;
			Pp = (2 * ValA) % ValF;
			for (index = NbrFactors - 1; index >= 0; index--) {
				P[index] = Pp;
				Tmp *= factor[index];
				mod[index] = Tmp;
				Pp = MultMod(MultMod(Pp, factor[index], ValF), factor[index], ValF);
			}
			Modulo = factor[NbrFactors - 1];  // get largest prime factor
			ValAM = (ValA + Modulo) % Modulo;
			ValBM = (ValB + Modulo) % Modulo;
			ValCM = (ValC + Modulo) % Modulo;
			if (ValAM == 0) {  /* Linear equation: sol=-C/B */
				Sol1 = Sol2 = MultMod(Modulo - ValCM, ModInv(ValBM, Modulo), Modulo);
			}
			else {    /* Quadratic equation Ax^2+Bx+C=0 (mod F) */
				if (Modulo>2) {
					Sol1 = MultMod(ValBM, ValBM, Modulo) - MultMod(4 * ValAM, ValCM, Modulo);
					if (Sol1<0) { Sol1 += Modulo; }
					/* Find square root of Sol1 mod Modulo */
					if (Sol1 == 0) {                 /* if double root: sol = -t/2a */
						Sol1 = Sol2 = MultMod(ModInv((2 * ValAM + Modulo) % Modulo, Modulo), ((-ValBM) + Modulo) % Modulo, Modulo);
					}
					else {
						if (ModPow(Sol1, (Modulo - 1) / 2, Modulo) == 1) { /* if sols exist */
							if (Modulo % 8 == 5) {
								s = ModPow(2 * Sol1, (Modulo - 5) / 8, Modulo);
								Sol1 = MultMod(MultMod(MultMod(MultMod(2 * Sol1, s, Modulo), s, Modulo) - 1, Sol1, Modulo), s, Modulo);
							}
							else {
								if (Modulo % 8 != 1) {
									Sol1 = ModPow(Sol1, (Modulo + 1) / 4, Modulo);
								}
								else {
									VarR = 1;
									VarQ = Modulo - 1;
									while (VarQ % 2 == 0) {
										VarQ /= 2;
										VarR *= 2;
									}
									VarX = 2;
									while (true) {
										VarY = ModPow(VarX, VarQ, Modulo);
										if (ModPow(VarY, VarR / 2, Modulo) != 1) { break; }
										VarX++;
									}
									VarX = ModPow(Sol1, (VarQ - 1) / 2, Modulo);
									VarV = MultMod(Sol1, VarX, Modulo);
									VarW = MultMod(VarV, VarX, Modulo);
									while (VarW != 1) {
										VarK = 1; VarD = VarW;
										while (VarD != 1) {
											VarD = MultMod(VarD, VarD, Modulo);
											VarK *= 2;
										}
										VarD = ModPow(VarY, VarR / VarK / 2, Modulo);
										VarY1 = MultMod(VarD, VarD, Modulo);
										VarR = VarK;
										VarV = MultMod(VarV, VarD, Modulo);
										VarW = MultMod(VarW, VarY1, Modulo);
										VarY = VarY1;
									}   /* end while */
									Sol1 = VarV;
								}     /* end modulo 8 = 1 */
							}
							s = ModInv((2 * ValAM) % Modulo, Modulo);
							Sol2 = MultMod((Modulo + Sol1 - ValBM) % Modulo, s, Modulo);
							Sol1 = MultMod((2 * Modulo - Sol1 - ValBM) % Modulo, s, Modulo);
						}
						else {   /* No solution exists */
							Sol1 = Sol2 = -1;
						}
					}
				}
				else {         /* Modulo <= 2 */
					if (Modulo == 2) {
						switch ((int)ValBM * 2 + (int)ValCM) {
						case 0:           /* A = 1, B = 0, C = 0 */
							Sol1 = Sol2 = 0;    /* Solution only for s=0 */
							break;
						case 1:           /* A = 1, B = 0, C = 1 */
							Sol1 = Sol2 = 1;    /* Solution only for s=1 */
							break;
						case 2:           /* A = 1, B = 1, C = 0 */
							Sol1 = 0;         /* Solution for s=0 and s=1 */
							Sol2 = 1;
							break;
						default:          /* A = 1, B = 1, C = 1 */
							Sol1 = Sol2 = -1;   /* No solutions */
							break;
						}
					}                   /* End Modulo = 2 */
					else {                /* Modulo = 1 */
						Sol1 = Sol2 = 0;
					}
				}
			}               /* End Quadratic Equation */

			//std::cout << "**temp SolContFrac: Sol1 =" << Sol1 << "  Sol2=" << Sol2 << "\n";

			ValAM = (ValA + ValF) % ValF;
			ValBM = (ValB + ValF) % ValF;
			ValCM = (ValC + ValF) % ValF;
			if (teach) {
				UL1 = UL;
				VL1 = VL + (MagnifY == 1 ? "" : "'");
				std::cout << "Let " << UL1 << " = s" << VL1 << " - " << FP << "z, so [-(as" << sq << " + bs + c)/"
					<< FP << "]" << VL1 << sq << " + (2as + b)" << VL1 << "z - a" << FP << "z" << sq << " = 1.\nSo \n";
				ShowEq(A, 0, 0, B, 0, C, "s", "");
				std::cout << " should be multiple of " << ValF << "\n";
			}

			NbrEqs = EqNbr = 0;
			int sol = 0;
			t = Sol1;
			/* if Sol1 >= 0 execute loop twice, otherwise not at all */
			for (cont = (Sol1<0 ? 2 : 0); cont<2; cont++) {
				index = NbrFactors - 1;
				v = mod[index];
				dif = 0;
				q = (MultMod((MultMod(ValAM, t, ValF) + ValBM) % ValF, t, ValF) + ValCM) % ValF;  /* q%v = 0 */
				while (true) {
					if (q%v == 0) {
						if (index == 0) {          /* Solution found */
							NbrEqs++;
							if (teach) {
								s = t*mod[1];
								for (index2 = 1; index2<NbrFactors; index2++) {
									s += pos[index2] * mod[index2 + 1];
								}
								s = s%ValF;
								if (sol == 0) {
									sol = 1;
									std::cout << "This holds for s = " << s;
								}
								else {
									std::cout << ", " << s;
								}
							}
							else {
								if (sol == 0)
									sol = 1;
							}
						}
						else {  /* solution not found */
							pos[index] = t;
							t = 0;
							for (index2 = index; index2<NbrFactors; index2++) {
								t += pos[index2] * mod[index2 + 1];
							}
							t = t%ValF;
							Dif[index] = dif;
							Q[index] = q;
							dif = MultMod((MultMod((2 * t + v) % ValF, ValAM, ValF) + ValBM) % ValF, v, ValF);
							Pp = P[--index];
							t = 0; v = mod[index];
							continue;
						}
					}

					if (index != NbrFactors - 1 && ++t < factor[index]) {
						q = (q + dif) % ValF;
						dif = (dif + Pp) % ValF;
						continue;
					}

					else {
						while (++index < NbrFactors) {
							t = pos[index];
							v = mod[index];   /* Restore previous values */
							if (index < NbrFactors - 1 && ++t < factor[index]) {
								Pp = P[index];
								dif = Dif[index];
								q = (Q[index] + dif) % ValF;
								dif = (dif + Pp) % ValF;
								break;  // exit from 'while' loop
							}
						}
						if (index >= NbrFactors) {
							break;
						}
					}
				}

				if (Sol1 == Sol2) {
					break;   /* Do not process double root */
				}
				t = Sol2;
			}
			if (teach) {
				if (sol == 0) {
					printf("No values of s makes the previous assertion true. \n");
					SqrtDisc = SqrtDiscCopy;
					g_NUM = NUMcopy;    
					g_DEN = DENcopy;  
					g_Disc = DiscCopy;
					//std::cout << "**temp SolContFrac exit(1)\n";
					assert(_CrtCheckMemory());
					return;
				}
				printf(". \n");
			}     /* end if (teach) */

			t = Sol1;
			/* if Sol1 >= 0 execute loop twice, otherwise not at all */
			for (cont = (Sol1<0 ? 2 : 0); cont<2; cont++) {
				index = NbrFactors - 1; v = mod[index];
				dif = 0;
				q = (MultMod((MultMod(ValAM, t, ValF) + ValBM) % ValF, t, ValF) + ValCM) % ValF;  /* q%v = 0 */
				while (true) {
					//std::cout << "**temp SolContFrac(2)\n";
					if (q%v == 0) {
						if (index == 0) {          /* Solution found */
							EqNbr++;
							s = t*mod[1];
							for (index2 = 1; index2<NbrFactors; index2++) {
								s += pos[index2] * mod[index2 + 1];
							}
							s = s%ValF;
							//  Calculate Dp_A = -((As+B)s+C)/F
							//            Dp_B = 2As+B
							//            Dp_C = -AF
							//            Dp_R = gcd(Dp_A, Dp_B, Dp_C)

							Dp_C = A * s;             // Dp_C = A*s
							Dp_A = (Dp_C + B) * s;    // Dp_A = (A*s+B)s
							Dp_T = Dp_A + C;          // Dp_T = (A*s+B)s+C
							Dp_R = -g_F;            
							DivideDoublePrecLong(Dp_T, Dp_R, &Dp_A);   // Dp_A = -((As+B)s+C)/F

							Dp_C *= 2;                // Dp_C = 2As
							Dp_B = B;               
							Dp_B += Dp_C;             // Dp_B = 2As+B
							Dp_C = Dp_R * A;          // Dp_C = -AF

							gcd(Dp_A, Dp_B, &Dp_T);   // Dp_T = gcd(Dp_A, Dp_B)
							gcd(Dp_T, Dp_C, &Dp_R);   // Dp_R = gcd(Dp_A, Dp_B, Dp_C)
							/* temporary */
							/*std::cout << "**temp SolContFrac  Dp_A=" << Dp_A << "  Dp_B=" << Dp_B;
							std::cout << "  Dp_C=" << Dp_C << "   DP_R=" << Dp_R << "\n";*/
							/* end temporary */
							if (teach) {
								std::cout << "Let s = " << s << ". Replacing in the above equation:\n";
								ShowBigEq(Dp_A, Dp_B, Dp_C, VL1, "z");
								printf(" = 1\n");
							}

							if (Dp_R > 1) {
								if (teach) {
									std::cout << "Since the gcd of the three coefficients is " << Dp_R <<
										" there are no integer solutions.\n";
								}
							}
							else {
								GetRoot(Dp_A, Dp_B, Dp_C);   // sneaky!! values returned in global vars used by ContFrac
								if (SCFstr == "") {
									if (ContFrac(Dp_A, 3, 1, s, T, MagnifY, A)) {
										if (Bi_L2 < 0) {
											Bi_L1 = -Bi_L1; //ChangeSign
											Bi_L2 = -Bi_L2;//ChangeSign
										}
										//std::cout << "**temp SolContFrac (2C) copy sol to ArrayX and ArrayY";
										Bi_Xcopy = Bi_L1; 
										Bi_Ycopy = Bi_L2; 
									}

									if (ContFrac(Dp_A, 3, (-1), s, T, MagnifY, A)) {
										if (Bi_L2 < 0) {
											Bi_L1 = -Bi_L1;  //ChangeSign;
											Bi_L2 = -Bi_L2;  //ChangeSign
										}
										if (teach) {
											NbrSols++;
											printf("Choosing the solution with minimum y:  \n%d solutions\n", NbrSols);

										}
										//if (Compare(Bi_L2, Bi_Ycopy) > 0) {
										if (Bi_L2>  Bi_Ycopy)  {
											//std::cout << "**temp SolContFrac (2A)\n";
											ShowLargeXY("X", "Y", Bi_Xcopy, Bi_Ycopy, true, "", "");
										}
										else {
											//std::cout << "**temp SolContFrac (2B)\n";
											ShowLargeXY("X", "Y", Bi_L1, Bi_L2, true, "", "");
										}
										//w("</TABLE>");
									}
									else {
										if (Bi_Ycopy != 0) {
											//std::cout << "**temp SolContFrac(3)\n";
											ShowLargeXY("X", "Y", Bi_Xcopy, Bi_Ycopy, true, "", "");
										}
									}
								}
								else {
									ContFrac(Dp_A, 4, 1, s, T, MagnifY, A);
									ContFrac(Dp_A, 4, -1, s, T, MagnifY, A);
								}
							}
						}
						else {
							pos[index] = t;
							t = 0;
							for (index2 = index; index2<NbrFactors; index2++) {
								t += pos[index2] * mod[index2 + 1];
							}
							t = t%ValF;
							Dif[index] = dif;
							Q[index] = q;
							dif = MultMod((MultMod((2 * t + v) % ValF, ValAM, ValF) + ValBM) % ValF, v, ValF);
							Pp = P[--index];
							t = 0;
							v = mod[index];
							continue;
						}
					}
					if (index < NbrFactors - 1 && ++t < factor[index]) {
						q = (q + dif) % ValF;
						dif = (dif + Pp) % ValF;
						//std::cout << "**temp SolContFrac(4)\n";
						continue;
					}
					else {
						while (++index < NbrFactors) {
							//std::cout << "**temp SolContFrac(5)\n";
							t = pos[index]; v = mod[index];   /* Restore previous values */
							if (index < NbrFactors - 1 && ++t < factor[index]) {
								Pp = P[index];
								dif = Dif[index];
								q = (Q[index] + dif) % ValF;
								dif = (dif + Pp) % ValF;
								break;
							}
						}
						if (index >= NbrFactors) {
							//std::cout << "**temp SolContFrac(6)\n";
							break;   // break out of while loop
						}
					}
				}
				//std::cout << "**temp SolContFrac(7)\n";
				if (Sol1 == Sol2) {
					break;  /* Do not process double root */
				}
				//std::cout << "**temp SolContFrac(8)\n";
				t = Sol2;
			}
			if (teach) {
				putchar('\n');
			}
			SqrtDisc = SqrtDiscCopy;
			g_NUM  = NUMcopy;   // restore original values     
			g_DEN  = DENcopy;  
			g_Disc = DiscCopy;
			//std::cout << "**temp SolContFrac(9)\n";
		} while (MagnifY*MagnifY > gcdAF);
		//std::cout << "**temp SolContFrac(10)\n";
	}   // end of for loop

	//std::cout << "**temp SolContFrac exit(11) \n";
	assert(_CrtCheckMemory());
}

/* calculate base^exp mod (Mod)*/
long long ModPow(long long Base, long long Exp, long long Mod) {
	long long Pot, Pwr, mask, value;
	if (Exp == 0) { return 1LL; }

	assert(Exp >= 0);

	mask = 1LL;
	Pot = 1LL;
	Pwr = Base;
	value = 0;
	while (true) {
		if ((Exp & mask) != 0) {
			Pot = MultMod(Pot, Pwr, Mod);
			value += mask;
			if (value == Exp) { return Pot; }
		}
		mask *= 2LL;
		Pwr = MultMod(Pwr, Pwr, Mod);
	}
}

long long ModInv(long long Val, long long Mod) {
	long long U1, U3, V1, V3, Aux, Q;
	U1 = 1; U3 = Val; V1 = 0; V3 = Mod;
	while (V3 != 0) {
		Q = U3 / V3;
		Aux = U1 - V1*Q;
		U1 = V1;
		V1 = Aux;
		Aux = U3 - V3*Q;
		U3 = V3;
		V3 = Aux;
	}
	return (U1 + Mod) % Mod;
}

/* Calculate Num1*Num2 mod Mod.
/* removes overflow risk. sign of result is same as sign of Num1*Num2 */
long long MultMod(long long Num1, long long Num2, long long Mod) {
	mpz_int N1, Prod;
	long long x;

	Mod = abs(Mod);   // ensure Mod is +ve 
	N1 = Num1; 

	Prod = (N1 * Num2)%Mod ;  
	x = DoublePrecToLong(Prod);   // magnitude of Prod  < Mod, therfore cannot overflow.
	return x;
}