#ifndef __LU__
#define __LU__
#include "matriz.h"
#include "dados.h"

typedef struct Sist_LU
{
    double **L;
    double **U;
    double **I;
    double *X;
    int n;
} LU;

void matriz_residuo(double **result, double **m, double **mInv, int tam);

void print_saida(LU *sis);

void iniSisLU(LU * sis, double coef_max);

LU *aloca_LU();

void preenche_mId(LU *sis, double **m);

void preenche_LU(LU *sis, double **mI, double **A);

void trocalinhaLU(LU *sis, int linha_1, int linha_n);

void pivoteamentoLU(int colun_inicial, LU *sis);

void FatoracaoLU(LU *sis);

void convert_M_to_V(double **Iden, double *v, int coluna, int tam);

void convert_LU_SL(LU *lu, SL *sis, double **m, double *v);

void retrossubsINV(SL *sis);

double **resolveLU(LU *lu);

void matriz_Inversa(SL *sisU, int col, double **matriz);

double **multiplica_matriz(double **m, double **mInv, int tam);

double Norma_LU(double **m, int tam, int it, FILE *out);

void refLU(args *argumentos);

void print_tempo(FILE *out);

#endif