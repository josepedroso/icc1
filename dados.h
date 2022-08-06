#ifndef __DADOS__
#define __DADOS__

#include <stdio.h>
#include <unistd.h>

/**
 * @brief 
 * 
 */
typedef struct args{
    FILE *IN;
    FILE *OUT;
    int N;
    int K;
}args;

void trata_args(int argc, char **argv, args *argumentos);

double *aloca_vetor(int n);

double **aloca_matriz(int n);

void lee_matriz(double **m, int n, FILE *out);

void lee_vetor(double *v, int n);

SL *aloca_sist(int n);

void preenche_sis(SL *sis);

void lee_sis(SL *sis, FILE *out);

double **copia_matriz(double **matriz_original, int tam);

double *copia_vetor(double *vetor, int tam);

double *inicializa_vetorx(SL *sis);

SL *copia_sl(SL *sis, double *vetor);

#endif