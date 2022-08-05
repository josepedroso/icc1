#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include "matriz.h"
#include "LU.h"
#include <math.h>

#define EPLISON1 0.001
#define EPLISON2 0.001
#define c_max 10


void iniSisLU(LU * sis, double coef_max)
{
    unsigned int n = sis->n;
    // para gerar valores no intervalo [0,coef_max]
    double invRandMax = ((double)coef_max / (double)RAND_MAX);

    // inicializa sistema normal 
    // inicializa a matriz A
    for (unsigned int i = 0; i < n; ++i)
    {
        for (unsigned int j = 0; j < n; ++j)
        {
            sis->L[i][j] = (double)rand() * invRandMax;
        }
        sis->X[i] = 0.0;
        sis->L[i][i] = 1;
    }
}

LU *aloca_LU(args *argumentos)
{   
    int tam;
    argumentos->N ? tam = argumentos->N : fscanf(argumentos->IN, "%d ", &tam);
    LU *sist = malloc(tam * tam * sizeof(LU *));
    sist->n = tam;

    sist->L = aloca_matriz(sist->n);
    sist->U = aloca_matriz(sist->n);
    sist->I = aloca_matriz(sist->n);
    sist->X = aloca_vetor(sist->n);

    return sist;
}
double **matriz_inicial(int tam)
{
    double **m = aloca_matriz(tam);
    for (int i = 0; i < tam; i++)
    {
        for (int j = 0; j < tam; j++)
        {
            if (i == j)
            {
                m[i][j] = 1.0;
            }
            else
            {
                m[i][j] = 0.0;
            }
        }
    }
    return m;
}
void preenche_mId(LU *sis, double **m)
{
    sis->I = copia_matriz(m, sis->n);
}

void preenche_LU_Inicial(LU *sis, double **m, args *argumentos)
{
    if(argumentos->N){
        iniSisLU(sis,c_max);
        preenche_mId(sis, m);
        return;
    }
    for (int i = 0; i < sis->n; i++)
    {
        for (int j = 0; j < sis->n; j++)
        {
            fscanf(argumentos->IN, "%lf ", &(sis->U[i][j])); 
        }
        sis->X[i] = 0.0;
        sis->L[i][i] = 1;
    }
    preenche_mId(sis, m);
}

void preenche_LU(LU *sis, double **mI, double **A)
{

    for (int i = 0; i < sis->n; i++)
    {
        for (int j = 0; j < sis->n; j++)
        {
            if (i == j)
            {
                sis->L[i][j] = 1;
            }
            // printf("coluna %d  linha %d \n",j,i);
            A[i][j] = sis->U[i][j];
        }
        // printf("vetor linha %d \n",i);
        sis->X[i] = 0.0;
    }
    preenche_mId(sis, mI);
}
void trocalinhaLU(LU *sis, int linha_1, int linha_n)
{
    double aux;
    for (int k = 0; k < sis->n; k++)
    {
        aux = sis->U[linha_1][k];
        sis->U[linha_1][k] = sis->U[linha_n][k];
        sis->U[linha_n][k] = aux;
    }
}
void pivoteamentoLU(int colun_inicial, LU *sis)
{
    int maior_i = colun_inicial;
    // printf("maior i antes  %d \n ",maior_i);
    for (int i = colun_inicial; i < sis->n; i++)
    {
        if (sis->U[i][colun_inicial] > sis->U[maior_i][colun_inicial])
        {
            maior_i = i;
            // printf("maior i depois  %d \n ",maior_i);
        }
    }
    trocalinhaLU(sis, colun_inicial, maior_i);
}
void FatoracaoLU(LU *sis)
{

    for (int i = 0; i < sis->n; i++)
    {
        pivoteamentoLU(i, sis);
        for (int j = i + 1; j < sis->n; j++)
        {
            double m = sis->U[j][i] / sis->U[i][i];
            sis->U[j][i] = 0.0;
            sis->L[j][i] = m;
            for (int k = i + 1; k < sis->n; k++)
            {
                sis->U[j][k] -= sis->U[i][k] * m;
            }
        }
    }
}
void convert_M_to_V(double **Iden, double *v, int coluna, int tam)
{
    for (int i = 0; i < tam; i++)
    {
        v[i] = Iden[i][coluna];
    }
}

void convert_LU_SL(LU *lu, SL *sis, double **m, double *v)
{ // copia o sistema seja L ou U para sistema linear e seus vetores
    sis->n = lu->n;
    sis->A = copia_matriz(m, sis->n);
    sis->B = copia_vetor(v, sis->n);
    // lee_vetor(sis->X,sis->n);
    for (int i = 0; i < lu->n; i++)
    {
        sis->X[i] = 0.0;
    }
}
void retrossubsINV(SL *sis)
{
    // linha
    for (int i = 0; i < sis->n; i++)
    {
        sis->X[i] = sis->B[i];
        //  printf("sis->X[i] = sis->B[i]; %lf \n", sis->B[i]);
        for (int j = 0; j < i; j++)
        {
            //  printf("sis->A[i][j] %lf \n", sis->A[i][j]);
            // printf("sis->X[j] %lf \n", sis->X[j]);

            sis->X[i] -= sis->A[i][j] * sis->X[j];
            // printf("sis->X[i] %lf \n", sis->X[i]);
        }
        sis->X[i] /= sis->A[i][i];
        //  printf("sis->X[i] %lf \n", sis->X[i]);
    }
}

double **resolveLU(LU *lu)
{
    double **matriz_Inv = aloca_matriz(lu->n);
    SL *sisL = aloca_sist(lu->n);
    SL *sisU = aloca_sist(lu->n);
    double *v = aloca_vetor(lu->n);
    //double *r = aloca_vetor(lu->n);
    for (int colunaI = 0; colunaI < lu->n; colunaI++)
    { // coluna da matriz identidade
        // printf("convert m to v \n");
        convert_M_to_V(lu->I, v, colunaI, lu->n);
        convert_LU_SL(lu, sisL, lu->L, v);

        // printf("entrou");
        eliminacaoGauss(sisL);
        retrossubs(sisL);
        convert_LU_SL(lu, sisU, lu->U, sisL->X);

        // printf("SISu \n");
        // lee_sis(sisU);
        eliminacaoGauss(sisU);
        // printf("retro \n");
        retrossubs(sisU);
        // printf("SISu \n");
        // lee_sis(sisU);
        matriz_Inversa(sisU, colunaI, matriz_Inv);
    }
    return matriz_Inv;
}
void matriz_Inversa(SL *sisU, int col, double **matriz)
{

    for (int i = 0; i < sisU->n; i++)
    {
        matriz[i][col] = sisU->X[i];
    }
}
double **multiplica_matriz(double **m, double **mInv, int tam)
{
    double **result = aloca_matriz(tam);

    for (int i = 0; i < tam; i++)
    {
        for (int j = 0; j < tam; j++)
        {
            result[i][j] = 0;
            for (int k = 0; k < tam; k++)
            {
                //printf("%lf \n", mInv[i][j]);

                result[i][j] += m[i][k] * mInv[k][j];
            }
        }
    }
    return result;
}
double Norma_LU(double **m, int tam)
{
    double soma = 0;
    for (int i = 0; i < tam; i++)
    {
        for (int j = 0; j < tam; j++)
        {
            //printf("%lf \n", m[i][j]);
             soma = soma + (m[i][j] * m[i][j]);
        }
        //printf("%lf \n", soma);
    }
    soma = sqrt(soma);
    return soma;
}
void refLU(args *argumentos)
{
    int it = 0;
    LU *lu = aloca_LU(argumentos);
    double **matriz = aloca_matriz(lu->n);
    double **matriz_Inv = aloca_matriz(lu->n);
    double **result = aloca_matriz(lu->n);
    double **Identidade = aloca_matriz(lu->n);
    double norma = 1;

    Identidade = matriz_inicial(lu->n);
    preenche_LU_Inicial(lu, Identidade, argumentos);
    matriz = copia_matriz(lu->U, lu->n);
    FatoracaoLU(lu);
    matriz_Inv = resolveLU(lu);
    lee_matriz(matriz_Inv, lu->n);
    result = multiplica_matriz(matriz, matriz_Inv, lu->n);
    norma = Norma_LU(result, lu->n);
    printf(" %lf ", norma);
    it = it + 1;
    printf(" %d", it);
    while ((norma > EPLISON1) && (it < 10))
    {
        preenche_LU(lu, result, Identidade);
        printf(" SEGUNDA VEZ");
        matriz_Inv = resolveLU(lu);
        result = multiplica_matriz(matriz, matriz_Inv, lu->n);
        norma = Norma_LU(result, lu->n);
        it = it + 1;
        printf(" %d", it);
    }
}

/**
 * @brief 
 * 
 * @param sis 
 */
void print_saida(LU *sis){
    printf("# Tempo LU: ");
    printf("# Tempo iter: ");
    printf("# Tempo residuo: ");
    lee_matriz(sis->I, sis->n);
}