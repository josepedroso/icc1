#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include "matriz.h"
#include "LU.h"
#include <math.h>

#define EPLISON1 0.001
#define EPLISON2 0.001
#define c_max 10

void iniSisLU(LU *sis, double coef_max)
{
    int n = sis->n;
    // para gerar valores no intervalo [0,coef_max]
    double invRandMax = ((double)coef_max / (double)RAND_MAX);

    // inicializa sistema normal
    // inicializa a matriz A
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            sis->U[i][j] = (double)rand() * invRandMax;
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
    if (argumentos->N)
    {
        iniSisLU(sis, c_max);
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
            // fprintf(argumentos->OUT, "coluna %d  linha %d \n",j,i);
            A[i][j] = sis->U[i][j];
        }
        // fprintf(argumentos->OUT, "vetor linha %d \n",i);
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
    // fprintf(argumentos->OUT, "maior i antes  %d \n ",maior_i);
    for (int i = colun_inicial; i < sis->n; i++)
    {
        if (sis->U[i][colun_inicial] > sis->U[maior_i][colun_inicial])
        {
            maior_i = i;
            // fprintf(argumentos->OUT, "maior i depois  %d \n ",maior_i);
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
void  convert_V_to_M(double**m,double *v,int coluna,int tam){
    for (int i = 0; i < tam;i++){
        m[i][coluna] = v[i];
    }

}
double ** resolve_sisL(LU * lu,args *argumentos){
    SL *sis = aloca_sist(lu->n);
    double *v = aloca_vetor(lu->n);
    double **m = aloca_matriz(lu->n);
    for (int colunaI = 0; colunaI < lu->n; colunaI++){
        convert_M_to_V(lu->I, v, colunaI, lu->n);
        convert_LU_SL(lu, sis, lu->L, v);
         fprintf(argumentos->OUT, " antes de resolver \n");
        lee_sis(sis, argumentos->OUT);
        eliminacaoGauss(sis);
        retrossubs(sis);
         fprintf(argumentos->OUT, " dps de resolver \n");
        lee_sis(sis, argumentos->OUT);
        convert_V_to_M(m, sis->X, colunaI, lu->n);
        
    }
    lee_matriz(m,lu->n,argumentos->OUT);
    return m;
}
double ** resolve_sisU(LU * lu,args *argumentos){
    SL *sis = aloca_sist(lu->n);
    double *v = aloca_vetor(lu->n);
    double **m = aloca_matriz(lu->n);
    for (int colunaI = 0; colunaI < lu->n; colunaI++){
        convert_M_to_V(lu->I, v, colunaI, lu->n);
        convert_LU_SL(lu, sis, lu->U, v);
         fprintf(argumentos->OUT, " antes de resolver \n");
        lee_sis(sis, argumentos->OUT);
        eliminacaoGauss(sis);
        retrossubs(sis);
         fprintf(argumentos->OUT, " dps de resolver \n");
        lee_sis(sis, argumentos->OUT);
        convert_V_to_M(m, sis->X, colunaI, lu->n);
        
    }
    lee_matriz(m,lu->n,argumentos->OUT);
    return m;
}
double **resolveLU(LU *lu,args *argumentos)
{
    double *v = aloca_vetor(lu->n);
    // double *r = aloca_vetor(lu->n);
    fprintf(argumentos->OUT, "sisL \n");
    lu->I = copia_matriz(resolve_sisL(lu,argumentos), lu->n);
    fprintf(argumentos->OUT, "sisU \n");
    double ** matriz_Inv = resolve_sisU(lu,argumentos);
    return matriz_Inv;
    fprintf(argumentos->OUT, "\n //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// ");
}

void matriz_residuo(double **result, double **m, double **mInv, int tam)
{

    for (int i = 0; i < tam; i++)
    {
        for (int j = 0; j < tam; j++)
        {
            result[i][j] = 0;
            for (int k = 0; k < tam; k++)
            {
                result[i][j] += m[i][k] * mInv[k][j];
            }
            result[i][j] = (i == j ? 1.0 - (result[i][j]) : - (result[i][j]));// ao fazer assim dispensa a identidade
        }
    }
}
void matriz_mult(double **result, double **m, double **mInv, int tam)
{

    for (int i = 0; i < tam; i++)
    {
        for (int j = 0; j < tam; j++)
        {
            result[i][j] = 0;
            for (int k = 0; k < tam; k++)
            {
                result[i][j] += m[i][k] * mInv[k][j];
            }
        }
    }
}
double Norma_LU(double **m, int tam, int it, FILE *out)
{
    double soma = 0;
    for (int i = 0; i < tam; i++)
    {
        for (int j = 0; j < tam; j++)
        {
            // fprintf(argumentos->OUT, "%lf \n", m[i][j]);
            soma = soma + (m[i][j] * m[i][j]);
        }
        // fprintf(argumentos->OUT, "%lf \n", soma);
    }
    soma = sqrt(soma);
    fprintf(out, "# iter %d: <%.15g>\n", it, soma);
    return soma;
}
void refLU(args *argumentos)
{

    // double tLU;
    // double tSL;
    // double tR;

    int it = 0;

    LU *lu = aloca_LU(argumentos);

    double **matriz = aloca_matriz(lu->n);
    double **matriz_Inv = aloca_matriz(lu->n);
    double **result = aloca_matriz(lu->n);
    double **Identidade = aloca_matriz(lu->n);
    double **matriz_Inv_Ant = aloca_matriz(lu->n);

    double norma = 1;

    Identidade = matriz_inicial(lu->n); // matriz identidade inicial
    preenche_LU_Inicial(lu, Identidade, argumentos); // preenche matriz U antes do pivoteamento
    lee_matriz(lu->U, lu->n, argumentos->OUT);

    matriz = copia_matriz(lu->U, lu->n); // salva matriz de entrada
    FatoracaoLU(lu); // converte a matriz L e U EM DOIS SISTEMAS LINEARES
    matriz_Inv = resolveLU(lu, argumentos);
    /*matriz_residuo(result, matriz, matriz_Inv, lu->n);
    lee_matriz(matriz_Inv, lu->n, argumentos->OUT);

    fprintf(argumentos->OUT, "#\n");
    norma = Norma_LU(result, lu->n, 0, argumentos->OUT);
    it = it + 1;
    */
    while (it < argumentos->K)
    {
            if (norma < EPLISON1){
                break;
            }
        matriz_residuo(result, matriz, matriz_Inv, lu->n);
        norma = Norma_LU(result, lu->n, it, argumentos->OUT);
        it = it + 1;
        // matriz_Inv_Ant = copia_matriz(matriz_Inv, lu->n);
        matriz_mult(result,matriz, matriz_Inv, lu->n);        
        preenche_LU(lu, result, Identidade);
        matriz_Inv = resolveLU(lu, argumentos);
       
        //soma_matriz(matriz_Inv, matriz_Inv_Ant, lu->n);
        //lee_matriz(matriz_Inv, lu->n, argumentos->OUT);

        
    }
    print_tempo(argumentos->OUT);
    lee_matriz(matriz_Inv, lu->n, argumentos->OUT);
}

/**
 * @brief
 *
 * @param sis
 */
void print_tempo(FILE *out)
{
    fprintf(out, "# Tempo LU: \n");
    fprintf(out, "# Tempo iter: \n");
    fprintf(out, "# Tempo residuo: \n");
}