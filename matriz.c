#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include "matriz.h"
#include "time.h"
#include "dados.h"
#include "math.h"

void retrossubs(SL *sis)
{
    for (int i = (sis->n) - 1; i >= 0; --i)
    {
        sis->X[i] = sis->B[i];
        for (int j = i + 1; j < sis->n; j++)
        {
            sis->X[i] -= sis->A[i][j] * sis->X[j];
        }
        sis->X[i] /= sis->A[i][i];
    }
}

void trocalinha(SL *sis, int linha_1, int linha_n)
{
    double aux;
    for (int k = 0; k < sis->n; k++)
    {
        aux = sis->A[linha_1][k];
        sis->A[linha_1][k] = sis->A[linha_n][k];
        sis->A[linha_n][k] = aux;
    }
}

void pivoteamento(int colun_inicial, SL *sis)
{
    int maior_i = colun_inicial;
    // printf("maior i antes  %d \n ",maior_i);
    for (int i = colun_inicial; i < sis->n; i++)
    {
        if (sis->A[i][colun_inicial] > sis->A[maior_i][colun_inicial])
        {
            maior_i = i;
        }
    }
    trocalinha(sis, colun_inicial, maior_i);
}

void eliminacaoGauss(SL *sis)
{

    for (int i = 0; i < sis->n; i++)
    {
        pivoteamento(i, sis);
        for (int j = i + 1; j < sis->n; j++)
        {
            double m = sis->A[j][i] / sis->A[i][i];
            sis->A[j][i] = 0.0;
            for (int k = i + 1; k < sis->n; k++)
            {
                sis->A[j][k] -= sis->A[i][k] * m;
            }
            sis->B[j] -= sis->B[i] * m;
        }
    }
}

SL *pregauss(SL *sis)
{
    double aux;

    for (int i = 0; i < sis->n; i++)
    {
        for (int j = 0; j < sis->n; j++)
        {
            if (i == j)
            {
                aux = sis->A[i][j];
                sis->A[i][j] = sis->B[i];
                sis->B[i] = aux;
            }
            else
            {
                sis->A[i][j] = -sis->A[i][j];
            }
        }
    }
    lee_matriz(sis->A, sis->n);
    lee_vetor(sis->B, sis->n);
    return sis;
}
//devia devolver int mas nao devolve
void GaussSeidel(SL *sis, double eplison1, double eplison2)
{
    //double erro_maior = 0.0;
    pregauss(sis);
    double erro_temp = 0.0;
    double aux;
    int k = 0;

    while ((k < 20) && (soma_residuo(residuo(sis), sis->n) < eplison1))
    {
        for (int i = 0; i < sis->n; i++)
        { // lembrar de dividir o coeficiente
            aux = sis->X[i];
            for (int j = 0; j < sis->n; j++)
            {
                sis->X[i] = +sis->A[i][j] * sis->X[j] / aux; // Troca de linha  isolando a diagonal e divide pelo numerador de x
                if (fabs(sis->X[i]) > erro_temp)
                {
                    erro_temp = sis->X[i];
                }
            }
        }

        k++;
        printf("Timestamp: %d\n", (int)time(NULL));
    }
}

double *sub_vetor(double *v1, double *v2, int tam)
{ // subtrai vetor
    double *v_resul = aloca_vetor(tam);
    for (int i = 0; i < tam; i++)
    {
        v_resul[i] = v1[i] - v2[i];
    }
    return v_resul;
} // ax-b

double *residuo(SL *sis)
{
    double *v_AX = aloca_vetor(sis->n);
    double *v_resul = aloca_vetor(sis->n);
    for (int i = 0; i < sis->n; i++)
    {
        double soma = 0.0;
        for (int j = i; j < sis->n; j++)
        {
            soma += sis->A[i][j] * sis->X[j];
        }
        v_AX[i] = soma;
    }
    v_resul = sub_vetor(sis->B, v_AX, sis->n);
    return v_resul;
}

double criterio_parada1(double *residuo, int tam)
{
    double norma = 0.0;
    for (int i = 0; i < tam; i++)
    {
        norma += fabs(residuo[i]);
    }
    return norma;
}

double criterio_parada2(double *r_ant, double *r_dps, int tam)
{ // compara o x anterior com o atual
    double maior_diff = 0.0;
    for (int i = 0; i < tam; i++)
    {
        if (fabs(r_ant[i] - r_dps[i]) > maior_diff)
        {
            maior_diff = fabs(r_ant[i] - r_dps[i]);
        }
    }
    return maior_diff;
}
double ** soma_matriz(double ** mA, double ** mB, int tam)
{ // soma matriz
    double **m_result = aloca_matriz(tam);
    for (int i = 0; i < tam; i++)
    {
        for( int j = 0 ; j< tam; j++){
        m_result[i][j] = mA[i][j] + mB[i][j];
        }
    }
    return m_result;
} // ax-b
double *soma_vetor(double *v1, double *v2, int tam)
{ // soma o x anterior ao w  x(i)+w
    double *v_result = aloca_vetor(tam);
    for (int i = 0; i < tam; i++)
    {
        v_result[i] = v1[i] + v2[i];
    }
    return v_result;
}

double soma_residuo(double *v1, int tam)
{ // soma o x anterior ao w  x(i)+w
    double soma = 0.0;
    for (int i = 0; i < tam; i++)
    {
        soma += v1[i];
    }
    return soma;
}


double* refinamento(SL *sis, double *vetor, double eplison1, double eplison2)
{
    SL *novo_sis;
    //double *v_result = aloca_vetor(sis->n); //  SOMA VETOR W OBTIDO A X(i)
    if (criterio_parada1((residuo(sis)), sis->n) < eplison1)
    {
        return sis->X;
    }
    else
    {
        novo_sis = copia_sl(sis, vetor); //  matriz A  X vetor (w) = Vetor (Residuo)
        eliminacaoGauss(novo_sis);
        retrossubs(novo_sis); // obtem o valor do novo X(W)
        soma_vetor(sis->X, novo_sis->X, sis->n);
        if (criterio_parada2(sis->X, novo_sis->X, sis->n) < eplison2)
        {
            return novo_sis->X;
        }
        else
        {
            refinamento(novo_sis, residuo(novo_sis), eplison1, eplison2);
        }
    }
    return novo_sis->X;
    printf("Timestamp: %d\n", (int)time(NULL));
}
