#ifndef __MATRIZ__
#define __MATRIZ__

typedef struct Sist_Linear{
    double **A;
    double *B;
    double *X;
    int n;
}SL;

void retrossubs(SL * sis);

void trocalinha(SL * sis,int j,int i);

void pivoteamento(int colun_inicial,SL * sis);

void eliminacaoGauss(SL * sis);

void preenche_vetorX( SL * sis);

int GaussSeidel(SL * sis,double eplison1,double eplison2);

SL *  pregauss(SL * sis);

double * sub_vetor (double * v1,double *v2,int tam);

double *  residuo(SL * sis);

SL *copia_sl(SL * sis,double * vetor);

double criterio_parada1(double * residuo,int tam);

double criterio_parada2(double * r_ant,double * r_dps,int tam);

double * soma_vetor(double * v1,double * v2,int tam);

double soma_residuo(double * v1,int tam);

double * refinamento (SL *sis, double * vetor,double eplison1,double eplison2);

#endif