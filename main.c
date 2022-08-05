#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include "matriz.h"
#include "dados.h"
#include "LU.h"

#define EPLISON1 0.00001
#define EPLISON2 0.00001
// f == function

int main(int argc, char **argv)
{

    //double tLU;
    //double tSL;
    //double tR;
    // inicializa gerador de números aleatóreos
    srand(202201);
    args *argumentos = (args *)malloc(sizeof(args *));

    /*
    SL *s, *n;
    double *r, *t;
    s = aloca_sist();

    preenche_sis(s);
    r = aloca_vetor(s->n);
    pregauss(s);
    GaussSeidel(s, EPLISON1, EPLISON2);

    r = residuo(s);
    refinamento(s, r, EPLISON1, EPLISON2);
    lee_vetor(r, s);
    n = copia_sl(s, r);
    inicializa_vetorx(s);
    */

    trata_args(argc, argv, argumentos);
    refLU(argumentos);
}