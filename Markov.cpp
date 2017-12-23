//     DTMC implements the paper A general Markov chain approach for disease 
//     and rumour spreading in complex networks, Guilherme Ferraz de Arruda, 
//     Francisco Aparecido Rodrigues, Pablo Martín Rodríguez, Emanuele Cozzo,
//     Yamir Moreno, Journal of Complex Networks, cnx024, 
//     https://doi.org/10.1093/comnet/cnx024
//     Copyright (C) 2017  Guilherme Ferraz de Arruda
//     
//     This program is free software; you can redistribute it and/or modify
//     it under the terms of the GNU General Public License as published by
//     the Free Software Foundation; either version 2 of the License, or
//     (at your option) any later version.
//     This program is distributed in the hope that it will be useful,
//     but WITHOUT ANY WARRANTY; without even the implied warranty of
//     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//     GNU General Public License for more details.
//     You should have received a copy of the GNU General Public License
//     along with this program; if not, write to the Free Software
//     Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

#include "Markov.h"

int ContactProbabilitiesRP ( long N, long t, igraph_vector_t *adjl, double lambda, double alpha, double delta1, double delta2, double gamma, double eta, double beta, double **P_ignorant, double **P_in, double **P_stiflerifler, double *a, double *b, double *c )
{
    #pragma omp parallel for
    for ( long i=0; i<N; i++ ) {
        long j = -1;
        for ( long ii=0; ii<igraph_vector_size ( &adjl[i] ); ii++ ) {
            j = ( int ) VECTOR ( adjl[i] ) [ii];
            a[i] = a[i] * ( 1 - ( 1 - delta1 - delta2 ) * lambda * ( P_in[j][t] ) ); // RP
            b[i] = b[i] * ( 1 - alpha * ( P_stiflerifler[j][t]+P_in[j][t] ) );			// RP
            c[i] = c[i] * ( 1 - beta * ( P_ignorant[j][t] ) );					// RP
        }
    }


    return 1;
}


int ContactProbabilitiesCP ( long N, long t, igraph_vector_t *adjl, double lambda, double alpha, double delta1, double delta2, double gamma, double eta, double beta, double **P_ignorant, double **P_in, double **P_stiflerifler, double *a, double *b, double *c, igraph_vector_t *k )
{
    #pragma omp parallel for
    for ( long i=0; i<N; i++ ) {
        long j = -1;
        for ( long ii=0; ii<igraph_vector_size ( &adjl[i] ); ii++ ) {
            j = ( int ) VECTOR ( adjl[i] ) [ii];
            double p = 1.0/ ( VECTOR ( *k ) [i] );
            a[i] = a[i] * ( 1 - ( 1 - delta1 - delta2 ) * lambda * 1.0/ ( VECTOR ( *k ) [j] ) * ( P_in[j][t] ) ); // CP
            b[i] = b[i] + ( alpha * p * ( P_stiflerifler[j][t]+P_in[j][t] ) );				// CP
            c[i] = c[i] + ( beta * p * ( P_ignorant[j][t] ) );					// CP

        }
        b[i] = 1.0 - b[i];		// CP
        c[i] = 1.0 - c[i];		// CP
        
    }


    return 1;
}

int TotalProbability ( long N, long t, igraph_vector_t *adjl, double lambda, double alpha, double delta1, double delta2, double gamma, double eta, double beta, double **P_ignorant, double **P_in, double **P_stiflerifler, double *a, double *b, double *c  )
{


    #pragma omp parallel for
    for ( long i=0; i<N; i++ ) {
        P_ignorant[i][t+1] = P_ignorant[i][t] * a[i] + P_in[i][t] * delta1 + P_stiflerifler[i][t] * gamma;
        P_in[i][t+1] = P_ignorant[i][t] * eta * ( 1.0 - a[i] ) + P_in[i][t] * b[i] * ( 1.0 - delta1 - delta2 ) + P_stiflerifler[i][t] * ( 1.0 - c[i] ) * ( 1.0 - gamma );
        P_stiflerifler[i][t+1] = P_ignorant[i][t] * ( 1.0 - eta ) * ( 1.0 - a[i] ) + P_in[i][t] * ( ( 1.0 -b[i] ) * ( 1.0 - delta1 - delta2 ) + delta2 ) + P_stiflerifler[i][t] * c[i] * ( 1.0 - gamma );
    }

    return 1;
}
