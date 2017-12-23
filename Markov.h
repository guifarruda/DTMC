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

#ifndef Markov_H
#define Markov_H

#include <stdlib.h>
#include <igraph/igraph.h>

int ContactProbabilitiesRP ( long N, long t, igraph_vector_t *adjl, double lambda, double alpha, double delta1, double delta2, double gamma, double eta, double beta, double **P_ignorant, double **P_spreader, double **P_stifler, double *a, double *b, double *c );

int ContactProbabilitiesCP ( long N, long t, igraph_vector_t *adjl, double lambda, double alpha, double delta1, double delta2, double gamma, double eta, double beta, double **P_ignorant, double **P_spreader, double **P_stifler, double *a, double *b, double *c, igraph_vector_t *k );

int TotalProbability ( long N, long t, igraph_vector_t *adjl, double lambda, double alpha, double delta1, double delta2, double gamma, double eta, double beta, double **P_ignorant, double **P_spreader, double **P_stifler, double *a, double *b, double *c );

#endif
