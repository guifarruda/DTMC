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

#ifndef MonteCarlo_H
#define MonteCarlo_H

#define IGNORANT 0
#define SPREADER 1
#define STIFLER 2

#include <stdlib.h>
#include <igraph/igraph.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

typedef struct {
    long *Ignorant, *Spreader, *Stifler;
}rates;

void FreeRate(rates *R);

rates* MonteCarlo ( igraph_t *graph, long steps, double lambda, double alpha, double delta1, double delta2, double gamma, double eta, double beta, long ncont, igraph_vector_t* state, gsl_rng * randomic );

#endif
