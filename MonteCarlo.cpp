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

#include "MonteCarlo.h"


void FreeRate ( rates *R )
{
    free ( R->Ignorant );
    free ( R->Spreader );
    free ( R->Stifler );
    
    free(R);
}


rates* MonteCarlo ( igraph_t *graph, long steps, double lambda, double alpha, double delta1, double delta2, double gamma, double eta, double beta, long ncont, igraph_vector_t* state, gsl_rng * randomic ) {
    long N = igraph_vcount ( graph );

    rates *population = ( rates* ) calloc ( 1, sizeof ( rates ) );
    population->Ignorant = ( long* ) calloc ( steps, sizeof ( long ) );
    population->Spreader = ( long* ) calloc ( steps, sizeof ( long ) );
    population->Stifler = ( long* ) calloc ( steps, sizeof ( long ) );

    igraph_vector_t vertex;
    igraph_vector_init ( &vertex, N );

    for ( long i=0; i<N; i++ ) {
        if ( VECTOR ( *state ) [i] == IGNORANT ) {
            population->Ignorant[0]++;
        } else if ( VECTOR ( *state ) [i] == SPREADER ) {
            population->Spreader[0]++;
        } else if ( VECTOR ( *state ) [i] == STIFLER ) {
            population->Stifler[0]++;
        }

        VECTOR ( vertex ) [i] = i;
    }

    SETVANV ( graph,"S", state );

    long t;
    for ( t=1; t<steps; t++ ) {
        bool change_state = false;

        population->Ignorant[t] = population->Ignorant[t-1];
        population->Spreader[t] = population->Spreader[t-1];
        population->Stifler[t] = population->Stifler[t-1];

        igraph_vector_shuffle ( &vertex );

        if ( population->Spreader[t] > 0 || population->Stifler[t] > 0 ) {
            for ( long v=0; v<N; v++ ) {

                double rnd = gsl_rng_uniform ( randomic );
                long i = VECTOR ( vertex ) [v];
                long status_i = ( long ) VAN ( graph, "S", ( igraph_integer_t ) i );

                if ( status_i == SPREADER ) {

                    if ( rnd < delta1 ) {
                        VECTOR ( *state ) [i] = IGNORANT;
                        population->Spreader[t]--;
                        population->Ignorant[t]++;
                        change_state = true;

                    } else if ( rnd < delta1+delta2 ) {
                        VECTOR ( *state ) [i] = STIFLER;
                        population->Spreader[t]--;
                        population->Stifler[t]++;
                        change_state = true;

                    } else {
                        igraph_vector_t neis;
                        igraph_vector_init ( &neis, 0 );
                        igraph_neighbors ( graph, &neis, ( igraph_integer_t ) i,  IGRAPH_OUT );

                        if ( igraph_vector_size ( &neis ) == 0 ) {
                            igraph_vector_destroy ( &neis );
                            continue;
                        }

                        igraph_vector_shuffle ( &neis );

                        long ncontacts;
                        if ( ncont >= ( long ) igraph_vector_size ( &neis ) ) {
                            ncontacts = igraph_vector_size ( &neis );
                        } else {
                            ncontacts = ncont;
                        }

                        for ( long n=0; n<ncontacts; n++ ) {
                            long neighbor = ( long ) VECTOR ( neis ) [n];
                            long neighbor_status = ( long ) VAN ( graph, "S", ( igraph_integer_t ) neighbor );

                            if ( neighbor_status == IGNORANT ) {

                                if ( gsl_rng_uniform ( randomic ) < lambda ) {
                                    if ( gsl_rng_uniform ( randomic ) < eta ) {
                                        if ( VECTOR ( *state ) [neighbor] == IGNORANT ) {
                                            VECTOR ( *state ) [neighbor] = ( double ) SPREADER;
                                            population->Ignorant[t]--;
                                            population->Spreader[t]++;

                                            change_state = true;
                                        }
                                    } else {
                                        if ( VECTOR ( *state ) [neighbor] == IGNORANT ) {
                                            VECTOR ( *state ) [neighbor] = ( double ) STIFLER;
                                            population->Ignorant[t]--;
                                            population->Stifler[t]++;

                                            change_state = true;
                                        }
                                    }
                                }
                            } else if ( neighbor_status == SPREADER || neighbor_status == STIFLER ) {
                                if ( gsl_rng_uniform ( randomic ) < alpha ) {
                                    if ( VECTOR ( *state ) [i] != STIFLER ) {
                                        VECTOR ( *state ) [i] = ( double ) STIFLER;
                                        population->Spreader[t]--;
                                        population->Stifler[t]++;

                                        change_state = true;
                                    }
                                }
                            }


                        }

                        igraph_vector_destroy ( &neis );

                    }

                } else if ( status_i == STIFLER ) {
                    if ( rnd < gamma ) {
                        VECTOR ( *state ) [i] = IGNORANT;
                        population->Stifler[t]--;
                        population->Ignorant[t]++;
                        change_state = true;

                    } else {
                        igraph_vector_t neis;
                        igraph_vector_init ( &neis, 0 );
                        igraph_neighbors ( graph, &neis, ( igraph_integer_t ) i,  IGRAPH_OUT );

                        if ( igraph_vector_size ( &neis ) == 0 ) {
                            igraph_vector_destroy ( &neis );
                            continue;
                        }

                        igraph_vector_shuffle ( &neis );

                        long ncontacts;
                        if ( ncont >= ( long ) igraph_vector_size ( &neis ) ) {
                            ncontacts = igraph_vector_size ( &neis );
                        } else {
                            ncontacts = ncont;
                        }

                        for ( long n=0; n<ncontacts; n++ ) {

                            long neighbor = ( long ) VECTOR ( neis ) [n];
                            long neighbor_status = ( long ) VAN ( graph, "S", ( igraph_integer_t ) neighbor );

                            if ( neighbor_status == IGNORANT ) {
                                if ( gsl_rng_uniform ( randomic ) < beta ) {
                                    if ( VECTOR ( *state ) [i] != SPREADER ) {
                                        VECTOR ( *state ) [i] = ( double ) SPREADER;
                                        population->Spreader[t]++;
                                        population->Stifler[t]--;

                                        change_state = true;
                                    }
                                }

                            }
                        }

                        igraph_vector_destroy ( &neis );
                    }
                }
            }
        }



        if ( change_state == true ) {
            SETVANV ( graph,"S", state );
        }

    }
    
    igraph_vector_destroy(&vertex);

    return population;
}
