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

#include <string.h>
#include <omp.h>
#include <unistd.h>

#include <time.h>

#include <igraph/igraph.h>

#include "Markov.h"
#include "MonteCarlo.h"
#include <igraph/igraph.h>
#include <igraph/igraph.h>

const int NEXP = 1000;

int PrintInfo ( int type )
{
    if ( type == 1 ) {
        printf ( "General Markov Chain model\n\n" );
        printf ( "Please cite:\n" );
    }

    if ( type >= 1 ) {
        printf ( "A general Markov chain approach for disease and rumour spreading in complex networks\n" );
        printf ( "Guilherme Ferraz de Arruda, Francisco Aparecido Rodrigues, Pablo Martín Rodríguez, Emanuele Cozzo, Yamir Moreno\n" );
        printf ( "Journal of Complex Networks, cnx024, https://doi.org/10.1093/comnet/cnx024\n\n" );
        printf ( "DTMC, Copyright (C) 2017  Guilherme Ferraz de Arruda\n" );
        printf ( "DTMC comes with ABSOLUTELY NO WARRANTY;\n" );
        printf ( "This is free software, and you are welcome to redistribute it under certain conditions.\n\n" );
    }

    if ( type == 1 ) {
        printf ( "Use:\n" );
        printf ( "   dtmc network.edgelist \\lambda \\alpha \\delta_1 \\delta_2 \\gamma \\eta \\beta Steps [cp/rp] [markov/montecarlo/both] \n" );
        printf ( "   dtmc network.edgelist \\lambda \\alpha \\delta_1 \\delta_2 \\gamma \\eta \\beta Steps [P_0] [cp/rp] [markov/montecarlo/both] \n" );
        printf ( "     markov: runs the numerical simularions solving the dynamical equations iterativelly\n" );
        printf ( "     montecarlo: runs the Monte Carlo simularions\n" );
        printf ( "     both: runs the markov and montecarlo options\n" );
        printf ( "     cp: contact processes\n" );
        printf ( "     rp: reactive processes\n" );
        printf ( "     \\lambda, \\alpha, \\delta_1, \\delta_2, \\gamma, \\eta, \\beta: are real parameters\n" );
        printf ( "     Steps: number of discrete time steps (integer)\n" );
        printf ( "     P_0: File with the intial nodal probability (optional parameter)\n" );
        printf ( "     network.edgelist: edgelist of the network under analysis\n" );

        printf ( "\n" );

        printf ( "Acknowledgement:\n" );
        printf ( "  This was coded using the igraph library:\n" );
        printf ( "  Gábor Csárdi, Tamás Nepusz: The igraph software package for complex network research. InterJournal\n" );
        printf ( "  Complex Systems, 1695, 2006." );
        printf ( "  For more, please see http://igraph.org/c/\n" );

        printf ( "\n" );
    }

    return 1;
}

int main ( int argc, char **argv )
{
    if ( strcmp ( argv[argc-1], "help" ) == 0 ) {
        PrintInfo ( 1 );

        return 1;
    }

    time_t start, end;
    time ( &start );

    char fname[256];

    igraph_i_set_attribute_table ( &igraph_cattribute_table );

    const gsl_rng_type * Tr;
    gsl_rng * randomic;
    gsl_rng_env_setup();
    Tr = gsl_rng_default;
    randomic = gsl_rng_alloc ( Tr );
    gsl_rng_set ( randomic, time ( NULL ) * getpid() );

    PrintInfo ( 1 );

    igraph_t graph;
    FILE *file = fopen ( argv[1], "r" );
    if ( file == NULL ) {
        printf ( "Could not find file %s\n", argv[1] );
        return 0;
    }
    igraph_read_graph_edgelist ( &graph, file, 0, IGRAPH_UNDIRECTED );
    fclose ( file );

    long N = ( long ) igraph_vcount ( &graph );

    double spreader0 = 0.05, stifler0 = 0;
    double ignorant0 = 1 - spreader0;
    double lambda = atof ( argv[2] ), alpha = atof ( argv[3] ), delta1 = atof ( argv[4] ), delta2 = atof ( argv[5] ), gamma = atof ( argv[6] );
    double eta = atof ( argv[7] ), beta = atof ( argv[8] );
    long Steps = atoi ( argv[9] );

    printf ( "Graph characterization:\n" );
    printf ( "Open: %s\n", argv[1] );
    printf ( "N = %ld\n", N );

    printf ( "\nRumor parameters:\n" );
    printf ( "lambda = %f\n", lambda );
    printf ( "alpha = %f\n", alpha );
    printf ( "delta_1 = %f\n", delta1 );
    printf ( "delta_2 = %f\n", delta2 );
    printf ( "delta_3 = %f\n", gamma );
    printf ( "beta = %f\n", beta );
    printf ( "eta = %f\n", eta );
    printf ( "T_max = %ld\n", Steps );

    // Check for inputs consitency
    if ( lambda < 0 || alpha < 0 || delta1 < 0 || delta2 < 0 || gamma < 0 || eta < 0 || beta < 0 ) {
        printf ( "Error: Negative probabilites\n" );
        return 0;
    }

    if ( abs ( lambda ) > 1 || abs ( alpha ) > 1 || abs ( delta1 ) > 1 || abs ( delta2 ) > 1 || abs ( gamma ) > 1 || abs ( eta ) > 1 || abs ( beta ) > 1 ) {
        printf ( "Error: Probabilites must be <= 1\n" );
        return 0;
    }

    if ( ( delta1 + delta2 ) > 1 ) {
        printf ( "Error: (\\delta_1 + \\delta_2) must be <= 1\n" );
        printf ( "For more see:\n" );
        PrintInfo(2);
        return 0;
    }

    if ( strcmp ( argv[argc-1], "markov" ) == 0 || strcmp ( argv[argc-1], "both" ) == 0 ) {
        printf ( "\nDTMC -- Discrete Time Markov Chain approach:\n" );

        if ( strcmp ( argv[argc-2], "cp" ) == 0 ) {
            printf ( "Contact Process (CP)\n" );
        } else {
            printf ( "Reactive Process (RP)\n" );
        }

        fflush ( stdout );

        igraph_vector_t k;
        igraph_vector_init ( &k, 0 );
        igraph_degree ( &graph, &k, igraph_vss_all(), IGRAPH_OUT, 0 );

        igraph_vector_t *adjl = ( igraph_vector_t* ) calloc ( N, sizeof ( igraph_vector_t ) );

        for ( long i=0; i<N; i++ ) {
            igraph_vector_init ( &adjl[i], 0 );
            igraph_neighbors ( &graph, &adjl[i], i, IGRAPH_OUT );
        }


        FILE *fp_in_nodes;
        if ( argc == 13 ) {
            fp_in_nodes = fopen ( argv[argc-3], "r" );
            printf ( "File with the intial nodal probability: %s\n", argv[argc-3] );
        } else if ( argc == 14 ) {
            stifler0 = atof ( argv[argc-4] );
            spreader0 = atof ( argv[argc-3] );
            ignorant0 = 1 - spreader0 - stifler0;

            printf ( "Initial conditions\n" );
            printf ( " Spreader: %.12le\n Ignorant: %.12le\n Stifler:  %.12le\n", spreader0, ignorant0, stifler0 );
        }


        double **P_ignorant = ( double ** ) calloc ( N, sizeof ( double* ) );
        double **P_spreader = ( double ** ) calloc ( N, sizeof ( double* ) );
        double **P_stifler = ( double ** ) calloc ( N, sizeof ( double* ) );

        for ( long i=0; i<N; i++ ) {
            P_ignorant[i] = ( double * ) calloc ( Steps, sizeof ( double ) );
            P_spreader[i] = ( double * ) calloc ( Steps, sizeof ( double ) );
            P_stifler[i] = ( double * ) calloc ( Steps, sizeof ( double ) );

            if ( argc == 13 ) {
                double tmp_ignorant = 0, tmp_stifler = 0, tmp_spreader = 0;

                fscanf ( fp_in_nodes, "%le\t%le\t%le\n", &tmp_ignorant, &tmp_stifler, &tmp_spreader );

                P_ignorant[i][0] = tmp_ignorant;
                P_spreader[i][0] = tmp_spreader;
                P_stifler[i][0] = tmp_stifler;

            } else {
                P_ignorant[i][0] = ignorant0;
                P_spreader[i][0] = spreader0;
                P_stifler[i][0] = stifler0;
            }
        }

        if ( argc == 13 ) {
            fclose ( fp_in_nodes );
        }


        for ( long t=0; t<Steps-1; t++ ) {
            double *a = ( double* ) calloc ( N, sizeof ( double ) );
            double *b = ( double* ) calloc ( N, sizeof ( double ) );
            double *c = ( double* ) calloc ( N, sizeof ( double ) );

            for ( long i=0; i<N; i++ ) {
                a[i] = 1;

                if ( strcmp ( argv[argc-2], "cp" ) == 0 ) {
                    b[i] = 0;		// CP
                    c[i] = 0;		// CP
                } else {
                    b[i] = 1;		// RP
                    c[i] = 1;		// RP
                }
            }

            if ( strcmp ( argv[argc-2], "cp" ) == 0 ) {
                ContactProbabilitiesCP ( N, t, adjl, lambda, alpha, delta1, delta2, gamma, eta, beta, P_ignorant, P_spreader, P_stifler, a, b, c, &k );
            } else {
                ContactProbabilitiesRP ( N, t, adjl, lambda, alpha, delta1, delta2, gamma, eta, beta, P_ignorant, P_spreader, P_stifler, a, b, c );
            }

            TotalProbability ( N, t, adjl, lambda, alpha, delta1, delta2, gamma, eta, beta, P_ignorant, P_spreader, P_stifler, a, b, c );

            free ( a );
            free ( b );
            free ( c );
        }

        fflush ( stdout );

        sprintf ( fname, "%s_Time_RumorMarkov_lambda_%.6f_alpha_%.6f_delta1_%.6f_delta2_%.6f_gamma_%.6f_eta_%.6f_beta_%.6f_tmax_%ld.dat", argv[1], lambda, alpha, delta1, delta2, gamma, eta, beta, Steps );
        FILE *fp_out_t = fopen ( fname, "w" );

        sprintf ( fname, "%s_Final_Nodal_RumorMarkov_lambda_%.6f_alpha_%.6f_delta1_%.6f_delta2_%.6f_gamma_%.6f_eta_%.6f_beta_%.6f_tmax_%ld.dat", argv[1], lambda, alpha, delta1, delta2, gamma, eta, beta, Steps );
        FILE *fp_out_tmax = fopen ( fname, "w" );

        for ( long t=0; t<Steps; t++ ) {
            double avg_ig = 0, avg_sp = 0, avg_st = 0;
            for ( long i=0; i<N; i++ ) {
                avg_ig += P_ignorant[i][t];
                avg_sp += P_spreader[i][t];
                avg_st += P_stifler[i][t];

                if ( t == Steps-1 ) {
                    fprintf ( fp_out_tmax, "%le\t%le\t%le\n", P_ignorant[i][t], P_stifler[i][t], P_spreader[i][t] );
                }

            }

            fprintf ( fp_out_t, "%e\t%e\t%e\n", avg_ig, avg_st, avg_sp );
        }

        fclose ( fp_out_t );
        fclose ( fp_out_tmax );

        for ( long i=0; i<N; i++ ) {
            free ( P_ignorant[i] );
            free ( P_spreader[i] );
            free ( P_stifler[i] );
        }

        free ( P_ignorant );
        free ( P_spreader );
        free ( P_stifler );

        for ( long i=0; i<N; i++ ) {
            igraph_vector_destroy ( &adjl[i] );
        }

        free ( adjl );


    }

    if ( strcmp ( argv[argc-1], "montecarlo" ) == 0 || strcmp ( argv[argc-1], "both" ) == 0 ) {
        igraph_vector_t initial[NEXP];
        rates *R[NEXP];

        long* MicroIgnorant = ( long * ) calloc ( N, sizeof ( long ) );
        long* MicroSpreader = ( long * ) calloc ( N, sizeof ( long ) );
        long* MicroStifler = ( long * ) calloc ( N, sizeof ( long ) );

        for ( long i=0; i<N; i++ ) {
            MicroIgnorant[i] = 0;
            MicroSpreader[i] = 0;
            MicroStifler[i] = 0;
        }

        igraph_t *g = ( igraph_t* ) calloc ( NEXP, sizeof ( igraph_t ) );

        // Probabilities over time
        sprintf ( fname, "%s_RumorMC_lambda_%.6f_alpha_%.6f_delta1_%.6f_delta2_%.6f_gamma_%.6f_eta_%.6f_beta_%.6f_tmax_%ld.dat", argv[1], lambda, alpha, delta1, delta2, gamma, eta, beta, Steps );
        FILE *fp = fopen ( fname, "w" );

        double* Ignorant = ( double * ) calloc ( Steps, sizeof ( double ) );
        double* Ignorant_2 = ( double * ) calloc ( Steps, sizeof ( double ) );
        double* Spreader = ( double * ) calloc ( Steps, sizeof ( double ) );
        double* Spreader_2 = ( double * ) calloc ( Steps, sizeof ( double ) );
        double* Stifler = ( double * ) calloc ( Steps, sizeof ( double ) );
        double* Stifler_2 = ( double * ) calloc ( Steps, sizeof ( double ) );


        #pragma omp parallel for schedule(dynamic,1)
        for ( long n0=0; n0<NEXP; n0++ ) {
            igraph_copy ( &g[n0], &graph );

            igraph_vector_init ( &initial[n0], N );

            for ( long i=0; i<N; i++ ) {
                if ( i< floor ( ignorant0*N ) ) {
                    VECTOR ( initial[n0] ) [i] = IGNORANT;
                } else if ( i< ( ignorant0+spreader0 ) *N ) {
                    VECTOR ( initial[n0] ) [i] = SPREADER;
                } else {
                    VECTOR ( initial[n0] ) [i] = STIFLER;
                }
            }

            igraph_vector_shuffle ( &initial[n0] );

            if ( strcmp ( argv[argc-2], "cp" ) == 0 ) {
                R[n0] = MonteCarlo ( &g[n0], Steps, lambda, alpha, delta1, delta2, gamma, eta, beta, 1, &initial[n0], randomic );
            } else {
                R[n0] = MonteCarlo ( &g[n0], Steps, lambda, alpha, delta1, delta2, gamma, eta, beta, N, &initial[n0], randomic );
            }

            printf ( "EXP: %ld/%d\n", n0+1, NEXP );
            printf ( "T = 0;     Spreader: %ld, Ignorant: %ld, Stifler: %ld;\n", R[n0]->Spreader[0], R[n0]->Ignorant[0], R[n0]->Stifler[0] );
            printf ( "T = T_max; Spreader: %ld, Ignorant: %ld, Stifler: %ld;\n", R[n0]->Spreader[Steps-1], R[n0]->Ignorant[Steps-1], R[n0]->Stifler[Steps-1] );
            fflush ( stdout );

            // Micro states
            for ( long i=0; i<N; i++ ) {
                if ( ( long ) VAN ( &g[n0], "S", ( igraph_integer_t ) i ) == IGNORANT ) {
                    #pragma omp atomic
                    MicroIgnorant[i]++;
                } else if ( ( long ) VAN ( &g[n0], "S", ( igraph_integer_t ) i ) == SPREADER ) {
                    #pragma omp atomic
                    MicroSpreader[i]++;
                } else {
                    #pragma omp atomic
                    MicroStifler[i]++;
                }
            }

            // Time
            for ( long i=0; i<Steps; i++ ) {
                #pragma omp atomic
                Ignorant[i] += R[n0]->Ignorant[i];
                Ignorant_2[i] += R[n0]->Ignorant[i]*R[n0]->Ignorant[i];

                Spreader[i] += R[n0]->Spreader[i];
                Spreader_2[i] += R[n0]->Spreader[i]*R[n0]->Spreader[i];

                Stifler[i] += R[n0]->Stifler[i];
                Stifler_2[i] += R[n0]->Stifler[i]*R[n0]->Stifler[i];
            }

            FreeRate ( R[n0] );
            igraph_vector_destroy ( &initial[n0] );

            igraph_destroy ( &g[n0] );
        }


        // Micro probabilites at steady-state
        sprintf ( fname, "%s_MicroFinal_RumorMC_lambda_%.6f_alpha_%.6f_delta1_%.6f_delta2_%.6f_gamma_%.6f_eta_%.6f_beta_%.6f_tmax_%ld.dat", argv[1], lambda, alpha, delta1, delta2, gamma, eta, beta, Steps );
        FILE *fpMicroFinal = fopen ( fname, "w" );


        for ( long i=0; i<Steps; i++ ) {
            fprintf ( fp, "%e\t%e\t%e\t%e\t%e\t%e\n",
                      ( double ) Ignorant[i]/ ( double ) ( NEXP ), Ignorant_2[i]/ ( double ) ( NEXP ),
                      ( double ) Stifler[i]/ ( double ) ( NEXP ), Stifler_2[i]/ ( double ) ( NEXP ),
                      ( double ) Spreader[i]/ ( double ) ( NEXP ), Spreader_2[i]/ ( double ) ( NEXP ) );
            fflush ( fp );
        }


        for ( long i=0; i<N; i++ ) {
            fprintf ( fpMicroFinal, "%e\t%e\t%e\n", ( double ) MicroIgnorant[i]/ ( double ) ( NEXP ), ( double ) MicroStifler[i]/ ( double ) ( NEXP ), ( double ) MicroSpreader[i]/ ( double ) ( NEXP ) );

            fflush ( fpMicroFinal );
        }


        // Free dynamic vectors
        free ( MicroSpreader );
        free ( MicroIgnorant );
        free ( MicroStifler );

        free ( Ignorant );
        free ( Ignorant_2 );
        free ( Stifler );
        free ( Stifler_2 );
        free ( Spreader );
        free ( Spreader_2 );

    }


    igraph_destroy ( &graph );

    time ( &end );
    printf ( "Finished in: %lf seconds\n\n", difftime ( end, start ) );

    return 1;
}

