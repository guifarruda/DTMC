# DTMC

DTMC,  (C) 2017  Guilherme Ferraz de Arruda
DTMC comes with ABSOLUTELY NO WARRANTY;
This iCopyrights free software, and you are welcome to redistribute it under certain conditions.

## General Markov Chain model

## Please cite:
A general Markov chain approach for disease and rumour spreading in complex networks
Guilherme Ferraz de Arruda, Francisco Aparecido Rodrigues, Pablo Martín Rodríguez, Emanuele Cozzo, Yamir Moreno
Journal of Complex Networks, cnx024, https://doi.org/10.1093/comnet/cnx024.
Alternativelly, see https://arxiv.org/abs/1609.00682.

## Depends on:
  1. igraph (http://igraph.org/c/)
  2. gls (https://www.gnu.org/software/gsl/)
  3. OpenMP (http://www.openmp.org/)
  4. Cmake (https://cmake.org/)

## Use:
<p>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;dtmc network.edgelist \lambda \alpha \delta_1 \delta_2 \gamma \eta \beta Steps [cp/rp] [markov/montecarlo/both]  <br/>
<p>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;dtmc network.edgelist \lambda \alpha \delta_1 \delta_2 \gamma \eta \beta Steps [P_0] [cp/rp] [markov/montecarlo/both]  <br/>
<p>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;markov: runs the numerical simularions solving the dynamical equations iterativelly <br/>
<p>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;montecarlo: runs the Monte Carlo simularions <br/>
<p>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;both: runs the markov and montecarlo options <br/>
<p>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;cp: contact processes <br/>
<p>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;rp: reactive processes <br/>
<p>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;\lambda, \alpha, \delta_1, \delta_2, \gamma, \eta, \beta: are real parameters <br/>
<p>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Steps: number of discrete time steps (integer) <br/>
<p>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;P_0: File with the intial nodal probability (optional parameter) <br/>
<p>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;network.edgelist: edgelist of the network under analysis <br/>

## Examples:
  1. Run the compile_run_examples.sh to run the 3 examples on real networks. 
  2. The networks used in those examples were also provided here and their appropriate citation are given on the paper.
  3. For more, see "A general Markov chain approach for disease and rumour spreading in complex networks", Journal of Complex Networks, cnx024, https://doi.org/10.1093/comnet/cnx024. Alternativelly, see https://arxiv.org/abs/1609.00682.

## Acknowledgement:
  This was coded using the igraph library:
  Gábor Csárdi, Tamás Nepusz: The igraph software package for complex network research. InterJournal
  Complex Systems, 1695, 2006.  For more, please see http://igraph.org/c/

