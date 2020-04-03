# sim-micro-genome-evo
Forward-simulation of microbial genome evolution

Microbial genomes (e.g., those of bacterial and viral pathogens) differ from eukaryotic genomes in having mostly protein-coding sequences. As a result, simulation software assuming mostly non-coding, neutrally-evolving genome sequences are not suitable for understanding microbial genome evolution.

Specifically, microbial genome evolution has the following characteristics that warantee a specialized simulation tool:

1. Haploid, mostly protein coding sequences. No introns.

2. Strong selection & weak drift due to large population sizes. Pesence of strong purifying as well as adaptive natural selection preclude the use of neutral, efficient coalescence-based simulators like ms.

3. Recombination rates at the order of mutation rates (not as frequent as in sexually-reproducing eukaryotes).

4. Recombination occurs by gene conversion (not exchange).

This repository consists of Perl modules (made with Moose) for forward-simulation of microbial genome evolution at population level. It also include calling scripts and post-processing scripts. It has previously been used to produce results of this publication:  http://www.genetics.org/content/189/3/951.long
