

#### slim set up
if (!require("devtools")) install.packages(devtools)
devtools::install_github("rdinnager/slimr")  #downloads the latest version

slim_setup() ### this will install slim if you don't have it already
slim_is_avail() ### really makes sure slim is there.

####The purpose of this script is to write a slim function where I can vary:
# - k number of populations
# - Fst (i.e. the length of time k populations are diverged)
# - number of markers (this can happen outside of slim with the vcf file and potentially just random thining)
# - the rate of migration between the populations
# - one thing to consider is whether we can vary the likelihood of the genotypes that are called 
### -- this would be great to test entropy against everything else, but would take some thinking about to jitter the calls realistically


### let's start by writing a k=2 hybridizing script that 1) estimates Fst and 2) prints a vcf
library(slimr)
# package adegenet as we are using genlight objects
#if (!require(adegenet)) install.packages("adegenet")
library(adegenet)

slim_script(
  slim_block (initialize(),
              {
                defineConstant("L", 1e6);
         initializeSLiMOptions(nucleotideBased=T);
         initializeAncestralNucleotides(randomNucleotides(L));
         initializeRecombinationRate(1e-8);
initializeMutationTypeNuc("m1",0.5,"f",0.0);
 initializeGenomicElementType("g1", m1, 1.0, mmJukesCantor(1e-8));
initializeGenomicElement(g1, 0, L-1);
 m1.color = "yellow";
 m1.colorSubstitution = "yellow";
 initializeMutationTypeNuc("m2", 0.5, "f", 0.01);
 m2.color = "red";
 m2.colorSubstitution = "red";
 initializeMutationTypeNuc("m3", 0.5, "f", 0.0);
 m3.color = "blue";
 m3.colorSubstitution = "blue";
 initializeMutationTypeNuc("m4", 0.5, "f", 0.01);
 m4.color = "green";
 m4.colorSubstitution = "green";
}),

  slim_block(1,
             {
               ## create 2 populations of 500 
               sim.addSubpop("p1", 500); ### add more populations here
               sim.addSubpop("p2", 500);
             }),
  
  slim_block(2999, late(), ### change here for initial divergence
             {
               #calculate Fst
              print(calcFST(p1.genomes, p2.genomes));
             }),
  
  
  slim_block(3500, late(),
             
             {
               p1.setMigrationRates(p2, 0.002); ### change here to change the rate of hybridization between populations
               p2.setMigrationRates(p1, 0.002);
             }),

  slim_block(3600, late(), ## change here for number of generations of hybridization
             
             {
               print(calcFST(p1.genomes, p2.genomes));
               g=c(p1.individuals, p2.individuals);
               g.genomes.outputVCF(filePath="working.vcf", simplifyNucleotides=T);
               sim.simulationFinished();
               
  })) %>% slim_run(capture_output = TRUE, show_output = TRUE) ->script_2_run

system("pwd")
system("cp /var/folders/s2/9217dxv55ds4_1szvkfz5qm40000gn/T//RtmpRS5K3z/file11e2825289576.txt ./slim_out.txt") ### this moves the random output file to the working directory

