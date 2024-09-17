

#### slim set up
if (!require("devtools")) install.packages(devtools)
devtools::install_github("rdinnager/slimr")  #downloads the latest version
library(slimr)

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
# package adegenet as we are using genlight objects
#if (!require(adegenet)) install.packages("adegenet")
library(adegenet)

migration_choices<-c(0.002, 0.02, 0.2) ### these are the three levels that I used in my Mol Ecol paper
number_of_pops<-c(2,3) ### need to extend this and all of the loops out to k=5
pop_size<-c(1000, 10000, 100000)


for(i in 1:length(migration_choices)){
  
  for(j in 1:length(number_of_pops)){
    
    for(k in 1:length(pop_size)){
    
slim_script(
  slim_block (initialize(),
              {
                defineConstant("L", 1e9);
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
               Npop<-r_inline(number_of_pops[j])
               popN<-r_inline(pop_size[k])
               if(Npop==2){
               sim.addSubpop("p1", popN); ### add more populations here - This is probably going to be an ifelse statement from the loop
               sim.addSubpop("p2", popN);}
              if(Npop==3){  
                print("3!")
                sim.addSubpop("p3", popN); ### add more populations here - This is probably going to be an ifelse statement from the loop
               sim.addSubpop("p4", popN);
               sim.addSubpop("p5", popN)}### this works, just need to figure out what to do with it!
              }),
  
  slim_block(2999, late(), ### change here for initial divergence
             {
               #calculate Fst
               Npop<-r_inline(number_of_pops[j])
               if(Npop==2)
              print(calcFST(p1.genomes, p2.genomes));
               if(Npop==3){
                 print("yay 3!")
                print(calcFST(p3.genomes, p4.genomes));
               print(calcFST(p3.genomes, p5.genomes));
               print(calcFST(p4.genomes, p5.genomes))}
             }),
  
  
  slim_block(3500, late(),
             
             {
               m = r_inline(migration_choices[i])
               Npop<-r_inline(number_of_pops[j])
               
               if(Npop==2){
               p1.setMigrationRates(p2, m); ### caveat here is that the migration between the populations is equal. It doesn't have to be, but it is
               p2.setMigrationRates(p1, m);}
               
               if(Npop==3){
                p3.setMigrationRates(p4, m); 
               p3.setMigrationRates(p5, m); 
               p4.setMigrationRates(p3, m); 
               p4.setMigrationRates(p5, m); 
               p5.setMigrationRates(p3, m); 
               p5.setMigrationRates(p3, m)}
             }),

  slim_block(3505, late(), ## 10 generations of hybridization - can change if needed, but likely won't vary this.
             
             {
               Npop<-r_inline(number_of_pops[j])
               if(Npop==2){
               print(calcFST(p1.genomes, p2.genomes));
               g=c(p1.individuals, p2.individuals);}
               if(Npop==3){
                 print(calcFST(p3.genomes, p4.genomes));
                 print(calcFST(p3.genomes, p5.genomes));
                 print(calcFST(p4.genomes, p5.genomes));
                 g=c(p3.individuals, p4.individuals, p5.individuals);}
               file_path = r_inline(paste0("test", "m", migration_choices[i], "pop", number_of_pops[j],"N",pop_size[k], ".vcf")) ### also need to change here
               g.genomes.outputVCF(filePath=file_path, simplifyNucleotides=T);
               sim.simulationFinished();
               
  })) %>% slim_run(capture_output = TRUE, show_output = TRUE) ->script_2_run

outputcommand<-paste("cp", script_2_run$output_file, paste0("./", "m", migration_choices[i],"pop", number_of_pops[j],"N", pop_size[k], ".txt"))
system(outputcommand) ### this moves the random output file to the working directory

  }
}
}

#slim_function("test_name", 1)

              