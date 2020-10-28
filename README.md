# quantgen
code based on Quantitative Genetics (Caballero A. 2020)

based on chance, should be used for entertainment only and should not be used for investment purposes


### stochastic 2-allele mutation/selection/drift/migration model
`model.py`: stochastic 2-allele frequency dynamics

#### input: 
initial population parameters

   `--p `  initial observed frequency of allele A1 (defaults to 1 heterozygous individual unless otherwise specified)
  
   `--N `  population size (default = 10000)
  
   `--t `  number of generations (default = 100)
  
   `--reps ` number of replicate populations (default = 10)
  
   `--s `  selection coefficient (default = 0)
  
   `--h `  heterozygous effect (default = 0)
  
   `--m `  probability of migration (assumes immigration of individuals homozygous for allele A2, random emigration, default = 0)
  
   `--u `  probability of mutation (default = 2 * 10e-8)
  
   `--out `  output (appended with .png, .txt)
   
#### output:
plot of allele frequency over generations, .tsv tracking allele frequency across replicates over generations

   
#### example: drift
   python model.py --out test --N 100 --p 0.5 --out test
