import numpy as np
import argparse
from matplotlib import pyplot as plt
import pandas as pd

parser=argparse.ArgumentParser()
parser.add_argument('--p', type = float, help = 'initial observed frequency of allele A1, defaults to 1 heterozygous individual unless otherwise specificed')
parser.add_argument('--N', type = int, help = 'population size', default = 10000)
parser.add_argument('--t', type = int, help = 'number of generations', default = 100)
parser.add_argument('--reps', type = int, help = 'number of replicates', default = 10)
parser.add_argument('--s', type = float, help = 'selection coefficient', default = 0)
parser.add_argument('--h', type = float, help = 'heterozygote effect', default = 0)
parser.add_argument('--m', type = float, help = 'probability of migration (assumes immigration of individuals without allele, random emigration, equal population sizes)', default = None)
parser.add_argument('--u', type = float, help = 'probability of mutation', default = 2*(10**-8))
parser.add_argument('--out', type = str, help = 'output (appended with ".png, .txt")')
args=parser.parse_args()

def initialize (p, N): ## initialize GTs based on input allele frequency and population size
    gts = np.random.choice(a=[False, True], size = 2*N, p = [1-p, p])
    return gts 

def drift(p, gts): ## random drift
    next_gen = np.empty(len(gts))
    for i in range(len(gts)):
        next_gen[i] = gts[np.random.randint(len(gts))]
    return next_gen

def selection(gts, s, h): ## based on probability of repoduction (denoted by s,  h)
    p = np.sum(gts) / len(gts)    
    N = len(gts) / 2
    inds = np.array_split(gts, N)
    AA = 1
    Aa = 1 - h*s
    aa = 1-s
    next_gen = []
    while len(next_gen) < len(gts):
        genotype = inds[np.random.randint(len(inds))]
        if genotype.all() == True:
            if np.random.choice(a = [False, True], p = [1-AA, AA]) == True:
                next_gen.append(True)
        elif genotype.any() == True:
            if np.random.choice(a = [False, True], p = [1-Aa, Aa]) == True:
                next_gen.append(np.random.choice(a = [False, True], p = [.5, .5]))
        elif genotype.all() == False:
            if np.random.choice(a = [False, True], p = [1-aa, aa]) == True:
                next_gen.append(False)
    return np.array(next_gen)

def migration(gts, m):
    emm = round(len(gts) * m)    
    if (emm % 2) != 0:
        emm = emm - 1
    next_gen = list(np.delete(gts, np.random.randint(len(gts), size = emm)))
    imm = round(m*len(gts))
    for hap in range(imm):
        next_gen.append(False)
    return np.array(next_gen)

def mutation(gts, u): ## probabiliy of back/forward mutation
    mutants = np.empty(len(gts))
    for i in range(len(gts)):
        mute = np.random.choice(a = [False, True], p = [(1-u), u])
        if mute == False:
            mutants[i] = gts[i]
        elif mute == True and gts[i] == False:
            mutants[i] = True
        elif mute == True and gts[i] == True:
            mutants[i] = False
    return mutants

df = pd.DataFrame(columns = range(args.reps), index = range(args.t)) ## initialize dataframe

for i in range(args.reps):
    ## initialize p
    if args.p == None:
        p = 1/(2*args.N)
    else:
        p = args.p
    ## initialize gts
    gts = initialize(p, args.N)
    df.loc[0, i] = np.sum(gts) / (2 * args.N)
    count = 1
    for g in range(args.t): ## for generation
        gts = gts
        p = df.loc[g, i] ## determine allele frequency
        gts = drift(p, gts) ## drift
        gts = selection(gts, args.s, args.h) ## selection
        if args.m != None:
            gts = migration(gts, args.m) ## migration
            gts = mutation(gts, args.u) ## mutation
        else:
            gts = mutation(gts, args.u)
        df.loc[count, i] = np.sum(gts) / len(gts) ## return allele frequency
        count += 1
df.to_csv(args.out + '.txt', sep = '\t')

fig, ax = plt.subplots() ## plot allele frequency over generations
colors = ["#d1bce3","#c49bbb","#a1867f","#585481","#19297c","#35524a","#a2e8dd","#32de8a","#292f36","#4ecdc4"]
for i in range(args.reps):
    ax.plot(df[i], c = colors[np.random.randint(len(colors))])
plt.savefig(args.out + '.png')

        
        
        
