# automatically generated command file for RRTree

file of aligned sequences: exple.nexus
format: nexus
#specific to PHYLIP files
#name length: 10
lineages: 3

#lineage affectation of sequences
Human: 1
Mouse: 2
Sheep: 3
Pig: 3
Mink: 3
Guinea_pig: 2
Horse: 3
Chicken: 0

#lineage names
outgroup: Outgroup
lineage1: Human
lineage2: Rodents
lineage3: Ferungulates
sequence type: prot
topology: 1
tree file: exple.tree
threshold: 50
out text file: 0
out table file: 0
verbose: 0
