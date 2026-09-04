# Phylogenetic tree files

`01_primary_tree_grafen.tre` is Supplementary Data 3 and is the topology-based
Grafen primary tree used for the reported comparative analyses. Its branch
lengths represent relative node depth, not geological time.

Supplementary Data 4 contains the 14 alternatives used with that primary tree
in the 15-tree sensitivity analysis. Twelve alternatives are stored directly
in this directory. The two Grafen-power variants (`power = 0.5` and `power = 2`)
are generated deterministically from the primary topology by
`run_phylogenetic_comparative_analyses.R`. The release asset includes those
generated Newick files plus a manifest that maps every analysis `tree_id` to
its source.
