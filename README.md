----------------------------
# 2023-npc-grafting-distance

### Synopsis: 
Analyses of grafting distance between FG Nup anchor domains based on Kim et all. 2018 paper 

### Manuscript:
Dynamic molecular mechanism of the nuclear pore complex permeability barrier
Kozai et al., 2023 https://www.biorxiv.org/content/10.1101/2023.03.31.535055v1

### Author: 
Barak Raveh https://ravehlab.org

### Last updated: 
April 13, 2023

----------------------------
Reproducing the analysis:
-------------------------

1. Install Python
2. Install IMP (https://integrativemodeling.org/), either from source code or from packages such as conda (https://anaconda.org/conda-forge/imp), or an IMP installation over google colab (see tutorials on IMP website)
3. Download this folder
4. Using a python that has IMP properly installed, just run from this folder:

$ python analyze_anchor_distance.py 47-35_1spoke.rmf3 1 distance.csv

--------------
General Notes:
--------------

1. We used the average distance to the nearest anchoring site within a certain category within a single spoke (e.g. all FG nups, or only cytoplasmic subunits). For example, if we have three anchoring sites in a single spoke, the nearest neighbor (NN) of the first being 2 nm away, the NN of the second being 4 nm away, and the NN of the third being 5 nm away, then s in the calculation is 3.66 nm (just as a fake example, obviously).
2. We omitted Nup42 - it is too dynamic to be considered here
3. Nup2 is not in the map (we assume it colocalizes with Nup60, based on PDB of complex)
4. We used the following anchoring sites (the residue ranges correspond to a single sphere in the Kim et al. 2018 model for each Nup subunit). We also indicate what copies (subunits) are assigned to what parts of the NPC
Nsp1.601-636 (copies 1/2 - cytoplasmic; copies 3/4 - central)
Nup1.301-350 (nuclear)
Nup49.201-269 (both central)
Nup57.201-286 (both central)
Nup60.351-398 (both nuclear)
Nup100.551-575_ (copy 1 - central; copy 2 - cytoplasmic)
Nup116.751-775 (both cytoplasmic)
Nup145.201-225 (copy 1 - nuclear; copy 2 - central)
Nup159.1082-1116 (both cytoplasmic)
