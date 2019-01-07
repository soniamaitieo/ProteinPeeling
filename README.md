
# Protein Peeling

The objective of this tool is to divide a protein into small compact units that compose protein three-dimensional structures. These fragments, called protein units (PU), are a new level of description to well understand and analyze the organization of protein structures defined by Jean-Christophe Gelly's team.  

The method only works from the contact probability matrix of inter-Calpha distances. Distances are transformed into probabilities with a logistic function. The original programm split the protein by hierarchical clustering.  


> [Reference 1](https://www.ncbi.nlm.nih.gov/pubmed/16301202) :
Jean-Christophe Gelly, Alexandre G. de Brevern, and Serge Hazout. ’Protein Peeling’ :
an approach for splitting a 3d protein structure into compact fragments.
Bioinforma-tics (Oxford, England), 22(2) :129–133, January 2006

The original version of Protein Peeling developped by Jean-Christophe Gelly's team is here :  [Protein Peeling 3D](http://www.dsimb.inserm.fr/dsimb_tools/peeling3/)



The aim of this university project is to propose a new method of Protein Peeling without hierarchical segmentation but by searching iteratively for the regions of the most compact and independent proteins.  

The compactness (κ) and the separation(σ) are defined in this paper :

> [Reference 2](http://advances.sciencemag.org/content/3/1/e1600552.full) :
Guillaume Postic, Yassine Ghouzam, Romain Chebrek, and Jean-Christophe Gelly.
An ambiguity principle for assigning protein structural domains.
Science Advances,3(1), 2017.



## Requirements

* PyMol

```
 sudo apt-get install pymol
 ```
 * Python library

```
pip install numpy
pip install pandas
pip install matplotlib
pip install biopython

```
## Run

Please launch the tool at base.

Run the program
```
 python3 Protein_Peeling.py PDBfile CHAIN MIN MAX
 ```    

With :

- PDBfile : the path to your PDBfile, we recommand you to put your pdb file in `data/` folder

- CHAIN : the Chain of your protein (ex: A)

- MIN : min size of PU (ex:10)

- MAX: max size of PU (ex:50)


## Example usage


```
 python3 Protein_Peeling.py /datas/1atn.pdb A 10 50
 ```    

 The tool will create a folder in `results/` with the name of pdb (ex: `1atn` ).  
 Inside `results/1atn`,  different folders "PUx" (x: number of the splitting) with all the possibilities to split the protein are created ( `results/1atn/PU1/`).  
 Inside each possibilities, pdbfiles and pymol reprensentation are generated.
 (the representation is made with `src/visu_all_pdbs.py` and Pymol )   
 It will print in the terminal the fold wich countain the best splitting and plot contact map.
The file `ranking_spitting.txt ` contains the rank and the information of each splitting.

 <center>
 **Protein 1atn with the best spliiting (foldname:PU14) colored by it's PUs**  

![Illustration](/results/1atn/PU1/PU1.png)  

**PU example (residues:202-251) from splitting 12 (foldname:PU14)**  

![Illustration](/results/1atn/PU1/PU1.289-334.pdb.png)  

**Contact map divided according to the all spliiting combinations**  

![Illustration](/results/1atn/contactmap_all.gif)  

</center>


## Documentation

The Documentation is made with Pydoc.
```
pydoc -w ProteinPeeling
```
