## Physicochemical and Medicinal Chemistry data extraction

This post is a step-by-step implementation of my approach to calculating the Physicochemical properties and Medicinal Chemistry and it is meant for sharing, studying, and critiquing by fellow researchers who are new and interested in this topic.

### **Physicochemical Properties:**

*   [x] **Molecular Weight**
*   [ ] **Volume**
*   [ ] **Density**
*   [x] **Heavy Atoms**
*   [x] **Aromatic Heavy Atoms**
*   [x] **Fraction Csp3**
*   [x] **Rotatable Bonds**
*   [x] **Hetro Atoms**
*   [x] **H-Bond Acceptors**
*   [x] **H-Bond Donors**
*   [x] **Ring Count**
*   [x] **Aromatic Ring Count**
*   [x] **Stereo Centers**
*   [x] **Molar Refractivity**
*   [x] **tPSA**
*   [x] **LogP**
*   [ ] **pKa**
*   [x] **LogS**
*   [ ] **LogD7.4**

---

### **Medicinal Chemistry:**

*   [x] **QED**
*   [x] **SAscore**
*   [x] **NPScore**
*   [x] **Fsp3**
*   [x] [**Lipinski Violations**](https://en.wikipedia.org/wiki/Lipinski%27s_rule_of_five)
<!-- *   [ ] [**Pfizer Violations**](https://github.com/santuchal/adme_predection/blob/master/ref/hughes2008.pdf)
*   [ ] [**GSK Violations**](https://github.com/santuchal/adme_predection/blob/master/ref/gleeson2008.pdf)
*   [ ] [**Golden Triangle**](https://github.com/santuchal/adme_predection/blob/master/ref/johnson2009.pdf) -->
*   [x] [**Ghose Violations**](https://github.com/santuchal/adme_predection/blob/master/ref/ghose1999.pdf)
*   [x] [**Veber Violations**](https://github.com/santuchal/adme_predection/blob/master/ref/veber2002.pdf)
*   [x] [**Muegge Violations**](https://github.com/santuchal/adme_predection/blob/master/ref/muegge2001.pdf)
*   [x] [**Egan Violations**](https://github.com/santuchal/adme_predection/blob/master/ref/egan2000.pdf)
*   [x] [**PAINS Alerts**](https://en.wikipedia.org/wiki/Pan-assay_interference_compounds)
*   [x] **BRENK Alerts**
*   [ ] **ALARM NMR Rule**
*   [ ] **Chelator Rule**

---

### **Installation**

You need to install following library:

* RdKit
* Sqlite3
* Numpy
* Pandas
* Tqdm
* docopt

### **Run**

final.py --classification <MOLECULE> --in <INFILE_NAME> --out <OUTFILE_NAME> --save_as <TYPE>

* --classification type small_molecule | peptide | big_molecule
* --in INFILE_NAME  input SMILES(default),Amino Acid sequence or Nuclic Acid Sequence in single column CSV file
* --save_as csv(default) | excel | sqlite3 | 
* --out OUTFILE_NAME  output csv 

* Example : 

python commandline_test.py --classification peptide --in ../data/test_peptide.csv --out output --save_as sqlite3

python commandline_test.py --classification peptide --in ../data/test_peptide.csv --out output --save_as excel

python commandline_test.py --classification peptide --in ../data/test_peptide.csv --out output --save_as csv

python commandline_test.py --classification small_molecule --in ../data/temp.smi --out output --save_as sqlite3

python commandline_test.py --classification small_molecule --in ../data/temp.smi --out output --save_as excel

python commandline_test.py --classification small_molecule --in ../data/temp.smi --out output --save_as csv


