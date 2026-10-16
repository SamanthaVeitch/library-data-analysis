created 
- environment.yml
- README.md
- .gitignore

<div style="margin-top: 50px;"></div>

from base file changed breast-cancer-rna-analysis.py
- file name to use underscores
- line 8/9 pathing to
  - parser.add_argument("--input", default="data/METABRIC_RNA_Mutation.csv", help="Kaggle breast cancer dataset")
- line 10 pathing to 
  - parser.add_argument("--output", default="figures/breast-cancer-rna_analysis.png")

<div style="margin-top: 50px;"></div>

__folder structure__<br>
library_data_analysis/<br>
├── .gitignore<br>
├── environment.yml<br>
├── README.md<br>
├── data/<br>
├── figures/<br>
└── src/<br>
&emsp;├── breast_cancer_rna_analysis.py<br>
&emsp;├── gen_gene_expr.py<br>
&emsp;└── plot_pkg_stars.py<br>

***

_Need pycodestyle workflow / git action?_

