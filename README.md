# ACAM

The source code and results of performance comparison of ACAM against Garnett, CellAssign, SingleR and scmap.

## Notes

- Codes of our method are in the fold [ACAM](./ACAM/). Codes of comparisons are in the fold [compares](./compares/). 
- Since codes are the same among datasets, only the example of the kidney dataset is given. Note that we split Wu dataset randomly into five equal-size subsets. 
- Packages you may need to install and functions we construct to conduct ACAM are in the fold [function](./function/).
- Markers of different forms to cater for the use of different methods. They are given in the fold begin with markers.
- In the fold [results](./results/):
  - 'Y_xxx.raw': The original labels.
  - 'Y_xxx': The numeric form of 'Y_xxx.raw'.
  - 'xxx.comb': Clustering results for checking convenience, since the time cost of clustering may be long.
  - 'umap_xxx': The umap dimensional reduction form of the data.
  - 'cds1_xxx', 'Y_xxx_cellassign', 'Y_xxx_singleR', 'scmap_xxx': comparison results of garnett, cellassign, singler and scmap respectively.
- Consensus tables of cell types and numbers between 'Y_xxx' and 'Y_xxx.raw' are in the fold [celltype](./celltype/) for some dataset needed.
