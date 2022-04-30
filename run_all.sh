#!/bin/bash

#===========================================================
# create environment
#===========================================================

. conda_init.sh

./bin/create_conda_envs.sh
./bin/add_conda_pkg.sh

#===========================================================
# alignments
#===========================================================

cd bin/align/
./download_programs.sh
./trim.sh
./create_index.sh
./align.py
cd ../../

#===========================================================
# external data sets
#===========================================================

cd data-raw
./get_aizarani2019.sh
./get_sharma2019.sh
cd ../

#===========================================================
# ref data
#===========================================================

# download ref data
cd bin/ref/
./download_gencode26.sh
Rscript gencode26_info.R
./download_rb_genes.sh
cd ../../

#===========================================================
# Filtering with DIEM
#===========================================================

cd bin/diem/

qsub -sync y run_diem.1.R
qsub -sync y run_diem.2.R
qsub -sync y run_diem.3.R
qsub -sync y run_diem.4.R
qsub -sync y run_diem.match.R
qsub -sync y run_diem.5.R

cd ../../

#===========================================================
# Aizarani clustering
#===========================================================

cd bin/aizarani2019/
qsub -sync y sct_clust.sh
Rscript sct_clust.plot.R
qsub -sync y single_r.ref.sh
Rscript single_r.plot.R
qsub -sync y sct_clust.path_enr.sh
Rscript path_enr.plot.R
Rscript sct_clust.assign_ct.R
Rscript sct_clust.plot_ct.R
cd ../../

#===========================================================
# NASH-HCC clustering
#===========================================================

cd src/d_sct_cca_all

qsub -sync y cca.init.R
qsub -sync y find_markers.R
Rscript plot_clust.init.R
Rscript assign_cell_types.init.R
Rscript npc.1.R
Rscript npc.2.R
Rscript hep.1.R
Rscript hep.2.R
Rscript plot_mrk.npc.R
Rscript plot_mrk.hep.R
Rscript assign_cell_types.npc.R
Rscript assign_cell_types.hep.R
Rscript merge_init.R
qsub -sync y single_r.aiz.sh
Rscript plot_singleR.R # plot singleR results
Rscript plot_clust.merge.R
Rscript reactome_enr.CellType.R
Rscript reactome_enr.CellType.plot.R
Rscript CellTypeProp.R

# For supp figs
Rscript plot_clust.merge.raw.R
Rscript plot_tumor.sfig.R

Rscript split_sample_sct.R
qsub -sync y split_sample_sct.cc.sh

cd ../../

#===========================================================
# Sharma clustering
#===========================================================

cd bin/sharma2019/
qsub -sync y pr_raw_data.sh # process raw data to seurat
qsub -sync y clust_sharma.sh # cluster with no integration
qsub -sync y sct_pat_sharma.sh # run SCT on each patient
qsub -sync y single_r.aiz.sh
cd ../../


#===========================================================
# Sharma and Aizarani integration
#===========================================================

cd bin/sharma_aiz/

qsub -sync y int_datasets.rand.sh
Rscript anno_int.R
Rscript plt_int.R
Rscript plot_mrk.1.R
qsub -sync y find_markers.sh
qsub -sync y path_enr.res1.sh
Rscript path_enr.res1.plot.R
Rscript assign_ct.R
qsub -sync y find_markers.main.sh
Rscript plt_int.2.R # plot cell type assignments
Rscript tum_diff.R
Rscript export_umd.R
qsub -sync y nonprol_ref.sh

# include this one?
# Rscript mrk_dotplot.R

cd ../../

#===========================================================
# TCGA
#===========================================================

cd src/tcga_hcc

# process raw data
Rscript process_data.R
Rscript condense_pheno.R
Rscript norm_counts.R

Rscript pca.R

Rscript sample_map.R # map various IDs

# DE
Rscript de.tcga.R # run genome-wide DE for tumors
Rscript de.mrk.R # check tumor DE of cell-type marker genes
Rscript sc.tum_de.score.R # score droplets for tumor DE gene expression

# gene-based survival tests
Rscript surv.gene.R
Rscript surv.gene.mrk.R
Rscript sc.surv.score.R # droplet scoring of survival genes

# for main cell-types
Rscript main.est_ctp.R
Rscript main.ctp_pheno.R
Rscript main.mrk_coexpr_ctp.R
Rscript main.plot.R
Rscript main.surv.R
Rscript main.surv.plot.R
Rscript main.surv.table.R
Rscript main.mol.R
Rscript main.mol.plot.R # TODO: fix this plotting script

# fine cell-types
Rscript fine.est_ctp.R
Rscript fine.co_abund.plot.R

# mutation-expression association
qsub -sync y mut_expr.sh
Rscript sc.mut.score.R

cd ../../

#===========================================================
# LCI (GSE14520)
#===========================================================

cd src/GSE14520
Rscript get_data.R

Rscript sample_map.R

# DE
Rscript de.LCI.R
Rscript de.mrk.R # check tumor DE of cell-type marker genes
Rscript sc.tum_de.score.R

# gene-based survival tests
Rscript surv.gene.R
Rscript surv.gene.mrk.R
Rscript sc.surv.score.R # droplet scoring of survival genes

# main cell-types
Rscript main.est_ctp.R
Rscript main.ctp_pheno.R
Rscript main.mrk_coexpr_ctp.R
Rscript main.plot.R
Rscript main.surv.R
Rscript main.surv.plot.R
Rscript main.surv.table.R

cd ../../

#===========================================================
# Manuscript items
#===========================================================

cd bin/manuscript
Rscript fig1.R
Rscript fig2.R
Rscript fig3.R
Rscript fig4.R
Rscript supp_fig3.R
Rscript supp_fig4.R
Rscript supp_fig5.R
Rscript supp_fig6.R
Rscript supp_fig7.R
Rscript supp_fig8.R
Rscript supp_fig9.R
python supp_table1.sc_marker.py
python supp_table2.path_enr.py
python supp_table3.bisque.py
python supp_table4.tum_propdiff.py
python supp_table5.surv_propdiff.py
python supp_table6.mut_ctp.py
cd ../../




