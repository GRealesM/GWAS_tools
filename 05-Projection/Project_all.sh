#!/bin/bash

# This script will run fivee Rscripts to project all sumstat files onto their respective bases we have so far.
scriptpath=(~/rds/rds-cew54-basis/GWAS_tools/05-Projection/)

echo "This script will reduce all sumstat files in this directory by the manifests in all three bases"
echo "Projecting for IMD basis..."
echo ""
Rscript "$scriptpath"Projecting_dbs_IMD_basis.R
echo ""
echo "Projecting for cell basis..."
echo ""
Rscript "$scriptpath"Projecting_dbs_cell_basis.R  
echo ""
echo "Projecting for cytokine basis..."
echo ""
Rscript "$scriptpath"Projecting_dbs_cytokine_basis.R
echo "Done!"
echo "Projecting for milmetacytokine basis..."
echo ""
Rscript "$scriptpath"Projecting_dbs_milmetacytokine_basis.R
echo "Done!"
echo "Projecting for cytokine 13 basis..."
echo ""
Rscript "$scriptpath"Projecting_dbs_cytokine13_basis.R
echo "Done!"
