#!/bin/bash

# This script will run five Rscripts to filter all sumstat files in this directory by the manifests in the bases we have so far.

scriptpath=(~/rds/rds-cew54-basis/GWAS_tools/04-Reduction/)

echo "This script will reduce all sumstat files in this directory by the manifests in all three bases"
echo "Reducing for IMD basis..."
echo ""
Rscript "$scriptpath"Reducing_for_IMDbasis.R
echo ""
echo "Reducing for cell basis..."
echo ""
Rscript "$scriptpath"Reducing_for_cellbasis.R  
echo ""
echo "Reducing for cytokine basis..."
echo ""
Rscript "$scriptpath"Reducing_for_cytokinebasis.R
echo "Done!"
echo "Reducing for milmetacytokine basis..."
echo ""
Rscript "$scriptpath"Reducing_for_milmetacytokinebasis.R
echo "Done!"
echo "Reducing for cytokine 13 basis..."
echo ""
Rscript "$scriptpath"Reducing_for_cytokine13basis.R
echo "Done!"
