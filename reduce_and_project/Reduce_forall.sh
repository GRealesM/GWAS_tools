#!/bin/bash

# This script will run three Rscripts to filter all sumstat files in this directory by the manifests in the three bases we have so far (IMD, cell and cytokine)


echo "This script will reduce all sumstat files in this directory by the manifests in all three bases"
echo "Reducing for IMD basis..."
echo ""
Rscript Reducing_for_IMDbasis.R
echo ""
echo "Reducing for cell basis..."
echo ""
Rscript Reducing_for_cellbasis.R  
echo ""
echo "Reducing for cytokine basis..."
echo ""
Rscript Reducing_for_cytokinebasis.R
echo "Done!"
