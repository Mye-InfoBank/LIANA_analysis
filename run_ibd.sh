#!/bin/bash

for splitting_key in condition; do
    sbatch submit_liana_analysis.sh \
        "/nfs/data/COST_IBD/downstream_tasks/interactions/output/IBD_v11/sub_atlas_cellxgene_tier1_tier2_removed_colon.h5ad" \
        "cell_type:tier_2" \
        "/nfs/data/COST_IBD/downstream_tasks/interactions/output/IBD_v11/colon/${splitting_key}/" \
        "$splitting_key"
done

for splitting_key in condition; do
    sbatch submit_liana_analysis.sh \
        "/nfs/data/COST_IBD/downstream_tasks/interactions/output/IBD_v11/sub_atlas_cellxgene_tier1_tier2_removed_ileum.h5ad" \
        "cell_type:tier_2" \
        "/nfs/data/COST_IBD/downstream_tasks/interactions/output/IBD_v11/ileum/${splitting_key}/" \
        "$splitting_key"
done