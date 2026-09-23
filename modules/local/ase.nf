process PREPARE_ASE_INPUT {
    tag 'ASE input'
    label 'process_low'

    input:
    path samplesheet
    path hs_counts
    path copy_number
    path priors
    path tx2gene

    output:
    path 'ASE_input.tsv', emit: ase_input

    script:
    """
    prepare_ase_input.py \
        --samplesheet ${samplesheet} \
        --hs-counts ${hs_counts} \
        --copy-number ${copy_number} \
        --priors ${priors} \
        --tx2gene ${tx2gene} \
        --output ASE_input.tsv
    """

    stub:
    """
    printf 'sample\tdna_id\tgroup\tID\tH1_counts\tH2_counts\tNonHS_counts\tCN_H1\tCN_H2\tH1_prior\tH2_prior\nA1\tDNA1\tA\tGENE1\t10\t10\t5\t1\t1\t0.8\t0.8\nB1\tDNA1\tB\tGENE1\t10\t10\t5\t1\t1\t0.8\t0.8\n' > ASE_input.tsv
    """
}

process PREPARE_COMBINED_ASE_INPUT {
    tag 'ASE input'
    label 'process_low'

    input:
    path samplesheet
    path counts
    path priors

    output:
    path 'ASE_input.tsv', emit: ase_input

    script:
    """
    prepare_ase_input.py \
        --samplesheet ${samplesheet} \
        --counts ${counts} \
        --priors ${priors} \
        --output ASE_input.tsv
    """

    stub:
    """
    printf 'sample\tdna_id\tgroup\tID\tH1_counts\tH2_counts\tNonHS_counts\tCN_H1\tCN_H2\tH1_prior\tH2_prior\nA1\tDNA1\tA\tGENE1\t10\t10\t5\t1\t1\t0.8\t0.8\nB1\tDNA1\tB\tGENE1\t10\t10\t5\t1\t1\t0.8\t0.8\n' > ASE_input.tsv
    """
}

process SPLIT_ASE_INPUT {
    tag 'ASE shards'
    label 'process_low'

    input:
    path ase_input
    val contrast
    val shards

    output:
    path 'shards/*.shard_*.tsv', emit: shards
    path 'shards/manifest.tsv', emit: manifest

    script:
    """
    split_ase_input.py \
        --input ${ase_input} \
        --contrast '${contrast}' \
        --shards ${shards} \
        --output-dir shards
    """

    stub:
    """
    # One shard per requested contrast, so a stub run still exercises the fan
    # out that several contrasts produce.
    mkdir -p shards
    printf 'contrast\tgroup_a\tgroup_b\tshard\tgenes\ttotal_counts\tfile\n' > shards/manifest.tsv
    echo '${contrast}' | tr ',' '\n' | while IFS=: read -r group_a group_b; do
        [ -n "\${group_a}" ] || continue
        label="\${group_a}_VS_\${group_b}"
        printf 'sample\tgroup\tID\tH1_counts\tH2_counts\tNonHS_counts\tCN_H1\tCN_H2\tH1_prior\tH2_prior\ttest\tgroupA\tgroupB\n' > "shards/\${label}.shard_0001.tsv"
        printf 'A1\t%s\tGENE1\t10\t10\t5\t1\t1\t0.8\t0.8\t%s\t%s\t%s\n' "\${group_a}" "\${label}" "\${group_a}" "\${group_b}" >> "shards/\${label}.shard_0001.tsv"
        printf 'B1\t%s\tGENE1\t10\t10\t5\t1\t1\t0.8\t0.8\t%s\t%s\t%s\n' "\${group_b}" "\${label}" "\${group_a}" "\${group_b}" >> "shards/\${label}.shard_0001.tsv"
        printf '%s\t%s\t%s\t1\t1\t50\t%s.shard_0001.tsv\n' "\${label}" "\${group_a}" "\${group_b}" "\${label}" >> shards/manifest.tsv
    done
    """
}

process COMPILE_ASE_MODEL {
    tag 'Stan model'
    label 'process_compile'

    input:
    path stan_model

    output:
    path 'ASE_model.rds', emit: model

    script:
    """
    Rscript -e 'suppressPackageStartupMessages(library(rstan)); saveRDS(rstan::stan_model(file="${stan_model}"), "ASE_model.rds")'
    """

    stub:
    """
    touch ASE_model.rds
    """
}

process RUN_ASE_SHARD {
    tag "${shard.baseName}"
    label 'process_ase'

    input:
    path shard
    path compiled_model
    val iterations
    val warmup
    val seed

    output:
    path "${shard.baseName}.results.tsv", emit: results

    script:
    """
    ASE.R \
        --input ${shard} \
        --model ${compiled_model} \
        --output ${shard.baseName}.results.tsv \
        --iterations ${iterations} \
        --warmup ${warmup} \
        --seed ${seed}
    """

    stub:
    """
    printf 'ID\ttest\tgroupA\tgroupB\tgroupA_priorH1\tgroupA_priorH2\tgroupA_totalH1\tgroupA_totalH2\tgroupA_totalNonHS\tgroupA_meanH1\tgroupA_meanH2\tgroupA_meanNonHS\tgroupB_priorH1\tgroupB_priorH2\tgroupB_totalH1\tgroupB_totalH2\tgroupB_totalNonHS\tgroupB_meanH1\tgroupB_meanH2\tgroupB_meanNonHS\tgroupA_alpha_mean\tgroupA_theta_mean\tgroupA_theta_q025\tgroupA_theta_q975\tgroupB_alpha_mean\tgroupB_theta_mean\tgroupB_theta_q025\tgroupB_theta_q975\tgroupA_alphaAI_pvalue\tgroupA_thetaAI_pvalue\tgroupB_alphaAI_pvalue\tgroupB_thetaAI_pvalue\tdelta_theta_mean\tdelta_alpha_mean\tdiffAI_pvalue\trope_value\tanalysis_flag\nGENE1\tA_VS_B\tA\tB\t0.8\t0.8\t10\t10\t5\t10\t10\t5\t0.8\t0.8\t10\t10\t5\t10\t10\t5\t1\t0.5\t0.4\t0.6\t1\t0.5\t0.4\t0.6\t1\t1\t1\t1\t0\t0\t1\t1\tSuccess\n' > ${shard.baseName}.results.tsv
    """
}

process MERGE_ASE_RESULTS {
    tag 'ASE results'
    label 'process_low'

    input:
    path results

    output:
    path 'diffASE_results.tsv', emit: results

    script:
    """
    merge_ase_results.py --inputs ${results} --output diffASE_results.tsv
    """

    stub:
    """
    cp ${results[0]} diffASE_results.tsv
    """
}
