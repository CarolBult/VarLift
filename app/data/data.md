# Data Sources for Variant-to-Allele Mapping
## MGI Data

1. data/MGI/VARIANT-ALLELE-MGI.tsv

    Description: MGI transcribed allele entries with genomic positions.

    Source: https://fms.alliancegenome.org/download/VARIANT-ALLELE_NCBITaxon10090.tsv.gz

    Usage: Provides allele-level genomic coordinates for mapping variants to mouse alleles.

2. data/MGI/HOM_ProteinCoding.rpt

    Description: MGI gene homology report.

    Format:
    ```
    MGI_ID    MusGeneSymbol    MusEntrezGeneID    MusHGNC_ID    HumGeneSymbol    HumEntrezGeneID
    ```
    
    Source: https://www.informatics.jax.org/downloads/reports/HOM_ProteinCoding.rpt

    Usage: Links mouse and human genes for cross-species mapping.

3. data/MGI/MGI_DiseaseMouseModel.rpt

    Description: MGI disease models.

    Source: https://www.informatics.jax.org/downloads/reports/MGI_DiseaseMouseModel.rpt

    Usage: Associates mouse models with disease phenotypes.

4. data/MGI/ALL_Phenotype.rpt

    Description: MGI alleles and phenotype information.

    Source: https://www.informatics.jax.org/downloads/reports/ALL_Phenotype.rpt
    
    Usage: Provides allele type and phenotype annotations.


## Disease Ontology

1. data/MONDO/mondo.obo

    Source: https://purl.obolibrary.org/obo/mondo.obo
    
2. data/MONDO/mondo_xref_map.json

    Description: MONDO xref IDs mapped to MONDO IDs.

    Created by: ```data_util.parse_mondo_obo(mondo_path)```

    Usage: Enables synonym resolution for disease terms.

3. data/MONDO/mondo_term_map.json

    Description: MONDO IDs mapped to MONDO terms.

    Created by: ```data_util.parse_mondo_obo(mondo_path)```

    Usage: Enables synonym resolution for disease terms.


## Ensembl Gene Data

1. data/Ensembl/human_gene_data.parquet

    Description: Human gene data extracted from Ensembl GFF3.

    Created by: ```parse_ensembl_gff3('data/Ensembl/Homo_sapiens.GRCh38.116.gff3.gz', species)``` (release from `ensemble_version` in app/config.yaml, via app/data_refresh.py)


2. data/Ensembl/mouse_gene_data.parquet

    Description: Mouse gene data extracted from Ensembl GFF3.

    Created by: ```parse_ensembl_gff3('data/Ensembl/Mus_musculus.GRCm39.116.gff3.gz', species)```



## ClinVar

1. data/NCBI/clinvar.grch38.vcf.gz

    Description: ClinVar entries for GRCh38 variants.

    Source: https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz

    Usage: Provides clinically relevant human variants for mapping.

2. clinvar.grch38.vcf.gz.tbi

    Description: Tabix file for quick indexing of data/NCBI/clinvar.grch38.vcf.gz