# MGI VarLift
The VarLift project is a tool providing a novel pipeline to identify equivalent genome sequence variants in protein-coding genes in mouse and human genomes for efficient mouse model discovery and development.

**New to VarLift?** The [User Guide](docs/USER_GUIDE.md) explains how to submit variants, how to read the results, how variant equivalence and phenotype similarity are scored, and which data versions VarLift uses. It is also available in the running app under **User Guide** (`/help`).


## Abstract
VarLift facilitates cross-species identification of equivalent single and multi-nucleotide genetic variants in protein coding sequences of orthologous genes. The tools is intended for two primary groups of users: (1) clinicians with human variant data of unknown phenotype or disease relevance who want to discover if there are equivalent variants in the laboratory mouse and what phenotype or disease annotations are associated with the mouse variants, and (2) researchers who have identified novel sequence variants in protein coding sequences of laboratory mice with abnormal phenotypes and want to determine if an equivalent variant exists in an orthologous human gene. VarLift integrates multiple authoritative resources into a unified pipeline, providing a reproducible framework for cross-species identification of equivalent sequence variants. Variant equivalence is measured by conservation of the variant’s molecular consequence, the position in the protein impacted by the sequence variant, the similarity of the amino acid change, the conservation of the protein domain in which the amino acid change occurs.


## Background
Genome sequencing has become increasingly prevalent in clinical and research settings; however, efficiently connecting identified sequence variants to clinically or biologically meaningful interpretation remains a significant challenge. Mouse Genome Informatics (MGI) is an expertly curated, authoritative resource that facilitates exploration of human genetics through integrated genomic, genetic, and biological data on the laboratory mouse. Despite the growing availability of large-scale sequencing and genetic variant and disease model data, there is currently no efficient approach for directly mapping equivalent human and mouse sequence variants for functional and clinical analysis through disease models. VarLift is a web-based application and software pipeline aimed at filling this gap through facilitating cross-species identification of equivalent sequence variants in orthologous genes through the discovery of disease models.


## Instructions
### Prerequisites
1. macOS, Linux, or Windows: VarLift was created for macOS and also runs on Linux and on Windows through WSL2. On Windows, follow [WINDOWS.md](WINDOWS.md) in place of the Lima steps below.

2. Homebrew: Used to download Lima

3. Conda: Used for package, dependency, and environment management.

### Tool Setup
1. Clone the repository and navigate to the project directory.

    ```
    git clone https://github.com/cjszem/VarLift # clone repo
    cd VarLift # navigate to project
    ```

2. Download Lima and create an Apptainer Linux virtual machine to run ensembl-vep.

    ```
    brew install lima # install lima with homebrew
    limactl start template:apptainer # create an Apptainer VM
    exit # leave limactl shell if inside

    limactl stop apptainer # stop instance
    limactl edit apptainer # edit apptainer configuration - Replace mount section.
    ```

    The mount section should read:
    ```yaml
    mounts:
    - location: "~/VarLift" # Should be absolute path to your VarLift directory
      mountPoint: "/varlift"
      writable: true
    ```

    Restart terminal to ensure the change has taken effect.

3. Create an Ensembl VEP image inside the apptainer VM.

    VarLift uses Ensembl VEP, which requires a Linux environment to be run. Lima provides the Linux environment and Apptainer creates the containerized image where vep can be run.

    ```
    limactl start apptainer # start apptainer instance
    limactl shell apptainer # enter the VM
    cd VarLift # go to repo directory if it is not already working directory
    apptainer build ensembl-vep.sif docker://ensemblorg/ensembl-vep # build a VEP image
    ```

4. Download VEP caches inside vep directory (remain inside of apptainer instance).

    ```
    mkdir -p app/vep/vep_cache # create cache directory

    apptainer exec --bind /varlift:/varlift ensembl-vep.sif perl /opt/vep/src/ensembl-vep/INSTALL.pl -a cf -s homo_sapiens --CACHEDIR /varlift/app/vep/vep_cache --DESTDIR /varlift/app/vep/vep_cache --ASSEMBLY GRCh38 # download human cache

    apptainer exec --bind /varlift:/varlift ensembl-vep.sif   perl /opt/vep/src/ensembl-vep/INSTALL.pl -a cf -s mus_musculus --CACHEDIR /varlift/app/vep/vep_cache --DESTDIR /varlift/app/vep/vep_cache --ASSEMBLY GRCm39 # download mouse cache
    ```

5. Create the conda environment using the provided yml file. All remaining steps occur within this environment.

    ```
    conda env create -f environment.yml # create conda env
    conda activate VarLift # activate conda env
    ```

6. Download all data files.

    ```
    python app/data_refresh.py
    ```

### Local Web App Setup

To use the web app locally follow the steps below:

1. Run the server using uvicorn.

    ```
    uvicorn server:app
    ```

2. Open the webpage at `http://127.0.0.1:8000` in your browser.

3. Insert a list of variants into the search field and submit the query.

### Example Variants

These variants exercise the main paths through VarLift and are useful for checking an installation. Human variants are GRCh38 and mouse variants GRCm39; choose the matching organism before submitting. Several human variants can be submitted together, one per line.

| Variant | Organism | Gene and change | What it shows |
|---|---|---|---|
| `2:157774114-157774114:C/T` | Human | ACVR1 R206H (fibrodysplasia ossificans progressiva) | Three Tier 1 knock-in alleles that model the same disease (Disease Match); one allele numbers the residue R176 on a shorter transcript but still aligns to R206 |
| `X:154030948-154030948:G/A` | Human | MECP2 R306\* (Rett syndrome) | A stop-gained variant: Tier 1 and Tier 2 alleles, and many alleles that are a different consequence (Tier 5) |
| `7:140753336-140753336:A/T` | Human | BRAF V600E | Three Tier 1 alleles; ClinVar lists many cancer conditions |
| `12:112477719-112477719:A/G` | Human | PTPN11 N308D (Noonan syndrome) | One Tier 1 allele among alleles at other residues of the same domain (Tier 3) |
| `X:73080081-73080081:G/A` | Mouse | Mecp2 R133C | Mouse query: nine ClinVar R>C candidates, of which only human R145C is at the aligned residue (Tier 1); pooled gene phenotype score |
| `2:105522966-105522966:G/T` | Mouse | Pax6 splice donor | Splice/Intronic scoring class: Intron and Intron Offset Match |

The web page also offers some of these as clickable examples above the search field.


### Local Python Use

The pipeline can also be run from Python, from the repository root with the conda environment active. These are the analysis steps the web app runs (see `run_pipeline` in `server.py`, which also orders the phenotype table for display).

1. Human variants to equivalent mouse alleles:

    ```python
    from app.processing_util import process_batch_query
    from app.main_util import hvar_query, mvar_fetch, hvar_query_score, phenotype_similarity

    # One variant per line; process_batch_query parses them (a tab and a label may follow each variant)
    variants = process_batch_query('2:157774114-157774114:C/T\n12:112477719-112477719:A/G')
    if variants.empty:
        raise ValueError('No variants recognized; use chromosome:start-end:ref/alt')
    variants = variants.drop(columns='Label')

    query_gene_df, query_prt_df, input_genes, unevaluated_df = hvar_query(variants)
    output_gene_df, output_prt_df, output_phenotype_df, input_genes, allele_coverage_df = mvar_fetch(input_genes)
    score_df = hvar_query_score(query_prt_df, output_prt_df, input_genes, mouse_phenotype_df=output_phenotype_df)

    # One row per evaluated allele: variant and gene phenotype scores, with the term mapping behind each
    phenotype_similarity_df = phenotype_similarity(query_prt_df, score_df[['Input Name', 'AlleleID']],
                                                   output_phenotype_df, detail=True)
    ```

2. Mouse variants to candidate human variants:

    ```python
    from app.processing_util import process_batch_query
    from app.main_util import mvar_query, hvar_fetch, mvar_query_score, gene_phenotype_similarity

    variants = process_batch_query('X:73080081-73080081:G/A')
    if variants.empty:
        raise ValueError('No variants recognized; use chromosome:start-end:ref/alt')
    variants = variants.drop(columns='Label')

    query_gene_df, query_prt_df, input_genes, unevaluated_df, allele_coverage_df = mvar_query(variants)
    output_gene_df, output_prt_df, output_phenotype_df, input_genes = hvar_fetch(query_prt_df, input_genes)
    score_df = mvar_query_score(query_prt_df, output_prt_df, input_genes)
    gene_score_df = gene_phenotype_similarity(allele_coverage_df, input_genes)
    ```

`unevaluated_df` lists submitted variants VarLift could not score (for example, intergenic or non-coding ones).

## Output

The web app returns one JSON object per search, which **Download All Results** saves and **Load Results** reopens. Its tables are:

| Key | Contents |
|---|---|
| `query_genes`, `query_variants` | The submitted variants and their genes, annotated by VEP on one transcript (consequence, exon or intron, Pfam domain, amino-acid change, and for human variants PolyPhen and ClinVar conditions with MONDO IDs) |
| `output_genes`, `output_variants` | The homologous gene and the variants compared with the submitted ones: MGI alleles with coordinates on every transcript (human queries), or ClinVar candidate variants (mouse queries) |
| `scores` | One row per compared variant and transcript: the Equivalency Tier and the criteria behind it (below) |
| `phenotypes` | MGI phenotypes of the compared alleles (human queries), or ClinVar HPO phenotypes of the candidates (mouse queries) |
| `phenotype_similarity` | Phenotype similarity scores with the term-by-term mapping behind each (below) |
| `allele_coverage` | Every MGI allele of the mouse gene, whether it could be evaluated (and why not), its phenotypes and disease models |
| `gene_mapping` | Each submitted variant's gene and its homolog |
| `unevaluated_variants` | Submitted variants that could not be scored |
| `hp_names`, `labels` | Names for the HPO IDs in the results, and the labels entered with the variants |

### Equivalency score columns

`Equivalency Tier` (1–5) ranks how equivalent the two variants are, from Tier 1 (identical change at the aligned residue) to Tier 5 (different consequence), and `Tier Basis` says what it rests on. The `... Match` columns are the evidence: `true`, `false`, or `null` where a criterion does not apply to the variant's `Scoring Class` (Coding or Splice/Intronic). `Aligned Position` is the residue that aligns with the submitted variant's residue in a protein alignment of the two transcripts, and `Local Alignment Identity` the share of identical residues within 10 residues either side. `Disease Match` (human queries) is `true` when the allele models one of the human variant's conditions. `Best Tier` is the variant's best tier on a protein-coding transcript; `Total Score` is kept for older tools and is superseded by the tier. The [User Guide](docs/USER_GUIDE.md#4-how-variant-equivalence-is-scored) defines each criterion and tier.

Example (ACVR1 R206H against Acvr1<sup>tm1Emsh</sup>):

```
Allele Symbol  Transcript ID       Scoring Class  Equivalency Tier  Biotype  Consequence  AA    AA Position  Domain  Aligned Position  Local Alignment Identity  Disease Match
Acvr1tm1Emsh   ENSMUST00000056376  Coding         1                 true     true         true  true         true    206               1.0                       true
```

### Phenotype score columns

Phenotype similarity asks how much of a mouse phenotype a human reference accounts for. Each mouse (MP) term is translated to HPO through curated SSSOM mappings, weighted by mapping type, and matched to the closest human term by semantic similarity; the score is the weighted share, from 0 to 1. It is asymmetric: human terms with no mouse counterpart do not lower it.

| Column | Meaning |
|---|---|
| `Phenotype Score (variant)` | Human queries: the allele's phenotypes against the HPO annotations of each of the human variant's diseases, keeping the best |
| `Matched Disease(s)` | The disease(s) giving that score; tied diseases are joined by `;` |
| `Phenotype Score (gene)` | Human queries: the allele's phenotypes against all HPO terms of the human gene. Mouse queries: the phenotypes of every MGI allele of the mouse gene, pooled, against all HPO terms of the human homolog |
| `MP Terms Mapped`, `MP Terms Total` | How many of the scored MP terms have an HPO translation, of all of them; a low ratio means much of the mouse phenotype cannot be compared |
| `Alleles Pooled`, `Alleles Total`, `Human HP Terms` | Mouse queries: alleles with phenotypes, of all the gene's alleles, and the size of the human reference |
| `Variant Mapping`, `Gene Mapping` | The term-by-term breakdown behind each score |

An empty score means no comparison was possible (no MP term translates, or the human side has no HPO annotations); it is not a score of zero. Section 5 of the [User Guide](docs/USER_GUIDE.md#5-how-phenotype-similarity-is-scored) describes the scoring in full.

Each search also writes its main tables as CSV files to `app/results`, overwritten by the next search.
