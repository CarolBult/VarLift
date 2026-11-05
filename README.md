# MGI VarLift
The VarLift project is a tool providing a novel pipeline to identify equivalent genome sequence variants in protein-coding genes in mouse and human genomes for efficient mouse model discovery and development.


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

3. Create a Ensembl VEP image inside the apptainer VM.

    VarLift uses Ensembl VEP, which requires a Linux environment to be run. Lima provides the Linux environment and Apptainer creates teh containerized image where vep can be run.

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


### Local Python Use

Create new a python file for your use.

1. Import processing function from `app/processing_util.py`. Import query and fetch functions from `app/main_util.py`.

    ```
    from app.processing_util import prepare_submission
    from app.main_util import hvar_query, mvar_fetch, hvar_query_score, mvar_query, hvar_fetch, mvar_query_score
    ```

2. Prepare variant formatting for pipeline.
   
    ```
    imput = '2:157774114-157774114:C/T\n1:39468726-39468726:T/G''
    variants = prepare_submission(input)
    ```

3. Submit variants to the pipeline.

    ```
    query_gene_df, query_prt_df, input_genes = hvar_query(variants)
    output_gene_df, output_prt_df, output_phenotype_df, input_genes = mvar_fetch(input_genes)
    score_df = hvar_query_score(query_prt_df, output_prt_df, input_genes)
    ```

## Output

1. Input gene table: Provides information on the given variant on a gene level.

    ```
    Gene Symbol Description Biotype Chromosome  Start	End	Strand	Ensembl_ID	Accession
    ACVR1 activin A receptor type 1	protein_coding	2	157736249	157876347	-	ENSG00000115170	HGNC:171
    ```

2. Input variant table: Provides information on the given variant on a transcript level.

    ```
    Gene Symbol	Transcript ID	Biotype	Exon Rank	Pfam Domain ID	Pfam Domain Name	Polyphen Prediction	Polyphen Score	Molecular Consequence	Codon Switch	Amino Acids	Associated Diseases	MONDO
    ACVR1	ENST00000434821	protein_coding	6/11	PF08515	Transforming growth factor beta type I GS-motif	probably_damaging	0.999	missense_variant	cGc/cAc	R206H	fibrodysplasia ossificans progressiva MONDO:0007606
    ```

3. Output gene table: Provides information on the given orthologous variant on a gene level.

    ```
    Gene Symbol	Description	Biotype	Chromosome	Start	End	Strand	Ensembl_ID	Accession
    Acvr1	activin A receptor%2C type 1	protein_coding	2	58278656	58457169	-	ENSMUSG00000026836	MGI:87911
    ```

4. Output variant table: Provides information on the given orthologous variant on a transcript level.

    ```
    Gene Symbol	AlleleID	AlleleSymbol	Transcript ID	Biotype	Exon Rank	Pfam Domain ID	Pfam Domain Name	Molecular Consequence	Codon Switch	Amino Acids	Associated Diseases	MONDO
    Acvr1	MGI:6140231	Acvr1tm1Glh	ENSMUST00000056376	protein_coding	6/11	PF08515	Transforming growth factor beta type I GS-motif	missense_variant	cGC/cAT	R206H	fibrodysplasia ossificans progressiva	MONDO:0007606
    Acvr1	MGI:6140231	Acvr1tm1Glh	ENSMUST00000090935	protein_coding	8/13	PF08515	Transforming growth factor beta type I GS-motif	missense_variant	cGC/cAT	R206H	fibrodysplasia ossificans progressiva	MONDO:0007606
    Acvr1	MGI:6140231	Acvr1tm1Glh	ENSMUST00000112599	protein_coding	5/10	PF08515	Transforming growth factor beta type I GS-motif	missense_variant	cGC/cAT	R206H	fibrodysplasia ossificans progressiva	MONDO:0007606
    Acvr1	MGI:6140231	Acvr1tm1Glh	ENSMUST00000112601	protein_coding	7/12	PF08515	Transforming growth factor beta type I GS-motif	missense_variant	cGC/cAT	R206H	fibrodysplasia ossificans progressiva	MONDO:0007606
    ...
    ```

5. Output score table: Boolean match values for each mouse transcript

    ```
    Allele ID	Allele Symbol	Transcript ID	Biotype Match	Consequence Match	AA Match	AA Position Match	Exon Match	Domain Match	Disease Match	total_score
    MGI:6140231	Acvr1tm1Glh	ENSMUST00000056376	true	true	true	true	true	true	true	100
    MGI:6140231	Acvr1tm1Glh	ENSMUST00000090935	true	true	true	true	false	true	true	85
    MGI:6140231	Acvr1tm1Glh	ENSMUST00000112599	true	true	true	true	false	true	true	85
    MGI:6140231	Acvr1tm1Glh	ENSMUST00000112601	true	true	true	true	false	true	true	85
    ...
    ```

6. Output phenotype table: Provides detailed phenotype information on all alleles in gene.

    ```
    Gene Symbol	AlleleID	AlleleSymbol	Allele Symbol	Phenotypes
    Acvr1	MGI:6140231	Acvr1tm1Glh	Acvr1	MP:0005390,skeleton phenotype
    Acvr1	MGI:1857711	Acvr1tm1Enl	Acvr1	MP:0001672,abnormal embryo development,MP:0001675,abnormal ectoderm development,MP:0001683,absent mesoderm,MP:0001695,abnormal gastrulation,MP:0001698,decreased embryo size,MP:0001710,absent amniotic folds,MP:0002230,abnormal primitive streak formation,MP:0003087,absent allantois,MP:0005030,absent amnion,MP:0009593,absent chorion,MP:0011098,embryonic lethality during organogenesis, complete penetrance,MP:0011186,abnormal visceral endoderm morphology,MP:0011190,thick embryonic epiblast
    ...
    ```

All outputs are temporarily written to corresponding CSV files in the app/results directory. They are downloadable as a complete JSON on the web app.
