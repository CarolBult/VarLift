# Running VarLift on Windows

VarLift runs on Windows through WSL2 (Windows Subsystem for Linux). Ensembl VEP runs directly in Apptainer inside WSL, so Lima and Homebrew are not needed. VarLift detects the platform and runs Apptainer directly on Linux, and through Lima on macOS.

These steps replace steps 2-5 of the Tool Setup in the [README](README.md). The remaining steps (downloading data and running the web app) are the same.

### Prerequisites
1. WSL2 with Ubuntu. In an administrator PowerShell, run the command below, then restart Windows and open the Ubuntu app to finish setup.

    ```
    wsl --install -d Ubuntu
    ```

2. Apptainer, installed inside Ubuntu.

    ```
    sudo add-apt-repository -y ppa:apptainer/ppa
    sudo apt update
    sudo apt install -y apptainer
    ```

3. Conda, installed inside Ubuntu. [Miniforge](https://github.com/conda-forge/miniforge) works well.

### Tool Setup
Run all commands in the Ubuntu terminal.

1. Clone the repository into your Linux home directory and navigate to it. Keep it in the Linux file system (for example `~/VarLift`), not under `/mnt/c`, which is much slower.

    ```
    cd ~
    git clone https://github.com/cjszem/VarLift
    cd VarLift
    ```

2. Build the Ensembl VEP image in `app/vep`.

    ```
    mkdir -p app/vep/vep_cache
    apptainer build app/vep/ensembl-vep.sif docker://ensemblorg/ensembl-vep
    ```

3. Download the VEP caches. `app/vep` is bound to `/main` inside the container, as VarLift does when it runs VEP.

    ```
    apptainer exec --bind "$PWD/app/vep:/main" app/vep/ensembl-vep.sif perl /opt/vep/src/ensembl-vep/INSTALL.pl -a cf -s homo_sapiens --CACHEDIR /main/vep_cache --DESTDIR /main/vep_cache --ASSEMBLY GRCh38 # download human cache

    apptainer exec --bind "$PWD/app/vep:/main" app/vep/ensembl-vep.sif perl /opt/vep/src/ensembl-vep/INSTALL.pl -a cf -s mus_musculus --CACHEDIR /main/vep_cache --DESTDIR /main/vep_cache --ASSEMBLY GRCm39 # download mouse cache
    ```

4. Create the conda environment. `environment.yml` pins macOS build strings, so remove them first to let conda pick the Linux builds.

    ```
    sed -E 's/^(  - [^=]+=[^=]+)=.*/\1/' environment.yml > environment-linux.yml
    conda env create -f environment-linux.yml
    conda activate VarLift
    ```

Then continue with step 6 of the README (`python app/data_refresh.py`) and the Local Web App Setup. Once uvicorn is running in WSL, open `http://127.0.0.1:8000` in a Windows browser.
