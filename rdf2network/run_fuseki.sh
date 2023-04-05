#!/bin/bash
#SBATCH --job-name=fuseki_server    # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=ana-claudia.sima@sib.swiss     # Where to send mail	
#SBATCH --ntasks=8                    # Run on 8 CPUs
#SBATCH --mem=32gb                     # Job memory request
#SBATCH --time=72:00:00               # Time limit hrs:min:sec
#SBATCH --output=fuseki_log_%j.log   # Standard output and error log
pwd; hostname; date

module load  gcc/10.4.0 jq openjdk/17.0.3_7

export JENA_HOME=/work/FAC/FBM/DBC/cdessim2/default/asima/apache-jena-fuseki-4.7.0

export  PATH=$PATH:$JENA_HOME/

export JVM_ARGS="-Xmx32G"

cd /work/FAC/FBM/DBC/cdessim2/default/asima/string_rdf_unzipped/

fuseki-server --tdb2 --loc=DB_xload/ /string_fuseki
