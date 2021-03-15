"""
Modules containing all configuration variables and paths.
To be used in a development environment.
"""
import os
from _version import __VERSION__, __SCRIPT_NAME__

SNAKEMAKE_VERSION = '5.7.4'  # Snakemake version used to tests the pipeline

# Databases and bedfiles
references = {'ADAPTER_FASTA': "/path/to/adapter.fasta",
              'NON_CG_CUTOFF': "0.2",
              'RESULTS_DIR': "/data/viyt/project/HiTransMeth_result/" + __VERSION__}

# Output directories of OMEN
# Directory containing the semaphores for each Run
semdir = os.path.join(references['RESULTS_DIR'], 'semaphores')
resdir = os.path.join(references['RESULTS_DIR'], 'runs')
calldir = os.path.join(references['RESULTS_DIR'], 'calls')
wkdir = "/scratch/viyt/SPC/project/HiTransMeth"

emails = {'SENDER': "test1@gmail.com",
          'RECIPIENTS': ["test1@gmail.com", "test2@gmail.com"]}