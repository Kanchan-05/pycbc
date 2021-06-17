# Two-stage Hierarchical Search Pipeline

This directory contains all the necessary configuration files for running a **Two-stage Hierarchical Pipeline** for a two-ifos. To know the working of this pipeline see: 
- [Soni et al. (2021)](https://arxiv.org/pdf/2106.08925.pdf)
- [Gadre  et al. (2019)](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.99.124035) 

The Hierarchical search pipeline looks for GW signals emitted from non-precessing CBCs with quasicircular orbits. The pipeline search over parameter space via matched filtering in two stages: Stage-1 and Stage-2
1. **Stage-1** searches for the signals using a *Coarse Bank* and data sampled at lower frequency. 
2. **Stage-2** search follows up the coincident triggers from stage-1 and again perform matched filtering analysis over a *Stage-2 Bank* and data sampled at higher frequency.


 > **NOTE**: 
 > - Currently the search is carried out by separately running two workflows one-after-another. First we run a workflow for stage-1 using the configuration files provided in **stage-1_configs**. After the workflow finishes, we run a workflow for stage-2 using the configuration files given in **stage-2_configs** directory. 
 > - The two workflows are executed in two separate directories inside a main local directory (`/home/albert.einstein/search/`) to avoid the confusion. Let's name these directories as: `coarse` for stage-1, and `hierarchical` for stage-2.
 > - All the banks used in the search are downloaded to `/home/albert.einstein/search/`. Let's call this directory as `banks`.
 > - **This setup is a temporary arrangement. We plan to combine both stages together in the future and make it as one complete workflow**

## Instructions to create a workflow for stage-1 search
To perform the stage-1 search, use the workflow `pycbc_make_coinc_search_workflow`, which is also used by [flat-search](https://git.ligo.org/ligo-cbc/pycbc-config/-/blob/master/O3C01/pipelineHL/README.md)

To run this workflow, first: 

- Create a directory `/home/albert.einstein/search/coarse` and download the `*ini` and `run.sh` files from **stage-1_configs** in it. A typical `run.sh` will looks like: 

```
WORKFLOW_NAME=coarse
OBSERVING_RUN='TYPE_OBSERVING_RUN_NUMBER'
DATA='DATA_TYPE'
CONFIG_TAG='1.16.13'
GPSSTART='YOUR_START_TIME'
GPSEND='YOUR_END_TIME'
CHUNKNUMBER='YOUR_ANALYSIS_NUMBER'
DISCRIPTION='INITIAL_OR_RERUN_OR_SOME_DESCRIPTIVE_WORD'

# Source the PyCBC version
source /cvmfs/oasis.opensciencegrid.org/ligo/sw/pycbc/x86_64_rhel_7/virtualenv/pycbc-v1.16.13/bin/activate

# Following paths are specific to Sarathi.
export PATH=/soft/intel/bin/:$PATH
export PATH=/soft/condor_mpi/bin:$PATH
export LD_LIBRARY_PATH=/soft/intel/compilers_and_libraries_2018.0.128/linux/mkl/lib/intel64:$LD_LIBRARY_PATH
export LIGO_DATAFIND_SERVER="ldr.gw.iucaa.in:80"

# To run the workflow, down the configuration files (*ini) from the gitrepo in your local directory.
# Give the path of this directory in GITLAB_URL
GITLAB_URL="ADD_PATH_TO_THE_LOCAL_DIRECTORY"

ecp-cookie-init LIGO.ORG https://git.ligo.org/users/auth/shibboleth/callback albert.einstein

pycbc_make_coinc_search_workflow \
  --workflow-name ${WORKFLOW_NAME} --output-dir output\
  --config-files \
  ${GITLAB_URL}/analysis.ini \
  ${GITLAB_URL}/executables.ini \
  ${GITLAB_URL}/injections_minimal.ini \
  ${GITLAB_URL}/plotting.ini \
  ${GITLAB_URL}/data.ini \
  ${GITLAB_URL}/gating.ini \
  ${GITLAB_URL}/gps_times.ini \
  --config-overrides workflow:start-time:${GPSSTART} workflow:end-time:${GPSEND} \
  'results_page:analysis-subtitle:"O'${OBSERVING_RUN}' Analysis '${WORKFLOW_NAME}' chunk-'${CHUNKNUMBER}', '${DATA}' data"' \
  results_page:output-path:"/home/${USER}/public_html/o${OBSERVING_RUN}/runs/hl/${DATA}/chunk${CHUNKNUMBER}/${WORKFLOW_NAME}/a${CHUNKNUMBER}_${DISCRIPTION}" \
  workflow:file-retention-level:all_triggers
```
- Then execute `./run.sh` to obtain an output directory as `/home/albert.einstein/search/coarse/output`
- And finally submit the workflow by: `pycbc_submit_dax --accounting-group ligo.prod.o1.cbc.bbh.pycbcoffline --dax output/coarse.dax --no-grid`

> **NOTE**: 
> - analysis.ini uses `bank-file`, which has to be the **Coarse Bank**. 
> - The `sample-rate` in analysis.ini is set **512** Hz.


## Instructions to create a workflow for stage-2 search
To perform the stage-2 search, use the workflow `pycbc_make_coinc_search_workflow_hierarchical` provided in the **stage-2_configs** directory.

To run the workflow, follow the steps: 

- Create another directory `/home/albert.einstein/search/hierarchical` and download all the **ini*, **py*, *pycbc**, and `run.sh` files from **stage-2_configs** directory.
- Execute `run.sh` inside this directory, until an output folder is generated as `/home/albert.einstein/search/hierarchical/output`.
- `cd output` and submit the workflow by `pycbc_submit_dax --accounting-group ligo.prod.o1.cbc.bbh.pycbcoffline --dax hierarchical.dax --no-grid`

**NOTE**: 
 - Give the correct paths to bank files (`/home/albert.einstein/search/banks/`) in the analysis.ini for:
    - `bank-file`
    - `connection-file`
    - `coarse-bank`
    - `fine-bank` 
 - Use **2048** as `sample-rate` in analysis.ini
 - `pycbc_make_coinc_search_workflow_hierarchical` workflow uses `pycbc_hierarchical_inspiral` for matched filtering, which modifies some of  the standard python files. These are: 
    - *segment.py*
    - *matched_filter.py*
    - *jobsetup.py*
    - *core.py*
    - *coincidence.py*
  
  Analyst should make sure that the modified files named as *hierarchical_segment.py*, *hierarchical_matched_filter.py*, *hierarchical_jobsetup.py*, *hierarchical_core.py*, *hierarchical_coincidence.py* are downloaded from **stage-2_configs** in the `/home/albert.einstein/search/hierarchical/` before running the workflow. 
 - As described in [Soni et al. (2021)](https://arxiv.org/pdf/2106.08925.pdf), the fit coefficients for the ranking statistics, is generated using stage-1 triggers. Therefore, we generate fit values by running `pycbc_hierarchical_fit_sngls`, and `pycbc_fit_over_param_combine` for each detector. 

### General Comments for stage-2 search


For the time being, the workflow for stage-2 is manually executed by entering right paths to the files from stage-1 in the following ini files. Analyst should provide the correct paths: 
- In analysis.ini: 

```
tmpltbank-pregenerated-bank = /home/albert.einstein/search/banks/bank.hdf
path_full_data = /home/albert.einstein/search/coarse/full_data
path_BNSINJ = /home/albert.einstein/search/coarse/BNSSTT2_INJ_coinc
path_BBHINJ = /home/albert.einstein/search/coarse/BBHSEOBNRV4_INJ_coinc
path_NSBHINJ = /home/albert.einstein/search/coarse/NSBHSEOBNRV4_INJ_coinc
connection-file = /home/albert.einstein/search/banks/connection_coarse_to_flat.hdf
statistic-files = /home/albert.einstein/search/statistic/H1L1-PHASE_TIME_AMP_v1.hdf
```

- In executable.ini:

```
fit_over_param =  /home/albert.einstein/search/hierarchical/pycbc_hierarchical_fit_sngls
fit_over_param_combine = /home/albert.einstein/search/hierarchical/pycbc_fit_over_param_combine
inspiral = /home/albert.einstein/search/hierarchical/pycbc_hierarchical_inspiral
plot_coinc_snrchi = /home/albert.einstein/search/hierarchical/pycbc_page_coinc_snrchi
plot_throughput = /home/albert.einstein/search/hierarchical/pycbc_plot_throughput
```

- In `pycbc_make_coinc_hierarchical_search_workflow`:

```
fitstage1 = ["/home/albert.einstein/search/coarse/output/fit_over_param/H1-FIT_OVER_PARAM-GPSSTART-GPSEND.hdf","/home/albert.einstein/search/coarse/output/fit_over_param/L1-FIT_OVER_PARAM-GPSSTART-GPSEND.hdf"]
```


For further details including run diagnosis and debugging see [LINK](https://pycbc.org/pycbc/latest/html/workflow/pycbc_make_coinc_search_workflow.html#monitor-and-debug-the-workflow-detailed-pegasus-documentation)

## License and Citation
We encourage to cite the paper and references to the codes: 
```
@article{Soni:2021XXX,
    author = "Soni, K. and Gadre, B. U. and Mitra, S. and Dhurandhar, S.",
    title = "{Hierarchical search for compact binary coalescences in the Advanced LIGO’s first two observing runs}",
    eprint = "2106.08925",
    archivePrefix = "arXiv",
    primaryClass = "gr-qc",
    month = "6",
    year = "2021"
}
```
