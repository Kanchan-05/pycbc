WORKFLOW_NAME=hierarchical
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

# To run the workflow, down the configuration files (*ini, *.py, pycbc* and pycbc_make_coinc_hierarchical_search_workflow) from the gitrepo in your local
# directory. Give the path of this directory in GITLAB_URL
GITLAB_URL="ADD_PATH_TO_THE_LOCAL_DIRECTORY"

ecp-cookie-init LIGO.ORG https://git.ligo.org/users/auth/shibboleth/callback albert.einstein

${LOCAL_PATH}/pycbc_make_coinc_hierarchical_search_workflow \
  --workflow-name ${WORKFLOW_NAME} --output-dir output \
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

