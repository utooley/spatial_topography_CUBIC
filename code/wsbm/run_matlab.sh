
set -euo pipefail
#$ -j y
#$ -o /cbica/projects/spatial_topography/output/job_output/
#$ -l h_vmem=10.5G,s_vmem=10.3G
#$ -V
#$ -cwd

code_dir='/cbica/projects/spatial_topography/code/wsbm'
matlab -nodisplay -r "cd ${code_dir}, run('apply_wsbm_consensus_part_to_test_sample.m'); exit"
