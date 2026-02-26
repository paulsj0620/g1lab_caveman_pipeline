import glob
import os, sys

cwd = os.getcwd()
rawdata_path = cwd + '/rawdata'
os.mkdir(rawdata_path)
sample_list = glob.glob("*gz")
print(sample_list)
sample_dict = {}

for sample in sample_list :
    sample_id = sample.split("_R1.fastq.gz")[0]

    sample_path = "{0}/{1}".format(rawdata_path, sample_id)
    os.mkdir(sample_path)

    sample_dict[sample_id] = {}
    os.system("mv {0}*gz {1}/".format(sample_id, sample_path))

files_list = glob.glob("{0}/*/*gz".format(rawdata_path))

output_config_file = open("samples.config.yaml", 'w')
output_config_file.write("samples:" + "\n")

for sample in sample_dict :
    for item in files_list :
        if sample in item :
            if ("_1.fq.gz" in item) or ("R1" in item) :
                fq1 = item
            elif ("_2.fq.gz" in item) or ("R2" in item) : 
                fq2 = item

    sample_dict[sample] = {"fastq1" : fq1, "fastq2" : fq2}

for sample in sample_dict :
    output_config_file.write("  {0}:".format(sample) + "\n")
    output_config_file.write("    fq1: {0}".format(sample_dict[sample]["fastq1"]) + "\n")
    output_config_file.write("    fq2: {0}".format(sample_dict[sample]["fastq2"]) + "\n")

output_config_file.close()

output_cmd_file = open("pipline.cmd.sh", 'w')
output_cmd_file.write('\n')
output_cmd_file.write('wkdir="{0}"\n'.format(cwd))
output_cmd_file.write('conffn="{0}/samples.config.yaml"\n'.format(cwd))
output_cmd_file.write('\n')
output_cmd_file.write("mkdir -p ${wkdir}/logs\n")
output_cmd_file.write('cmd=\"snakemake \\\n')
output_cmd_file.write("    --cores all \\\n")
output_cmd_file.write("    --printshellcmds \\\n")
output_cmd_file.write("    --snakefile /data/SJ/pipeline/onequeue_alignment/Snakefile \\\n")
output_cmd_file.write("    --config 'ref=/data/SJ/pipeline/onequeue_alignment/refs/ref.yaml' \\\n")
output_cmd_file.write("    --max-jobs-per-second 5 \\\n")
output_cmd_file.write("    --jobs 192\n")
output_cmd_file.write("    --configfile ${conffn} \\\n")
output_cmd_file.write("    --directory ${wkdir} \\\n")
output_cmd_file.write("    --executor slurm --default-resources slurm_account=pipeline slurm_partition=cpu --slurm-init-seconds-before-status-checks=60 \\\n")
output_cmd_file.write("    --use-singularity \\\n")
output_cmd_file.write("    --singularity-args '--bind /data,/tmp'\" \n")
output_cmd_file.write("echo $cmd")
