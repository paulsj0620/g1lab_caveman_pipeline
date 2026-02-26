CHR_TO_INDEX = {
    "1": 1, "2": 2, "3": 3, "4": 4, "5": 5,
    "6": 6, "7": 7, "8": 8, "9": 9, "10": 10,
    "11": 11, "12": 12, "13": 13, "14": 14, "15": 15,
    "16": 16, "17": 17, "18": 18, "19": 19, "20": 20,
    "21": 21, "22": 22,
    "X": 23, "Y": 24, "MT": 25
}
rule cvmn_setup:
    input:
        normal_bam=lambda wildcards: config["samples"][wildcards.sample]["normal_bam"],
        tumor_bam=lambda wildcards: config["samples"][wildcards.sample]["tumor_bam"] 
    output:
        cfg="analysis/{sample}/{sample}.cfg.ini",
        rslt=directory("analysis/{sample}/results"),
        algb="analysis/{sample}/{sample}.alg_bean"
    wildcard_constraints:
        sample="[^/]+"
    benchmark: "benchmarks/setup.{sample}.txt"
    log:
        stdout = "logs/cvmn_setup.{sample}.out",
        stderr = "logs/cvmn_setup.{sample}.err"
    threads: 1
    params:
        split="analysis/{sample}/split/splitList"
    resources:          
        tasks=2,        
        mem_mb=100,       
        mpi="srun",     
        runtime=2880,   
        tasks_per_node=1
    container: "docker://quay.io/wtsicgp/dockstore-cgpwgs:2.1.1"
    shell:
        "mkdir -p analysis/{wildcards.sample}/split/ && "
        "caveman setup"
        " -t {input.tumor_bam}"
        " -n {input.normal_bam}"
        " -r {config[ref_genome_idx]}"
        " -g {config[ref_ignore_regions]}"
        " -l {params.split}"
        " -f {output.rslt}"
        " -c {output.cfg}"
        " -a {output.algb}"
        " 1> {log.stdout} 2> {log.stderr}"

rule cvmn_split:
    input:
        cfg="analysis/{sample}/{sample}.cfg.ini"
    output:
        split_fn="analysis/{sample}/split/splitList.{chromosome}"
    wildcard_constraints:
        sample="[^/]+",
        chromosome="[0-9]+|X|Y|MT"
    benchmark: "benchmarks/split.{sample}.{chromosome}.txt"
    log:
        stdout = "logs/cvmn_split.{sample}_{chromosome}.out",
        stderr = "logs/cvmn_split.{sample}_{chromosome}.err"
    threads: 1
    resources:
        tasks=2,
        mem_mb=10000,
        mpi="srun",
        runtime=2880,
        tasks_per_node=1
    container: "docker://quay.io/wtsicgp/dockstore-cgpwgs:2.1.1"
    params : 
        index = lambda wildcards: CHR_TO_INDEX[wildcards.chromosome]
    shell:
        "caveman split"
        " -i {params.index}"
        " -f {input.cfg}"
        " 1> {log.stdout} 2> {log.stderr}"

CHROMOSOMES = list(CHR_TO_INDEX.keys())

checkpoint cvmn_spmg:
    input:
        split_files = expand("analysis/{{sample}}/split/splitList.{chromosome}", 
                             chromosome=CHROMOSOMES)
    output:
        merged = "analysis/{sample}/split/splitList",
        line_count = "analysis/{sample}/split/splitList.line_count"
    wildcard_constraints:
        sample = "[^/]+",
    threads: 1
    log:
        stdout = "logs/cvmn_spmg.{sample}.out",
        stderr = "logs/cvmn_spmg.{sample}.err"
    resources:
        tasks=1,
        mem_mb=4000
    shell:
        "cat {input.split_files} > {output.merged} && "
        "wc -l {output.merged} | awk '{{print $1}}' > {output.line_count}"
        " 2> {log.stderr}"

def get_mstep_jobs(wildcards):
    ckpt_output = checkpoints.cvmn_spmg.get(sample=wildcards.sample).output
    line_count_file = ckpt_output.line_count

    with open(line_count_file) as f:
        num_lines = int(f.read().strip())

    return expand("analysis/{sample}/results/flags/{index}.mstep.done",
                  sample=wildcards.sample,
                  index=range(1, num_lines + 1))

rule cvmn_mstep:
    input:
        cfg = "analysis/{sample}/{sample}.cfg.ini"
    output:
        flag = "analysis/{sample}/results/flags/{index}.mstep.done"
    wildcard_constraints:
        sample = "[^/]+",
        index = r"\d+"
    log:
        "logs/mstep_{sample}_{index}.log"
    benchmark:
        "benchmarks/mstep.{sample}_{index}.txt"
    threads: 1
    container: "docker://quay.io/wtsicgp/dockstore-cgpwgs:2.1.1"
    resources:
        tasks=1,
        mem_mb=4000
    shell:
        """
        caveman mstep -i {wildcards.index} -f {input.cfg} && \
        touch {output.flag}
        """

rule cvmn_merge:
    input:
        cfg = "analysis/{sample}/{sample}.cfg.ini",
        mstep_done = get_mstep_jobs
    output:
        covs = "analysis/{sample}/covs_arr",
        probs = "analysis/{sample}/probs_arr",
    container: "docker://quay.io/wtsicgp/dockstore-cgpwgs:2.1.1"
    wildcard_constraints:
        sample = "[^/]+"
    shell:
        "caveman merge -f {input.cfg} "
        "-c {output.covs} "
        "-p {output.probs}"

def get_estep_jobs(wildcards):
    ckpt_output = checkpoints.cvmn_spmg.get(sample=wildcards.sample).output
    line_count_file = ckpt_output.line_count

    with open(line_count_file) as f:
        num_lines = int(f.read().strip())

    return expand("analysis/{sample}/results/flags/{index}.estep.done",
                  sample=wildcards.sample,
                  index=range(1, num_lines + 1))

rule cvmn_estep:
    input:
        cfg = "analysis/{sample}/{sample}.cfg.ini",
        covs = "analysis/{sample}/covs_arr",
        probs = "analysis/{sample}/probs_arr"
    output:
        flag = "analysis/{sample}/results/flags/{index}.estep.done"
    wildcard_constraints:
        sample = "[^/]+",
        index = r"\d+"
    log:
        "logs/estep_{sample}_{index}.log"
    benchmark:
        "benchmarks/estep.{sample}_{index}.txt"
    threads: 1
    container: "docker://quay.io/wtsicgp/dockstore-cgpwgs:2.1.1"
    resources:
        tasks=1,
        mem_mb=4000
    shell:
        "caveman estep "
        "-i {wildcards.index} "
        "-f {input.cfg} "
        "-v {config[ref_version]} "
        "-w {config[ref_species]} "
        "-g {input.covs} "
        "-o {input.probs} && "
        "touch {output.flag}"

def get_mstep_output(wildcards):
    ckpt_output = checkpoints.cvmn_spmg.get(sample=wildcards.sample).output
    split_list_file = ckpt_output.merged

    regions = []

    with open(split_list_file, 'r') as f:
        for line in f:
            chrom, start, end = line.strip().split('\t')
            
            start_adj = int(start) + 1
            end_adj = int(end)
            regions.append({
                'sample': wildcards.sample,
                'chrom': chrom,
                'start': start_adj,
                'end': end_adj
            })
    
    print(regions)


rule vcf_merge_chr:
    input:
        estep_done = get_estep_jobs
    output:
        muts_vcf = "analysis/{sample}/vcf_merge/{chromosome}/merged.muts.vcf.gz",
        #snvs_vcf = "analysis/{sample}/vcf_merge/{chromosome}/merged.snvs.vcf.gz"
    wildcard_constraints:
        sample="[^/]+",
        chromosome="[0-9]+|X|Y|MT"
    benchmark: "benchmarks/vcf_merge.{chromosome}.{sample}.txt"
    log:
        "logs/vcf_mege.{sample}_{chromosome}.log"
    threads: 1
    resources:
        tasks=2,
        mem_mb=10000,
        mpi="srun",
        runtime=2880,
        tasks_per_node=1
    container: "docker://biocontainers/bcftools:v1.9-1-deb_cv1"
    shell:
        """
        ls -1 analysis/{wildcards.sample}/results/{wildcards.chromosome}/*muts.vcf.gz | sort -V > analysis/{wildcards.sample}/results/{wildcards.chromosome}/muts_files_sorted.txt && \
        ls -1 analysis/{wildcards.sample}/results/{wildcards.chromosome}/*snps.vcf.gz | sort -V > analysis/{wildcards.sample}/results/{wildcards.chromosome}/snvs_files_sorted.txt && \
        bcftools concat -a -O z -f analysis/{wildcards.sample}/results/{wildcards.chromosome}/muts_files_sorted.txt -o analysis/{wildcards.sample}/vcf_merge/{wildcards.chromosome}/merged.muts.vcf.gz && \
        bcftools concat -a -O z -f analysis/{wildcards.sample}/results/{wildcards.chromosome}/snvs_files_sorted.txt -o analysis/{wildcards.sample}/vcf_merge/{wildcards.chromosome}/merged.snvs.vcf.gz && \
        bcftools index analysis/{wildcards.sample}/vcf_merge/{wildcards.chromosome}/merged.muts.vcf.gz && \
        bcftools index analysis/{wildcards.sample}/vcf_merge/{wildcards.chromosome}/merged.snvs.vcf.gz
        """
