import os
import pyfiglet
from bin.config import Utils
from snakemake.utils import min_version
min_version('6.1.1')

utils = Utils(config)

utils = Utils(config)
utils.add_config(config)
utils.get_paths(config)
utils.get_ref(config) 
utils.print_config(config) 

print("="*100)
print(pyfiglet.figlet_format("G1 Lab Pipeline"))
print(pyfiglet.figlet_format("CaVEMan Automation"))

targets = list()
targets.append('caveman')

rule all:
    input: 
        utils.get_targets(targets),

include: 'modules/caveman.snakefile'
