split_line = "analysis/
print(.wildcards)
"""
    split_list_file = f"analysis/{wildcards.sample}/split/splitList"
    input_fn = open(split_list_file)
    index = 0
    for line in input_fn :
        split_line = line.rstrip("\r\n").split("\t")

        index += 1

        chrom = split_line[0]
        start = int(split_line[1]) + 1
        end   = int(split_line[0])

        config[wildcards.sample][index] = {"CHROM" : chrom, "START" : start, "END" : end}
        print(config)
        return config"""
