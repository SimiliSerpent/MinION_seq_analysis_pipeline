include: 'rules/config.py'

# Set constraints on wildcards: only alphanumeric characters are allowed
wildcard_constraints:
    sample='[A-Za-z0-9|_]+',
    species='[A-Za-z0-9]+',
    alignment='[A-Za-z0-9|_]+'


rule all:
    input: expand ('{out_file}', out_file=out_files)

# Trim sequences added by ONT protocol from long reads
rule porechopping_ONT:
    input: data_path + '/raw_fastq/{sample}.fastq'
    output: data_path + '/trimmed_fastq/{sample}_porechopped_ONT.fastq'
    threads: nthreads
    shell:
    	'mkdir -p $(dirname {output}) && '
        '../scripts/Porechops/Porechop_ONT/porechop-runner.py '
    	'--verbosity 1 '
    	'--thread {threads} '
    	'--end_size 200 '
    	'--min_trim_size 4 '
    	'--extra_end_trim 0 '
    	'--min_split_read_size 200 '
    	'-i {input} '
    	'-o {output}'

# Trim Template Switching Oligo (TSO) from long reads
rule porechopping_TSO:
    input: data_path + '/trimmed_fastq/{sample}_porechopped_ONT.fastq'
    output: data_path + '/trimmed_fastq/{sample}_porechopped_ONT_TSO.fastq'
    threads: nthreads
    shell:
        '../scripts/Porechops/Porechop_TSO/porechop-runner.py '
    	'--verbosity 1 '
    	'--thread {threads} '
    	'--end_size 150 '
    	'--min_trim_size 4 '
    	'--extra_end_trim 9 '
    	'--no_split '
    	'-i {input} '
    	'-o {output}'

# Trim poly-A tails
rule trimming_polyA:
    input: data_path + '/trimmed_fastq/{sample}_porechopped_ONT_TSO.fastq'
    output:
        fastq = data_path + '/trimmed_fastq/{sample}_porechopped_ONT_TSO_tail.fastq',
        png = analysis_path + '/tail_trimming/{sample}_tail_lengths_histo.png'
    threads: nthreads
    shell:
        'mkdir -p $(dirname {output.png}) && '
        'python3 -u ../scripts/trimm_polyA.py '
    	'-v 1 '
    	'-i {input} '
    	'-o {output.fastq} '
        '-p {output.png}'

# Filter long GA subsequences
rule filter_GA:
    input: data_path + '/trimmed_fastq/{sample}_porechopped_ONT_TSO_tail.fastq'
    output:
        fastq = data_path + '/trimmed_fastq/{sample}_porechopped_ONT_TSO_tail_GA.fastq',
        png = analysis_path + '/GA_filtering/{sample}_filtered_GA_reps_lengths_histo.png'
    threads: nthreads
    shell:
        'mkdir -p $(dirname {output.png}) && '
        'python3 -u ../scripts/filter_GA.py '
    	'-v 1 '
    	'-i {input} '
    	'-o {output.fastq} '
        '-p {output.png}'

# Compute nucleotides frequencies
rule computing_nuc_freq:
    input: data_path + '/{path_to_sample}.fastq'
    output: analysis_path + '/nuc_freq/{path_to_sample}_freq.txt'
    threads: 1
    shell:
        'mkdir -p $(dirname {output}) && '
        'python3 -u ../scripts/compute_nuc_freq.py -f {input} -o {output}'

# Plot nucleotides frequencies
rule plotting_nuc_freq:
    input: expand(analysis_path+'/nuc_freq/trimmed_fastq/{sample}_porechopped_ONT_TSO_tail_GA_freq.txt', sample=sample_names)
    output: analysis_path + '/nuc_freq/' + exp_id + '_nuc_freq.png'
    threads: 1
    shell:
        'python3 -u ../scripts/plot_nuc_freq.py '
        '-i {analysis_path}/nuc_freq '
        '-o {analysis_path}/nuc_freq/{exp_id}'

# Run FASTQC
rule computing_fastqc:
    input: data_path + '/trimmed_fastq/{sample}.fastq'
    output:
        html = analysis_path + '/fastqc/{sample}_fastqc.html',
        zip = analysis_path + '/fastqc/{sample}_fastqc.zip'
    threads: nthreads
    shell:
        'mkdir -p $(dirname {output.html}) && '
        'fastqc -o $(dirname {output.html}) -t {nthreads} {input}'

# Filtering references using specified regions of interest
rule filtering_reference:
    input:
        ref = PROJECT + '/Data/references/public/{species}.fa',
        roi = PROJECT + '/Data/references/ROI/{species}_regions.txt'
    output: PROJECT + '/Data/references/snakeref/{species}.fa'
    threads: nthreads
    shell:
        'samtools faidx {input.ref} '
        '-r {input.roi} '
        '--threads {nthreads} '
        '> {output}'

# Build reference from species sequences
rule mixing_references:
    input:
        base_ref = PROJECT + '/Data/references/snakeref/{prev_species}.fa',
        next_ref = PROJECT + '/Data/references/snakeref/{next_species}.fa',
        roi = PROJECT + '/Data/references/ROI/{next_species}_regions.txt'
    output: PROJECT + '/Data/references/snakeref/{prev_species}_{next_species}.fa'
    threads: nthreads
    shell:
        'cp {input.base_ref} {input.base_ref}.TEMP && '
        'samtools faidx {input.next_ref} '
        '-r {input.roi} '
        '--threads {nthreads} '
        '>> {input.base_ref}.TEMP && '
        'mv {input.base_ref}.TEMP {output}'

# Get sizes of species regions of interest
rule getting_regions_sizes:
    input:
        ref = PROJECT + '/Data/references/snakeref/{species}.fa',
        roi = PROJECT + '/Data/references/ROI/{species}_regions.txt'
    output: PROJECT + '/Data/references/sizes/{species}_sizes.txt'
    threads: 1
    shell:
        'mkdir -p $(dirname {output}) && '
        'python3 -u ../scripts/get_sizes.py '
        '-f {input.ref} '
        '-o {output} '
        '-r {input.roi}'

# Align reads
rule mapping:
    input:
        fq = data_path + '/trimmed_fastq/{sample}_porechopped_ONT_TSO_tail_GA.fastq',
        ref = PROJECT + '/Data/references/snakeref/' + reference
    output: data_path + '/alignments/{sample}_map2_snakeref.sam'
    threads: nthreads
    shell:
        'mkdir -p $(dirname {output}) && '
        'minimap2 -ax map-ont -t {nthreads} {input.ref} {input.fq} > {output}'

# Sort reads
rule sorting:
    input: data_path + '/alignments/{alignment}.sam'
    output: data_path + '/alignments/{alignment}_sorted.bam'
    threads: nthreads
    shell:
        'samtools sort --threads {nthreads} {input} -o {output}'

# Compute general statistics on alignments
rule computing_general_alignment_stats:
    input: data_path + '/alignments/{alignment}_sorted.bam'
    output: analysis_path + '/alignments_stats/{alignment}.samstats'
    threads: nthreads
    shell:
        'mkdir -p $(dirname {output}) && '
        'samtools stats --threads {nthreads} {input} > {output}'

# Compute target-specific statistics on alignments
rule computing_species_specific_alignment_stats:
    input:
        bam = data_path + '/alignments/{sample}_map2_snakeref_sorted.bam',
        sizes = PROJECT + '/Data/references/sizes/{species}_sizes.txt'
    output: analysis_path + '/alignments_stats/{species}/{sample}.{species}.samstats'
    threads: nthreads
    shell:
        'mkdir -p $(dirname {output}) && '
        'samtools stats --threads {nthreads} {input.bam} -t {input.sizes} '
        '> {output}'

# Extract summary of target-specific alignments stats
rule retrieving_species_stats_summary:
    input: analysis_path + '/alignments_stats/{species}/{sample}.{species}.samstats'
    output: analysis_path + '/alignments_stats/{species}/{sample}.{species}.summary'
    threads: 1
    shell:
        'cat {input} | grep ^SN | cut -f 2- > {output}'

# Retrieve read lengths from alignment statistics
rule retrieving_read_lengths:
    input: analysis_path + '/alignments_stats/{sample}_map2_snakeref.samstats'
    output: analysis_path + '/read_lengths/{sample}.read_lengths'
    threads: 1
    shell:
        'cat {input} | grep ^RL | cut -f 2- > {output}'

# Plot read lengths histogram
rule plotting_read_lengths:
    input: analysis_path + '/read_lengths/{sample}.read_lengths'
    output: analysis_path + '/read_lengths/{sample}_read_lengths_histo.png'
    threads: 1
    shell:
        'python3 -u ../scripts/plot_read_length_histo.py -d {input} -o {output}'

# Plot samples species composition
rule plotting_samples_composition:
    input: expand(analysis_path+'/alignments_stats/{species}/{sample}.{species}.summary', species=species_list, sample=sample_names)
    output:
        analysis_path + '/' + exp_id + '_samples_composition.png',
        analysis_path + '/' + exp_id + '_normalized_samples_composition.png'
    threads: 1
    shell:
        'python3 -u ../scripts/plot_composition.py -i {analysis_path}/alignments_stats -o {analysis_path}/{exp_id}'

# Compute sequences statistics and store in pickle files
rule gathering_sequence_stats:
    input:
        raw = expand(
            data_path+'/raw_fastq/{sample}.fastq',
            sample=sample_names
        ),
        trimmed = expand(
            data_path+'/trimmed_fastq/{sample}_porechopped_{trim_states}.fastq',
            sample=sample_names,
            trim_states=['ONT', 'ONT_TSO', 'ONT_TSO_tail', 'ONT_TSO_tail_GA']
        ),
        sam = expand(
            data_path+'/alignments/{sample}_map2_snakeref.sam',
            sample=sample_names
        )
    output: expand(\
        analysis_path+'/pickles/{sample}_{trim_states}.pkl',\
        sample=sample_names,\
        trim_states=['untrimmed', 'ONT', 'ONT_TSO', 'ONT_TSO_tail', 'ONT_TSO_tail_GA']\
    )
    threads: 10
    shell:
        'mkdir -p {analysis_path}/pickles && '
        'python3 -u ../scripts/compute_sequence_stats.py -e {data_path} -r {PROJECT}/Data/references/ROI -o {analysis_path}/pickles'

# Plot evolution of reads quality vs size across the different trimming states.
# Do so for every reasd, and for randomly sampled fraction of the reads.
rule plotting_seq_size_and_qual_vs_trimming:
    input: expand(\
        analysis_path+'/pickles/{sample}_{trim_states}.pkl',\
        sample=sample_names,\
        trim_states=['untrimmed', 'ONT', 'ONT_TSO', 'ONT_TSO_tail', 'ONT_TSO_tail_GA']\
    )
    output:
        all_reads = analysis_path + '/' + exp_id + '_read_size_quality_with_trimming.png',
        sub_1 = analysis_path + '/' + exp_id + '_read_size_quality_with_trimming_3000.png',
        sub_2 = analysis_path + '/' + exp_id + '_read_size_quality_with_trimming_1500.png',
        sub_3 = analysis_path + '/' + exp_id + '_read_size_quality_with_trimming_5_percent.png'
    threads: 10
    shell:
        'python3 -u ../scripts/plot_sequence_stats.py -p {analysis_path}/pickles -o {output.all_reads} && '
        'python3 -u ../scripts/plot_sequence_stats.py -p {analysis_path}/pickles -o {output.sub_1} -n 3000 && '
        'python3 -u ../scripts/plot_sequence_stats.py -p {analysis_path}/pickles -o {output.sub_2} -n 1500 && '
        'python3 -u ../scripts/plot_sequence_stats.py -p {analysis_path}/pickles -o {output.sub_3} -f 0.05'
