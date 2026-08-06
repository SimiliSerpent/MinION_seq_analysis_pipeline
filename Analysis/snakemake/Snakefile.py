include: 'rules/config.py'

# Set constraints on wildcards: only alphanumeric characters are allowed
wildcard_constraints:
    barcode='[A-Za-z0-9|_]+',
    species='[A-Za-z0-9]+',
    alignment='[A-Za-z0-9|_]+'


###############################################################
# DEMULTIPLEXING
###############################################################

# Perform demuxing
# WARNING: demuxing with porechop is known as flawed, avoid using at all cost
# TODO: demux using subset of utilized indexes / barcodes
if barcode_type == 'TWIST-LRLP-SHv2E':

    # Demultiplex reads barcoded with Illumina barcodes
    rule demuxing_Illumina:
        input: f'{data_path}/raw_fastq/all_reads.fastq'
        output:
            expand(
                f'{data_path}/demultiplexed/{{barcode}}.fastq',
                barcode=barcode_names
            )
        threads: nthreads
        shell:
            'mkdir -p {data_path}/demultiplexed && '
            'mkdir -p {analysis_path}/logs/demuxing && '
            'python3 -u ../scripts/Porechops/Porechop_modified/porechop-runner.py '
            '--verbosity 1 '
            '--thread {threads} '
            '--end_size 200 '
            '--barcode_threshold 75 '
            '--barcode_diff 5 '
            '--untrimmed '
            '-a Illumina_UDIs '
            '-i {input} '
            '-o {output} '
            '>> {analysis_path}/logs/demuxing/{barcode}_demuxing_Illumina.txt'

elif barcode_type == 'ONT-EXP-PBC001':

    # Demultiplex reads barcoded with ONT barcodes
    rule demuxing_ONT:
        input: f'{data_path}/raw_fastq/all_reads.fastq'
        output:
            expand(
                f'{data_path}/demultiplexed/{{barcode}}.fastq',
                barcode=barcode_names
            )
        threads: nthreads
        shell:
            'mkdir -p {data_path}/demultiplexed && '
            'mkdir -p {analysis_path}/logs/demuxing && '
            'python3 -u ../scripts/Porechops/Porechop_modified/porechop-runner.py '
            '--verbosity 1 '
            '--thread {threads} '
            '--barcode_threshold 75 '
            '--barcode_diff 5 '
            '--require_two_barcodes '
            '--untrimmed '
            '-a ONT_barcodes '
            '-i {input} '
            '-o {output} '
            '>> {analysis_path}/logs/demuxing/{barcode}_demuxing_ONT.txt'

    # Trim sequences added by ONT protocol from long reads
    # rule porechopping_ONT:
    #     input: data_path + '/demultiplexed/{barcode}.fastq'
    #     output: data_path + '/trimmed_fastq/{barcode}_porechopped_ONT.fastq'
    #     threads: nthreads
    #     shell:
    #         'mkdir -p $(dirname {output}) && '
    #         '../scripts/Porechops/Porechop_ONT/porechop-runner.py '
    #         '--verbosity 1 '
    #         '--thread {threads} '
    #         '--end_size 200 '
    #         '--min_trim_size 4 '
    #         '--extra_end_trim 0 '
    #         '--min_split_read_size 200 '
    #         '-i {input} '
    #         '-o {output} '
    #         '>> {analysis_path}/logs/porechop_logs.txt'


###############################################################
# TRIMMING
###############################################################

# Trim ONT ligation adapters from long reads
rule porechopping_sequencing_adapter:
    input: f'{data_path}/demultiplexed/{{barcode}}.fastq'
    output:
        fastq = f'{data_path}/trimmed_fastq/{{barcode}}_porechopped_1_seq_adapt.fastq',
        log = f'{analysis_path}/logs/cleaning/{{barcode}}_porechopping_1_seq_adapt.txt'
    threads: nthreads
    shell:
        'mkdir -p {analysis_path}/logs/cleaning && '
        'if [[ -s "{input}" ]]; then '
        'python3 -u ../scripts/Porechops/Porechop_modified/porechop-runner.py '
        '--verbosity 1 '
        '--thread {threads} '
        '--end_size 100 '
        '--min_trim_size 4 '
        '--extra_end_trim 0 '
        '--min_split_read_size 200 '
        '-a ONT_ligation '
        '-i {input} '
        '-o {output.fastq} '
        '> {output.log}; '
        'else '
        'touch {output.fastq} && touch {output.log}; '
        'fi'

# Trim sequences added during library preparation
if barcode_type == 'TWIST-LRLP-SHv2E':

    # Set var to last Twist trimmed seq so we know what file to use after
    # Twist seq trimming
    last_trimming_before_TSO = '3_twist_primers'
    # Set var to know the trimming steps
    prep_trim_states = ['2_twist_outer', '3_twist_primers']

    # Trim Twist outmost sequences
    rule porechopping_Twist_outer_seq:
        input: f'{data_path}/trimmed_fastq/{{barcode}}_porechopped_1_seq_adapt.fastq'
        output:
            fastq = f'{data_path}/trimmed_fastq/{{barcode}}_porechopped_2_twist_outer.fastq',
            log = f'{analysis_path}/logs/cleaning/{{barcode}}_porechopping_2_twist_outer.txt'
        threads: nthreads
        shell:
            'if [[ -s "{input}" ]]; then '
            'python3 -u ../scripts/Porechops/Porechop_modified/porechop-runner.py '
            '--verbosity 1 '
            '--thread {threads} '
            '--end_size 100 '
            '--min_trim_size 4 '
            '--extra_end_trim 0 '
            '--min_split_read_size 200 '
            '-a Twist_outer_seq '
            '-i {input} '
            '-o {output.fastq} '
            '> {output.log}; '
            'else '
            'touch {output.fastq} && touch {output.log}; '
            'fi'
    
    # TODO: trim Illumina's UDIs using subset of used indexes

    # Trim Twist primers
    rule porechopping_Twist_primers:
        input: f'{data_path}/trimmed_fastq/{{barcode}}_porechopped_2_twist_outer.fastq'
        output:
            fastq = f'{data_path}/trimmed_fastq/{{barcode}}_porechopped_3_twist_primers.fastq',
            log = f'{analysis_path}/logs/cleaning/{{barcode}}_porechopping_3_twist_primers.txt'
        threads: nthreads
        shell:
            'if [[ -s "{input}" ]]; then '
            'python3 -u ../scripts/Porechops/Porechop_modified/porechop-runner.py '
            '--verbosity 1 '
            '--thread {threads} '
            '--end_size 100 '
            '--min_trim_size 4 '
            '--extra_end_trim 0 '
            '--min_split_read_size 200 '
            '-a Twist_primers '
            '-i {input} '
            '-o {output.fastq} '
            '> {output.log}; '
            'else '
            'touch {output.fastq} && touch {output.log}; '
            'fi'

elif barcode_type == 'ONT-EXP-PBC001':

    # Set var to last Twist trimmed seq so we know what file to use after
    # Twist seq trimming
    last_trimming_before_TSO = '3_ont_unknown_seq'
    # Set var to know the trimming steps
    prep_trim_states = ['2_ont_primer_tails', '3_ont_unknown_seq']

    # TODO: trim ONT's barcodes using subset of used barcodes

    # Trim ONT primer tails
    rule porechopping_ONT_primer_tails:
        input: f'{data_path}/trimmed_fastq/{{barcode}}_porechopped_1_seq_adapt.fastq'
        output:
            fastq = f'{data_path}/trimmed_fastq/{{barcode}}_porechopped_2_ont_primer_tails.fastq',
            log = f'{analysis_path}/logs/cleaning/{{barcode}}_porechopping_2_ont_primer_tails.txt'
        threads: nthreads
        shell:
            'if [[ -s "{input}" ]]; then '
            'python3 -u ../scripts/Porechops/Porechop_modified/porechop-runner.py '
            '--verbosity 1 '
            '--thread {threads} '
            '--end_size 100 '
            '--min_trim_size 4 '
            '--extra_end_trim 0 '
            '--min_split_read_size 200 '
            '-a ONT_primer_tail '
            '-i {input} '
            '-o {output.fastq} '
            '> {output.log}; '
            'else '
            'touch {output.fastq} && touch {output.log}; '
            'fi'

    # Trim ONT unknown inner sequences
    rule porechopping_ONT_unknown_seq:
        input: f'{data_path}/trimmed_fastq/{{barcode}}_porechopped_2_ont_primer_tails.fastq'
        output:
            fastq = f'{data_path}/trimmed_fastq/{{barcode}}_porechopped_3_ont_unknown_seq.fastq',
            log = f'{analysis_path}/logs/cleaning/{{barcode}}_porechopping_3_ont_unknown_seq.txt'
        threads: nthreads
        shell:
            'if [[ -s "{input}" ]]; then '
            'python3 -u ../scripts/Porechops/Porechop_modified/porechop-runner.py '
            '--verbosity 1 '
            '--thread {threads} '
            '--end_size 100 '
            '--min_trim_size 4 '
            '--extra_end_trim 0 '
            '--min_split_read_size 200 '
            '-a ONT_unknown_seq '
            '-i {input} '
            '-o {output.fastq} '
            '> {output.log}; '
            'else '
            'touch {output.fastq} && touch {output.log}; '
            'fi'

elif barcode_type == 'UNBARCODED':

    # Set var to last trimmed seq so we know what file to use after
    # Only sequencing adapters trimming if no barcodes
    last_trimming_before_TSO = '1_seq_adapt'
    # Set var to know the trimming steps
    prep_trim_states = []

# Trim Template Switching Oligo (TSO) from long reads
rule porechopping_TSO:
    input: f'{data_path}/trimmed_fastq/{{barcode}}_porechopped_{last_trimming_before_TSO}.fastq'
    output:
        fastq = f'{data_path}/trimmed_fastq/{{barcode}}_porechopped_4_tso.fastq',
        log = f'{analysis_path}/logs/cleaning/{{barcode}}_porechopping_4_tso.txt'
    threads: nthreads
    shell:
        'if [[ -s "{input}" ]]; then '
        '../scripts/Porechops/Porechop_modified/porechop-runner.py '
        '--verbosity 1 '
        '--thread {threads} '
        '--end_size 150 '
        '--min_trim_size 4 '
        '--extra_end_trim 0 '
        '--no_split '
        '-a Takara_TSO '
        '-i {input} '
        '-o {output.fastq} '
        '> {output.log}; '
        'else '
        'touch {output.fastq} && touch {output.log}; '
        'fi'

# Trim poly-A tails
rule trimming_polyA:
    input: f'{data_path}/trimmed_fastq/{{barcode}}_porechopped_4_tso.fastq'
    output:
        fastq = f'{data_path}/trimmed_fastq/{{barcode}}_porechopped_5_polyA.fastq',
        png = f'{analysis_path}/read_cleaning/trimming_stats/polyA_trimming/{{barcode}}_polyA_lengths_histo.png',
        json = f'{analysis_path}/read_cleaning/trimming_stats/polyA_trimming/{{barcode}}_polyA_trimming_stats.json'
    threads: nthreads
    shell:
        'mkdir -p $(dirname {output.png}) && '
        'python3 -u ../scripts/pipeline/trimm_polyA.py '
        '-v 1 '
        '-i {input} '
        '-o {output.fastq} '
        '-p {output.png} '
        '-s {output.json}'

# Filter long GA subsequences
rule filter_GA:
    input: f'{data_path}/trimmed_fastq/{{barcode}}_porechopped_5_polyA.fastq'
    output:
        fastq = f'{data_path}/trimmed_fastq/{{barcode}}_porechopped_6_GA.fastq',
        png = f'{analysis_path}/read_cleaning/GA_filtering/{{barcode}}_filtered_GA_reps_lengths_histo.png',
        json = f'{analysis_path}/read_cleaning/GA_filtering/{{barcode}}_GA_filtering_stats.json'
    threads: nthreads
    shell:
        'mkdir -p $(dirname {output.png}) && '
        'python3 -u ../scripts/pipeline/filter_GA.py '
        '-v 1 '
        '-i {input} '
        '-o {output.fastq} '
        '-p {output.png} '
        '-s {output.json}'

# Filter escape sequences in porechop output (one file by porechopping operation)
rule filter_single_porechop_log:
    input: f'{analysis_path}/logs/cleaning/{{barcode}}_porechopping_{{trim_state}}.txt'
    output: f'{analysis_path}/logs/cleaning/filtered/{{barcode}}_porechopping_{{trim_state}}.txt'
    threads: 3
    shell:
        'mkdir -p $(dirname {output}) && '
        'python3 -u ../scripts/pipeline/filter_escape_seq.py '
        '-i {input} '
        '-o {output}'

# Extract porechopping statistics
rule extract_porechop_stats:
    input: 
        expand(
            f'{analysis_path}/logs/cleaning/filtered/{{barcode}}_porechopping_{{trim_state}}.txt',
            barcode=barcode_names,
            trim_state=['1_seq_adapt', '4_tso']+prep_trim_states
        )
    output: f'{analysis_path}/read_cleaning/trimming_stats/porechopping_stats.json'
    threads: 1
    shell:
        'mkdir -p $(dirname {output}) && '
        'python3 -u ../scripts/pipeline/extract_porechop_trim_stats.py '
        '-l {analysis_path}/logs/cleaning/filtered/ '
        '-o {output}'

# Plot the reads cleaning statistics for all barcodes
rule plot_barcodes_cleaning_stats:
    input:
        porechop = f'{analysis_path}/read_cleaning/trimming_stats/porechopping_stats.json',
        tail = expand(
            f'{analysis_path}/read_cleaning/trimming_stats/polyA_trimming/{{barcode}}_polyA_trimming_stats.json',
            barcode=barcode_names
        ),
        GA = expand(
            f'{analysis_path}/read_cleaning/GA_filtering/{{barcode}}_GA_filtering_stats.json',
            barcode=barcode_names
        )
    output: f'{analysis_path}/read_cleaning/barcodes_cleaning_stats.png'
    threads: 3
    shell:
        'python3 -u ../scripts/pipeline/plot_trim_stats.py '
        '-i {analysis_path} '
        '-o {output}'

# Plot the reads cleaning statistics for every used barcode
rule plot_samples_cleaning_stats:
    input:
        porechop = f'{analysis_path}/read_cleaning/trimming_stats/porechopping_stats.json',
        tail = expand(
            f'{analysis_path}/read_cleaning/trimming_stats/polyA_trimming/{{barcode}}_polyA_trimming_stats.json',
            barcode=barcode_names
        ),
        GA = expand(
            f'{analysis_path}/read_cleaning/GA_filtering/{{barcode}}_GA_filtering_stats.json',
            barcode=barcode_names
        )
    output: f'{analysis_path}/read_cleaning/samples_cleaning_stats.png'
    threads: 3
    shell:
        'python3 -u ../scripts/pipeline/plot_trim_stats.py '
        '-i {analysis_path} '
        '-s {data_path} '
        '-o {output}'


###############################################################
# COMPUTING SIMPLE STATS
###############################################################

# Compute nucleotides frequencies
rule computing_nuc_freq:
    input: f'{data_path}/{{path_to_barcode}}.fastq'
    output: f'{analysis_path}/nuc_freq/{{path_to_barcode}}_freq.txt'
    threads: 1
    shell:
        'mkdir -p $(dirname {output}) && '
        'python3 -u ../scripts/pipeline/compute_nuc_freq.py '
        '-f {input} '
        '-o {output}'

# Plot nucleotides frequencies
rule plotting_nuc_freq:
    input:
        expand(
            f'{analysis_path}/nuc_freq/trimmed_fastq/{{barcode}}_porechopped_6_GA_freq.txt',
            barcode=barcode_names
        )
    output: f'{analysis_path}/nuc_freq/{exp_id}_nuc_freq.png'
    threads: 1
    shell:
        'python3 -u ../scripts/pipeline/plot_nuc_freq.py '
        '-i {analysis_path}/nuc_freq '
        '-o {analysis_path}/nuc_freq/{exp_id}'

# Run FASTQC
rule computing_fastqc:
    input: f'{data_path}/trimmed_fastq/{{barcode}}.fastq'
    output:
        html = f'{analysis_path}/fastqc/{{barcode}}_fastqc.html',
        zip = f'{analysis_path}/fastqc/{{barcode}}_fastqc.zip'
    threads: nthreads
    shell:
        'mkdir -p $(dirname {output.html}) && '
        'fastqc -o $(dirname {output.html}) -t {nthreads} {input}'


###############################################################
# MAPPING
###############################################################

# Filtering references using specified regions of interest
rule filtering_reference:
    input:
        ref = f'{project_dir}/Data/references/public/{{species}}.fa',
        roi = f'{project_dir}/Data/references/ROI/{{species}}_regions.txt'
    output: f'{project_dir}/Data/references/snakeref/{{species}}.fa'
    threads: nthreads
    shell:
        'samtools faidx {input.ref} '
        '-r {input.roi} '
        '--threads {nthreads} '
        '> {output}'

# Build reference from species sequences
rule mixing_references:
    input:
        base_ref = f'{project_dir}/Data/references/snakeref/{{prev_species}}.fa',
        next_ref = f'{project_dir}/Data/references/snakeref/{{next_species}}.fa',
        roi = f'{project_dir}/Data/references/ROI/{{next_species}}_regions.txt'
    output: f'{project_dir}/Data/references/snakeref/{{prev_species}}_{{next_species}}.fa'
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
        ref = f'{project_dir}/Data/references/snakeref/{{species}}.fa',
        roi = f'{project_dir}/Data/references/ROI/{{species}}_regions.txt'
    output: f'{project_dir}/Data/references/sizes/{{species}}_sizes.txt'
    threads: 1
    shell:
        'mkdir -p $(dirname {output}) && '
        'python3 -u ../scripts/pipeline/get_sizes.py '
        '-f {input.ref} '
        '-o {output} '
        '-r {input.roi}'

# Align reads
rule mapping:
    input:
        fq = f'{data_path}/trimmed_fastq/{{barcode}}_porechopped_6_GA.fastq',
        ref = f'{project_dir}/Data/references/snakeref/{reference}'
    output: f'{data_path}/alignments/{{barcode}}_map2_snakeref.sam'
    threads: nthreads
    shell:
        'mkdir -p $(dirname {output}) && '
        'minimap2 -ax map-ont -t {nthreads} {input.ref} {input.fq} > {output}'

# Sort reads
rule sorting:
    input: f'{data_path}/alignments/{{alignment}}.sam'
    output: f'{data_path}/alignments/{{alignment}}_sorted.bam'
    threads: nthreads
    shell:
        'samtools sort --threads {nthreads} {input} -o {output}'

# Compute general statistics on alignments
rule computing_general_alignment_stats:
    input: f'{data_path}/alignments/{{alignment}}_sorted.bam'
    output: f'{analysis_path}/alignments_stats/{{alignment}}.samstats'
    threads: nthreads
    shell:
        'mkdir -p $(dirname {output}) && '
        'samtools stats --threads {nthreads} {input} > {output}'

# Compute target-specific statistics on alignments
rule computing_species_specific_alignment_stats:
    input:
        bam = f'{data_path}/alignments/{{barcode}}_map2_snakeref_sorted.bam',
        sizes = f'{project_dir}/Data/references/sizes/{{species}}_sizes.txt'
    output: f'{analysis_path}/alignments_stats/{{species}}/{{barcode}}.{{species}}.samstats'
    threads: nthreads
    shell:
        'mkdir -p $(dirname {output}) && '
        'samtools stats --threads {nthreads} {input.bam} -t {input.sizes} '
        '> {output}'

# Extract summary of target-specific alignments stats
rule retrieving_species_stats_summary:
    input: f'{analysis_path}/alignments_stats/{{species}}/{{barcode}}.{{species}}.samstats'
    output: f'{analysis_path}/alignments_stats/{{species}}/{{barcode}}.{{species}}.summary'
    threads: 1
    shell:
        'cat {input} | grep ^SN | cut -f 2- > {output}'

# Retrieve read lengths from alignment statistics
rule retrieving_read_lengths:
    input: f'{analysis_path}/alignments_stats/{{barcode}}_map2_snakeref.samstats'
    output: f'{analysis_path}/read_lengths/{{barcode}}.read_lengths'
    threads: 1
    shell:
        # Retrieve the lines corresponding to read lengths distribution in sam
        'grep -E ^RL {input} > {output}.tmp || true; '
        # If any (= not 0 reads), remove the leading column
        'if [[ -s "{output}.tmp" ]]; then '
        'cat {output}.tmp | cut -f 2- > {output} && rm {output}.tmp; '
        # Else, leave an empty stat file
        'else mv {output}.tmp {output}; '
        'fi'

# Plot read lengths histogram
rule plotting_read_lengths:
    input: f'{analysis_path}/read_lengths/{{barcode}}.read_lengths'
    output: f'{analysis_path}/read_lengths/{{barcode}}_read_lengths_histo.png'
    threads: 1
    shell:
        'python3 -u ../scripts/pipeline/plot_read_length_histo.py '
        '-d {input} '
        '-o {output}'

# Plot barcodes species composition
rule plotting_barcodes_composition:
    input:
        expand(
            f'{analysis_path}/alignments_stats/{{species}}/{{barcode}}.{{species}}.summary',
            species=species_list,
            barcode=barcode_names
        )
    output:
        f'{analysis_path}/composition/{exp_id}_barcodes_composition.png',
        f'{analysis_path}/composition/{exp_id}_barcodes_composition_with_total.png',
        f'{analysis_path}/composition/{exp_id}_barcodes_normalized_composition.png',
        f'{analysis_path}/composition/{exp_id}_barcodes_normalized_composition_with_total.png'
    threads: 1
    shell:
        'mkdir -p {analysis_path}/composition && '
        'python3 '
        '-u ../scripts/pipeline/plot_composition.py '
        '-i {analysis_path}/alignments_stats '
        '-o {analysis_path}/composition/{exp_id}'

# Plot samples species composition
rule plotting_samples_composition:
    input:
        expand(
            f'{analysis_path}/alignments_stats/{{species}}/{{barcode}}.{{species}}.summary',
            species=species_list,
            barcode=barcode_names
        )
    output:
        f'{analysis_path}/composition/{exp_id}_samples_composition.png',
        f'{analysis_path}/composition/{exp_id}_samples_composition_with_total.png',
        f'{analysis_path}/composition/{exp_id}_samples_normalized_composition.png',
        f'{analysis_path}/composition/{exp_id}_samples_normalized_composition_with_total.png'
    threads: 1
    shell:
        'mkdir -p {analysis_path}/composition && '
        'python3 '
        '-u ../scripts/pipeline/plot_composition.py '
        '-i {analysis_path}/alignments_stats '
        '-o {analysis_path}/composition/{exp_id} '
        '-s {samples_str}'

# Compute sequences statistics and store in pickle files
rule gathering_sequence_stats:
    input:
        demultiplexed = expand(
            f'{data_path}/demultiplexed/{{barcode}}.fastq',
            barcode=barcode_names
        ),
        trimmed = expand(
            f'{data_path}/trimmed_fastq/{{barcode}}_porechopped_{{trim_states}}.fastq',
            barcode=barcode_names,
            trim_states=['1_seq_adapt', '4_tso', '5_polyA', '6_GA']+prep_trim_states
        ),
        sam = expand(
            f'{data_path}/alignments/{{barcode}}_map2_snakeref.sam',
            barcode=barcode_names
        )
    output: 
        expand(\
            f'{analysis_path}/read_cleaning/pickles/{{barcode}}_{{trim_states}}.pkl',
            barcode=barcode_names,
            trim_states=['untrimmed', '1_seq_adapt', '4_tso', '5_polyA', '6_GA']+prep_trim_states
        )
    threads: 10
    shell:
        'mkdir -p {analysis_path}/read_cleaning/pickles && '
        'python3 -u ../scripts/pipeline/compute_sequence_stats.py '
        '-e {data_path} '
        '-r {project_dir}/Data/references/ROI '
        '-o {analysis_path}/read_cleaning/pickles'

# Plot evolution of reads quality vs size across the different trimming states.
# Do so for every read, and for randomly sampled fraction of the reads.
rule plotting_barcode_size_vs_quality:
    input:
        expand(
            f'{analysis_path}/read_cleaning/pickles/{{barcode}}_{{trim_states}}.pkl',
            barcode=barcode_names,
            trim_states=['untrimmed', '1_seq_adapt', '4_tso', '5_polyA', '6_GA']+prep_trim_states
        )
    output:
        all_reads = f'{analysis_path}/size_vs_quality/{exp_id}_barcode_size_vs_quality.png',
        sub_1 = f'{analysis_path}/size_vs_quality/{exp_id}_barcode_size_vs_quality_3000.png',
        sub_2 = f'{analysis_path}/size_vs_quality/{exp_id}_barcode_size_vs_quality_1500.png',
        sub_3 = f'{analysis_path}/size_vs_quality/{exp_id}_barcode_size_vs_quality_5_percent.png'
    threads: 10
    shell:
        'mkdir -p {analysis_path}/size_vs_quality && '
        'python3 -u ../scripts/pipeline/plot_sequence_stats.py '
        '-p {analysis_path}/read_cleaning/pickles '
        '-o {output.all_reads} && '
        'python3 -u ../scripts/pipeline/plot_sequence_stats.py '
        '-p {analysis_path}/read_cleaning/pickles '
        '-o {output.sub_1} '
        '-n 3000 && '
        'python3 -u ../scripts/pipeline/plot_sequence_stats.py '
        '-p {analysis_path}/read_cleaning/pickles '
        '-o {output.sub_2} '
        '-n 1500 && '
        'python3 -u ../scripts/pipeline/plot_sequence_stats.py '
        '-p {analysis_path}/read_cleaning/pickles '
        '-o {output.sub_3} '
        '-f 0.05'

# Same as above for named samples only.
rule plotting_sample_size_vs_quality:
    input:
        expand(
            f'{analysis_path}/read_cleaning/pickles/{{barcode}}_{{trim_states}}.pkl',
            barcode=barcode_names,
            trim_states=['untrimmed', '1_seq_adapt', '4_tso', '5_polyA', '6_GA']+prep_trim_states
        )
    output:
        all_reads = f'{analysis_path}/size_vs_quality/{exp_id}_sample_size_vs_quality.png',
        sub_1 = f'{analysis_path}/size_vs_quality/{exp_id}_sample_size_vs_quality_5_percent.png'
    threads: 10
    shell:
        'mkdir -p {analysis_path}/size_vs_quality && '
        'python3 -u ../scripts/pipeline/plot_sequence_stats.py '
        '-p {analysis_path}/read_cleaning/pickles '
        '-o {output.all_reads} '
        '-s {data_path} && '
        'python3 -u ../scripts/pipeline/plot_sequence_stats.py '
        '-p {analysis_path}/read_cleaning/pickles '
        '-o {output.sub_1} '
        '-s {data_path} '
        '-f 0.05'

# Plot evolution of reads quality vs size across the different species reads.
rule plotting_heatmaps:
    input:
        expand(
            f'{analysis_path}/read_cleaning/pickles/{{barcode}}_{{trim_states}}.pkl',
            barcode=barcode_names,
            trim_states=['untrimmed', '1_seq_adapt', '4_tso', '5_polyA', '6_GA']+prep_trim_states
        )
    output:
        barcodes = f'{analysis_path}/size_vs_quality/{exp_id}_barcode_heatmap.png',
        samples = f'{analysis_path}/size_vs_quality/{exp_id}_sample_heatmap.png'
    threads: 10
    shell:
        'mkdir -p {analysis_path}/size_vs_quality && '
        'python3 -u ../scripts/pipeline/plot_sequence_stats_heatmap.py '
        '-p {analysis_path}/read_cleaning/pickles '
        '-o {output.barcodes} && '
        'python3 -u ../scripts/pipeline/plot_sequence_stats_heatmap.py '
        '-p {analysis_path}/read_cleaning/pickles '
        '-o {output.samples} '
        '-s {data_path}'


###############################################################
# UMI DEDUPLICATION
###############################################################

# Filter out reads without a TSO trimmed (totally or partially).
# (other reads are not deemed sufficiently normal for UMI search)
rule filter_TSO_trimmed:
    input:
        pre_tso_trim = f'{data_path}/trimmed_fastq/{{barcode}}_porechopped_{last_trimming_before_TSO}.fastq',
        post_tso_trim = f'{data_path}/trimmed_fastq/{{barcode}}_porechopped_4_tso.fastq'
    output: f'{data_path}/umi_deduplication/{{barcode}}/{{barcode}}_00_tso-trimmed.fastq',
    log: f'{analysis_path}/logs/umi_deduplication/{{barcode}}_umi_dedup_log.txt'
    threads: 10
    shell:
        'mkdir -p {data_path}/umi_deduplication/{wildcards.barcode} && '
        'mkdir -p $(dirname {log}) && '
        'echo ">>>>> FILTERING READS WITH TSO TRIMMED <<<<<" > {log} && '
        'python3 -u ../scripts/pipeline/pick_touched_reads.py '
        '-i {input.pre_tso_trim} {input.post_tso_trim} '
        '-m {output} '
        '-v 3 '
        '>> {log}'

# Filter long GA repeats in reads with TSO trimmed
rule filter_GA_in_TSO_trimmed:
    input: f'{data_path}/umi_deduplication/{{barcode}}/{{barcode}}_00_tso-trimmed.fastq'
    output:
        fastq = f'{data_path}/umi_deduplication/{{barcode}}/{{barcode}}_01_GA-filtered.fastq',
        histo = f'{analysis_path}/umi_stats/{{barcode}}/{{barcode}}_GA_filt_histo.png',
        stats = f'{analysis_path}/umi_stats/{{barcode}}/{{barcode}}_GA_filt_stats.json'
    log: f'{analysis_path}/logs/umi_deduplication/{{barcode}}_umi_dedup_log.txt'
    threads: 10
    shell:
        'echo ">>>>> FILTERING GA REPEATS <<<<<" >> {log} && '
        'python3 -u ../scripts/pipeline/filter_GA.py '
        '-i {input} '
        '-o {output.fastq} '
        '-p {output.histo} '
        '-s {output.stats} '
        '-v 3 '
        '>> {log}'
    
# Remove reads shorter than UMI length
rule remove_short_reads:
    input: f'{data_path}/umi_deduplication/{{barcode}}/{{barcode}}_01_GA-filtered.fastq'
    output: f'{data_path}/umi_deduplication/{{barcode}}/{{barcode}}_02_longer-than-UMI.fastq'
    threads: 10
    shell:
        'seqkit seq -m {umi_length} {input} > {output}'

# Map reads to reference
rule mapping_TSO_trimmed:
    input:
        ref = f'{project_dir}/Data/references/snakeref/{reference}',
        fq = f'{data_path}/umi_deduplication/{{barcode}}/{{barcode}}_02_longer-than-UMI.fastq'
    output: f'{data_path}/umi_deduplication/{{barcode}}/{{barcode}}_03_mapped.sam'
    threads: nthreads
    shell:
        'minimap2 -ax map-ont -t {nthreads} {input.ref} {input.fq} > {output}'
    
# Sort alignment
rule sorting_TSO_trimmed:
    input: f'{data_path}/umi_deduplication/{{barcode}}/{{barcode}}_03_mapped.sam'
    output: f'{data_path}/umi_deduplication/{{barcode}}/{{barcode}}_04_sorted.bam'
    threads: nthreads
    shell:
        'samtools sort --threads {nthreads} {input} -o {output}'

# Filter primary alignment
# TODO: Reconsider using only primary aligned BAM records
rule filter_primary_aligned:
    input: f'{data_path}/umi_deduplication/{{barcode}}/{{barcode}}_04_sorted.bam'
    output: f'{data_path}/umi_deduplication/{{barcode}}/{{barcode}}_05_primary.bam'
    threads: nthreads
    shell:
        'samtools view --threads {nthreads} -b -F 0x900 {input} > {output}'

# Extract raw UMIs
checkpoint extract_raw_umis:
    input: f'{data_path}/umi_deduplication/{{barcode}}/{{barcode}}_05_primary.bam'
    output:
        bam = f'{data_path}/umi_deduplication/{{barcode}}/{{barcode}}_06_no-umis.bam',
        raw_umis = f'{analysis_path}/umi_stats/{{barcode}}/extracted_raw_umis.json',
        raw_umis_stats = f'{analysis_path}/umi_stats/{{barcode}}/extracted_raw_umis.tsv'
    log: f'{analysis_path}/logs/umi_deduplication/{{barcode}}_umi_dedup_log.txt'
    threads: 10
    shell:
        'mkdir -p {analysis_path}/umi_stats/{wildcards.barcode} && '
        'echo ">>>>> EXTRACT RAW UMIS <<<<<" >> {log} && '
        'python3 -u ../scripts/pipeline/extract_umis.py '
        '-b {input} '
        '-p {umi_pattern} '
        '-o {analysis_path}/umi_stats/{wildcards.barcode}/extracted '
        '-ob {output.bam} '
        '-v 3 '
        '>> {log}'

# Make fastq of UMIs-extracted reads
# (use checkpoint here to prevent later errors caused by sample with no reads
# mapping to the reference after TSO-trimming)
rule make_UMIs_less_fastq:
    input: f'{data_path}/umi_deduplication/{{barcode}}/{{barcode}}_06_no-umis.bam'
    output: f'{data_path}/umi_deduplication/{{barcode}}/{{barcode}}_07_no-umis.fastq'
    threads: 10
    shell:
        'samtools fastq --threads 10 {input} > {output}'

# Define helper function to modify expected files list accordingly
def umi_survivors(wildcards):
    barcodes_mapping_ref = []
    for name in barcode_names:
        bam = checkpoints.extract_raw_umis.get(barcode=name).output[0]
        if os.path.getsize(bam) > 0 and \
           int(subprocess.run(['samtools', 'view', '-c', '-F', '0x900', bam],
                              capture_output=True, text=True).stdout or 0) > 0:
            barcodes_mapping_ref.append(name)
    return barcodes_mapping_ref

# Map reads to reference
rule mapping_UMIs_less_reads:
    input:
        ref = f'{project_dir}/Data/references/snakeref/{reference}',
        fq = f'{data_path}/umi_deduplication/{{barcode}}/{{barcode}}_07_no-umis.fastq'
    output: f'{data_path}/umi_deduplication/{{barcode}}/{{barcode}}_08_mapped.sam'
    threads: nthreads
    shell:
        'minimap2 -ax map-ont -t {nthreads} {input.ref} {input.fq} > {output}'

# Sort alignment
rule sorting_UMIs_less:
    input: f'{data_path}/umi_deduplication/{{barcode}}/{{barcode}}_08_mapped.sam'
    output: f'{data_path}/umi_deduplication/{{barcode}}/{{barcode}}_09_no-UMIs.bam'
    threads: nthreads
    shell:
        'samtools sort --threads {nthreads} {input} -o {output}'
    
# Gather reads statistics from BAM
rule get_UMIs_less_reads_stats:
    input: f'{data_path}/umi_deduplication/{{barcode}}/{{barcode}}_09_no-UMIs.bam'
    output: f'{analysis_path}/umi_stats/{{barcode}}/reads_stats.json'
    log: f'{analysis_path}/logs/umi_deduplication/{{barcode}}_umi_dedup_log.txt'
    threads: 10
    shell:
        'echo ">>>>> GET ALIGNMENTS STATS <<<<<" >> {log} && '
        'python3 -u ../scripts/pipeline/extract_bam_stats.py '
        '-b {input} '
        '-o {output} '
        '-v 3 '
        '>> {log}'

# Split raw UMIs into connected read sets
rule split_raw_umis:
    input:
        raw_umis = f'{analysis_path}/umi_stats/{{barcode}}/extracted_raw_umis.json',
        stats = f'{analysis_path}/umi_stats/{{barcode}}/reads_stats.json'
    output: f'{analysis_path}/umi_stats/{{barcode}}/cc_umis.json'
    log: f'{analysis_path}/logs/umi_deduplication/{{barcode}}_umi_dedup_log.txt'
    threads: 10
    shell:
        'echo ">>>>> SPLITTING RAW UMIS <<<<<" >> {log} && '
        'python3 -u ../scripts/pipeline/split_raw_umis.py '
        '-r {input.raw_umis} '
        '-o {output} '
        '-s {input.stats} '
        '-v 3 '
        '>> {log}'

# Compute alignment and overlap distances between UMIs
rule compute_pairwise_dist:
    input:
        cc_umis = f'{analysis_path}/umi_stats/{{barcode}}/cc_umis.json',
        stats = f'{analysis_path}/umi_stats/{{barcode}}/reads_stats.json'
    output:
        align_dist_mat = f'{analysis_path}/umi_stats/{{barcode}}/umis_pairwise_dist_mat_alignment.npz',
        overlap_dist_mat = f'{analysis_path}/umi_stats/{{barcode}}/umis_pairwise_dist_mat_overlap.npz'
    log: f'{analysis_path}/logs/umi_deduplication/{{barcode}}_umi_dedup_log.txt'
    threads: nthreads
    shell:
        'echo ">>>>> COMPUTE PAIRWISE DISTANCES <<<<<" >> {log} && '
        'python3 -u ../scripts/pipeline/compute_umis_pairwise_dist.py '
        '-u {input.cc_umis} '
        '-o {analysis_path}/umi_stats/{wildcards.barcode}/umis_pairwise '
        '-s {input.stats} '
        '-c {nthreads} '
        '-v 3 '
        '>> {log}'

# Build raw-connected-UMIs graph
rule make_UMIs_graph:
    input:
        cc_umis = f'{analysis_path}/umi_stats/{{barcode}}/cc_umis.json',
        stats = f'{analysis_path}/umi_stats/{{barcode}}/reads_stats.json',
        align_dist_mat = f'{analysis_path}/umi_stats/{{barcode}}/umis_pairwise_dist_mat_alignment.npz',
        overlap_dist_mat = f'{analysis_path}/umi_stats/{{barcode}}/umis_pairwise_dist_mat_overlap.npz'
    output: f'{analysis_path}/umi_stats/{{barcode}}/clustering_UMIs_graph.pkl'
    log: f'{analysis_path}/logs/umi_deduplication/{{barcode}}_umi_dedup_log.txt'
    threads: 10
    shell:
        'echo ">>>>> MAKE UMIS GRAPH <<<<<" >> {log} && '
        'python3 -u ../scripts/pipeline/cluster_umis.py '
        '-u {input.cc_umis} '
        '-o {analysis_path}/umi_stats/{wildcards.barcode}/clustering '
        '-am {input.align_dist_mat} '
        '-om {input.overlap_dist_mat} '
        '-s {input.stats} '
        '-v 3 '
        '>> {log}'

# Compute graph statistics
rule get_graph_stats:
    input:
        cc_umis = f'{analysis_path}/umi_stats/{{barcode}}/cc_umis.json',
        stats = f'{analysis_path}/umi_stats/{{barcode}}/reads_stats.json',
        graph = f'{analysis_path}/umi_stats/{{barcode}}/clustering_UMIs_graph.pkl'
    output: f'{analysis_path}/umi_stats/{{barcode}}/cluster_stats.tsv'
    log: f'{analysis_path}/logs/umi_deduplication/{{barcode}}_umi_dedup_log.txt'
    threads: 10
    shell:
        'echo ">>>>> COMPUTE GRAPH STATS <<<<<" >> {log} && '
        'python3 -u ../scripts/pipeline/compute_cluster_stats.py '
        '-u {input.cc_umis} '
        '-o {output} '
        '-g {input.graph} '
        '-s {input.stats} '
        '-v 3 '
        '>> {log}'
    
# Plot top UMIs clusters
rule plot_top_clusters:
    input: f'{analysis_path}/umi_stats/{{barcode}}/clustering_UMIs_graph.pkl'
    output: f'{analysis_path}/umi_stats/{{barcode}}/top_clusters.png'
    log: f'{analysis_path}/logs/umi_deduplication/{{barcode}}_umi_dedup_log.txt'
    threads: 10
    shell:
        'echo ">>>>> PLOTTING TOP CLUSTERS <<<<<" >> {log} && '
        'python3 -u ../scripts/pipeline/plot_clusters.py '
        '-d {analysis_path}/umi_stats/{wildcards.barcode}/top_indiv_clusters '
        '-o {output} '
        '-g {input} '
        '-v 3 '
        '>> {log}'

# Split over-merged UMIs connected clusters
rule split_UMIs_clusters:
    input:
        graph = f'{analysis_path}/umi_stats/{{barcode}}/clustering_UMIs_graph.pkl',
        cc_umis = f'{analysis_path}/umi_stats/{{barcode}}/cc_umis.json',
        stats = f'{analysis_path}/umi_stats/{{barcode}}/reads_stats.json',
        fq = f'{data_path}/umi_deduplication/{{barcode}}/{{barcode}}_07_no-umis.fastq'
    output: f'{analysis_path}/umi_stats/{{barcode}}/{umi_method}_split_clustered_umis.json'
    log: f'{analysis_path}/logs/umi_deduplication/{{barcode}}_umi_dedup_log.txt'
    threads: 10
    shell:
        'echo ">>>>> SPLITTING CLUSTERS <<<<<" >> {log} && '
        'python3 -u ../scripts/pipeline/split_clusters.py '
        '-u {input.cc_umis} '
        '-o {analysis_path}/umi_stats/{wildcards.barcode}/{umi_method}_split '
        '-g {input.graph} '
        '-s {input.stats} '
        '-f {input.fq} '
        '-m {umi_method} '
        '-n 10 '
        '-v 3 '
        '>> {log}'

# Color alignments with assigned UMIs / molecules
rule color_BAM:
    input:
        bam = f'{data_path}/umi_deduplication/{{barcode}}/{{barcode}}_09_no-UMIs.bam',
        umis = f'{analysis_path}/umi_stats/{{barcode}}/{umi_method}_split_clustered_umis.json'
    output: f'{data_path}/umi_deduplication/{{barcode}}/{{barcode}}_10_{umi_method}_colored_with_UMIs.bam'
    log: f'{analysis_path}/logs/umi_deduplication/{{barcode}}_umi_dedup_log.txt'
    threads: 10
    shell:
        'echo ">>>>> COLORING BAM RECORDS <<<<<" >> {log} && '
        'python3 -u ../scripts/pipeline/color_bam_with_umis.py '
        '-u {input.umis} '
        '-b {input.bam} '
        '-o {output} '
        '-v 3 '
        '>> {log}'

# Deduplicate reads assigned to the same original molecule
rule deduplicate_reads:
    input:
        umis = f'{analysis_path}/umi_stats/{{barcode}}/{umi_method}_split_clustered_umis.json',
        stats = f'{analysis_path}/umi_stats/{{barcode}}/reads_stats.json',
        fq = f'{data_path}/umi_deduplication/{{barcode}}/{{barcode}}_07_no-umis.fastq'
    output: f'{data_path}/umi_deduplication/{{barcode}}/{{barcode}}_11_{umi_method}_deduplicated.fasta'
    log: f'{analysis_path}/logs/umi_deduplication/{{barcode}}_umi_dedup_log.txt'
    threads: nthreads
    shell:
        'echo ">>>>> DEDUPLICATING READS <<<<<" >> {log} && '
        'python3 -u ../scripts/pipeline/deduplicate_umis.py '
        '-u {input.umis} '
        '-o {analysis_path}/umi_stats/{wildcards.barcode}/deduplicated '
        '-s {input.stats} '
        '-f {input.fq} '
        '-c {nthreads} '
        '-v 3 '
        '>> {log} && '
        'mv {analysis_path}/umi_stats/{wildcards.barcode}/deduplicated_dedup.fasta '
        '{output}'

# Map deduplicated reads to reference
rule mapping_dedup_reads:
    input:
        ref = f'{project_dir}/Data/references/snakeref/{reference}',
        fa = f'{data_path}/umi_deduplication/{{barcode}}/{{barcode}}_11_{umi_method}_deduplicated.fasta'
    output: f'{data_path}/umi_deduplication/{{barcode}}/{{barcode}}_12_{umi_method}_dedup_mapped.sam'
    threads: nthreads
    shell:
        'minimap2 -ax map-ont -t {nthreads} {input.ref} {input.fa} > {output}'

# Sort deduplicated alignment
rule sorting_dedup:
    input: f'{data_path}/umi_deduplication/{{barcode}}/{{barcode}}_12_{umi_method}_dedup_mapped.sam'
    output: f'{data_path}/umi_deduplication/{{barcode}}/{{barcode}}_13_{umi_method}_dedup_sorted.bam'
    threads: nthreads
    shell:
        'samtools sort --threads {nthreads} {input} -o {output} && '
        'samtools index --threads {nthreads} {output}'

def all_targets(wildcards):
    # static targets computed at parse time (config's out_files,
    # minus the UMI-tail blocks — see below)
    targets = list(out_files)

    # dynamic UMI-tail targets, gated by the checkpoint
    umi_trimmed_bc_mapping_ref = umi_survivors(wildcards)
    for b in umi_trimmed_bc_mapping_ref:
        targets += [
            f'{data_path}/umi_deduplication/{b}/{b}_13_{umi_method}_dedup_sorted.bam',
            f'{data_path}/umi_deduplication/{b}/{b}_10_{umi_method}_colored_with_UMIs.bam',
            f'{analysis_path}/umi_stats/{b}/top_clusters.png',
            f'{analysis_path}/umi_stats/{b}/cluster_stats.tsv',
        ]
    return targets

rule all:
    input: all_targets # Use function to terminate without error even if some
    # samples do not have UMI-trimmed reads mapping the reference.
    # For previous behavior, switch to the following:
    # input: expand ('{out_file}', out_file=out_files)
    default_target: True