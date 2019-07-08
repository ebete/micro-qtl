# Snakemake pipeline definition file

# metadata
__author__ = "Thom Griffioen"
__copyright__ = "Copyright 2019 Thom Griffioen"
__email__ = "t.griffioen@nioo.knaw.nl"
__license__ = "MIT"

from snakemake.utils import min_version

min_version("5.0.0")

# configuration
configfile: "config.yml"  # global configuration
# configfile: "samples.yml"  # input files
localrules: all


# project constants
PROJECT = config["project"]
BOOTSTRAPS = list(range(500))

# output files
OUTFILES = list()
# output of GO term prioritisation
OUTFILES.append("{project}/prioritised_terms_{bootstrap_rounds}.txt")

# rule to target all output files
rule all:
    input:
         expand(OUTFILES, project=PROJECT, bootstrap_rounds=BOOTSTRAPS)

rule make_id_map:
    output:
          "{project}/db/id_map.pickle.xz"
    threads: 1
    conda:
         "environment.yml"
    params:
          map_file=config["static-files"]["idmap"],
          from_id=config["id-convert"]["from"],
          to_id="go"
    shell:
         'python3 make_uniprot_idmapping_db.py -vv --from "{params.from_id}" --to "{params.to_id}" "{params.map_file}" "{output}"'

rule make_random_peaks:
    output:
          "{project}/peaks_{bootstrap_rounds}.csv"
    threads: 1
    conda:
         "environment.yml"
    params:
          gen_file=config["static-files"]["genfile"],
          peaks=3,
          width_cm=10.0
    shell:
         'python3 read_genfile.py -i "{params.peaks}" -w "{params.width_cm}" "{params.gen_file}" > "{output}"'

rule genes_from_peaks:
    input:
         id_map=rules.make_id_map.output,
         peaks_file=rules.make_random_peaks.output
    output:
          "{project}/genes_in_peaks_{bootstrap_rounds}.csv"
    threads: 1
    conda:
         "environment.yml"
    params:
          genome_gff=config["static-files"]["refseq"]
          # peaks_file=config["static-files"]["lod-regions"]
    shell:
         'python3 get_go_in_region.py -vv -g "{params.genome_gff}" -m "{input.id_map}" "{input.peaks_file}" > "{output}"'

rule make_go_tree:
    output:
          "{project}/db/go_tree.pickle.xz"
    threads: 1
    conda:
         "environment.yml"
    params:
          obo_file=config["static-files"]["go"]
    shell:
         'python3 create_go_tree.py -vv -o "{output}" "{params.obo_file}"'

rule go_network_analysis:
    input:
         go_tree=rules.make_go_tree.output,
         peak_terms=rules.genes_from_peaks.output,
         peak_regions=rules.make_random_peaks.output
    output:
          "{project}/prioritised_terms_{bootstrap_rounds}.txt"
    threads: 1
    conda:
         "environment.yml"
    shell:
         #         'python3 go_network_analysis.py -vv "{input.go_tree}" "{input.peak_terms}" | dot -Tsvg > "{output}"'
         'python3 go_network_analysis.py -vv "{input.go_tree}" "{input.peak_terms}" "{input.peak_regions}" > "{output}"'
