from os.path import join as j
import numpy as np

configfile: "workflow/config.yaml"

#
# Paper
#
PAPER_DIR = "paper/current/"
PAPER_SRC = j(PAPER_DIR, "main.tex")
PAPER = j(PAPER_DIR, "main.pdf")

#
# Figures
#
FIG_DIR = config["fig_dir"]
FIGS = [j(FIG_DIR, f) for f in ("detected-cartel-stat.pdf", )]

#
# Data Dir
#
MAG_DATA_DIR = config["data_dir"]
MAG_SRC_DATA_DIR = j(MAG_DATA_DIR, "source")
MAG_CLEANED_DATA_DIR = j(MAG_DATA_DIR, "cleaned")

#
# Data base
#
MAG_DB_DIR = j(MAG_DATA_DIR, "database")
DB_CONF_DIR = "workflow"
DBNAME = "magdb"


#
# Download
#
MAG_CONTAINER_KEY = config["container_key"]

MAG_DATA_FILENAME = [
"Journals.txt",
"Authors.txt",
"PaperAuthorAffiliations.txt",
"PaperReferences.txt",
"Papers.txt",
"ConferenceSeries.txt",
]
MAG_DATA_FILE = j(MAG_SRC_DATA_DIR,"{mag_filename}")
MAG_DATA_FILE_ALL = expand(MAG_DATA_FILE, mag_filename = MAG_DATA_FILENAME)

#
# Cleaning
#
MAG_CLEANED_DATA_FILENAME = MAG_DATA_FILENAME + ["PaperJournalAffiliations.txt"]
MAG_CLEANED_DATA_FILE = j(MAG_CLEANED_DATA_DIR,"{mag_filename}")
MAG_CLEANED_DATA_FILE_ALL = expand(MAG_CLEANED_DATA_FILE, mag_filename = MAG_CLEANED_DATA_FILENAME)


# Paper count
NETWORK_DIR = j(MAG_DATA_DIR, "networks")
PAPER_COUNT_FILE = j(NETWORK_DIR, "paper_count.csv")
YEARLY_NODE_FILE = j(NETWORK_DIR, "nodes-{year}.csv")
YEARLY_EDGE_FILE = j(NETWORK_DIR, "edges-{year}.csv")
RAW_YEARLY_NODE_FILE = j(NETWORK_DIR, "raw-nodes-{year}.csv")
RAW_YEARLY_EDGE_FILE = j(NETWORK_DIR, "raw-edges-{year}.csv")

WINDOW_LENGTH = 2
YEARS = list(range(1998, 2010))
YEARLY_NODE_FILE_ALL = expand(YEARLY_NODE_FILE, year = YEARS)
YEARLY_EDGE_FILE_ALL = expand(YEARLY_EDGE_FILE, year = YEARS)
RAW_YEARLY_NODE_FILE_ALL = expand(RAW_YEARLY_NODE_FILE, year = YEARS)
RAW_YEARLY_EDGE_FILE_ALL = expand(RAW_YEARLY_EDGE_FILE, year = YEARS)


# Community detection
COMMUNITY_DIR = j(MAG_DATA_DIR, "community")
AGGREGATED_YEARS = list(range(2000, 2020))
DETECTED_COMMUNITY_FILE = j(COMMUNITY_DIR, "aggregated-community.csv")


# Cartel detection
CARTEL_DIR = j(MAG_DATA_DIR, "cartels")
CARTEL_YEARS = list(range(2000, 2020))
THETA_CIDRE = 0.15
ALPHA_CIDRE = 0.01
DETECTED_CARTEL_FILE = j(CARTEL_DIR, "cartels-{year}.csv")
DETECTED_CARTEL_FILE_ALL = expand(DETECTED_CARTEL_FILE, year=CARTEL_YEARS)

# Plot
FIG_DIR = config["fig_dir"]
FIG_DETECTED_CAETEL_STATS = j(FIG_DIR, "detected-cartel-stat.pdf")
FIG_CITATION_NET_CAETEL = j(FIG_DIR, "citation-net-cartels.pdf")

rule all:
    input: PAPER

rule paper:
    input: PAPER_SRC, FIGS
    params:
        paper_dir=PAPER_DIR
    output: PAPER
    shell: "cd {params.paper_dir}; make"

rule import_neo4j:
    input: MAG_CLEANED_DATA_FILE_ALL
    output: directory(MAG_DB_DIR)
    run:
        shell("bash workflow/mag2neo4j.sh {MAG_DB_DIR} {MAG_CLEANED_DATA_DIR} {DB_CONF_DIR} {DBNAME}")

rule cleanup_mag:
    input: MAG_DATA_FILE_ALL
    output: MAG_CLEANED_DATA_FILE_ALL
    run:
        shell("bash workflow/cleanup_mag_file.sh {MAG_SRC_DATA_DIR} {MAG_CLEANED_DATA_DIR}")

rule download_mag:
    output: MAG_DATA_FILE
    params:
        filename = lambda wildcards: wildcards.mag_filename
    run:
        shell("python3 workflow/get_mag_data.py {params.filename} {MAG_SRC_DATA_DIR} '{MAG_CONTAINER_KEY}'")

rule count_papers:
    input: directory(MAG_DB_DIR)
    output: PAPER_COUNT_FILE
    run:
        shell("python3 workflow/count_papers.py {output}")

rule construct_yearly_networks:
    input: PAPER_COUNT_FILE
    output:
        node = YEARLY_NODE_FILE,
        edge = YEARLY_EDGE_FILE
    params:
        year = lambda wildcards: wildcards.year
    run:
        shell("python3 workflow/construct_yearly_networks.py {input} {NETWORK_DIR} {params.year} {WINDOW_LENGTH} {output.node} {output.edge}")

rule construct_yearly_raw_networks:
    input: PAPER_COUNT_FILE
    output:
        node = RAW_YEARLY_NODE_FILE,
        edge = RAW_YEARLY_EDGE_FILE
    params:
        year = lambda wildcards: wildcards.year
    run:
        shell("python3 workflow/construct_yearly_networks.py {input} {NETWORK_DIR} {params.year} 9999 {output.node} {output.edge}")


rule detect_communities:
    input: YEARLY_NODE_FILE_ALL, YEARLY_EDGE_FILE_ALL, RAW_YEARLY_NODE_FILE_ALL, RAW_YEARLY_EDGE_FILE_ALL
    output: DETECTED_COMMUNITY_FILE
    params:
        years = " ".join(["%d" %d for d in AGGREGATED_YEARS])
    run:
        shell("python3 workflow/community_detection.py {params.years} {output}")

rule detect_cartels:
    input: DETECTED_CARTEL_FILE_ALL

rule detect_cartels_yearly:
    input: YEARLY_NODE_FILE, YEARLY_EDGE_FILE, RAW_YEARLY_NODE_FILE, RAW_YEARLY_EDGE_FILE, DETECTED_COMMUNITY_FILE
    output: DETECTED_CARTEL_FILE
    params:
        year = lambda wildcards : wildcards.year
    run:
        shell("python3 workflow/detect_cartels.py {params.year} {NETWORK_DIR} {THETA_CIDRE} {ALPHA_CIDRE} {DETECTED_COMMUNITY_FILE} {output}")

rule plot_cartel_stats:
    input: CARTEL_DIR
    output: FIG_DETECTED_CAETEL_STATS
    run:
        shell("python3 workflow/plot-cartel-stat.py {CARTEL_DIR} {output}")

