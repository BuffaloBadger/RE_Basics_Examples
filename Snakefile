include: "reb_3_1/Snakefile"
include: "reb_3_2/Snakefile"
include: "reb_3_3/Snakefile"
include: "reb_3_4/Snakefile"
include: "reb_4_1/Snakefile"
include: "reb_4_2/Snakefile"
include: "reb_4_3/Snakefile"
include: "reb_4_4/Snakefile"
include: "reb_9_1/Snakefile"
include: "reb_9_2/Snakefile"
include: "reb_9_3/Snakefile"
include: "reb_9_4/Snakefile"

# Specify that function  as the input for the all rule
rule all:
    input:
