
## Options

Required parameters:

| Parameters | Type | Default | Description                   |
|------------|------|---------|-------------------------------|
| --out_dir  | PATH | None    | output folder                 |
| --bam_fn   | FILE | None    | input bam file                |
| --bed_fn   | FILE | None    | input target regions bed file |
| --ref_fn   | FILE | None    | input reference fasta         |

Other parameters:

[EXPERIMENTAL] are for testing and advanced usage.

| Type                               | Parameters                    | Type  | Default | Description                                                                                                                                                                                                              |
|------------------------------------|-------------------------------|-------|---------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| General options                    | --threads                     | INT   | 48      | running threads, we recommend using 48 or above                                                                                                                                                                          |
| -                                  | --clair_ensemble_threads      | INT   | 16      | Clair-Ensemble threads, we recommend using 16                                                                                                                                                                            |
| -                                  | --subtype_parallel            | INT   | 3       | [EXPERIMENTAL] number of subtypes parallel run Clair                                                                                                                                                                      |
| initial filtering options          | --indel_l                     | INT   | 50      | filtering read with indel length > indel_l [50], set [0] to disable filtering                                                                                                                                            |
| clusterV options                   | --top_k                       | INT   | 25      | top k subtypes to output                                                                                                                                                                                                  |
| -                                  | --n_min_supports              | INT   | 50      | minimum read support for creating a subtype                                                                                                                                                                              |
| -                                  | --n_max_candidates            | INT   | 15      | [EXPERIMENTAL] number of selected candidates for clustering                                                                                                                                                              |
| -                                  | --min_af                      | FLOAT | 0.05    | [EXPERIMENTAL] minimum AF when clustering                                                                                                                                                                               |
| -                                  | --n_max_coverage              | INT   | 10000   | [EXPERIMENTAL] maximum read for clustering                                                                                                                                                                                   |
| consensus and HIVDB report options | --hivdb_url                   | STR   | ""      | hivdb URL default query from the internet, for localizing the HIVDB, please check https://github.com/hivdb/sierra, and update this setting accordingly, e.g. using --hivdb_url http://localhost:8111/sierra/rest/graphql |
| -                                  | --n_of_read_for_consensus | INT   | 1000    | [EXPERIMENTAL] number of original read for generating consensus                                                                                                                                                           |
| -                                  | --flye_genome_size            | STR   | 5k      | [EXPERIMEANTAL], flye --genome-size for generating consensus, we recommend using 5k for HIV genome                                                                                                                       |
| -                                  | --flye_nano_type              | STR   | nano-hq | [EXPERIMEANTAL], flye option for different ONT type, default --nano-hq, check https://github.com/fenderglass/Flye/blob/flye/docs/USAGE.md |
