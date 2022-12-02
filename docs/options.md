
## Options

Required parameters:

| Parameters | Type | Default | Description                   |
|------------|------|---------|-------------------------------|
| --out_dir  | PATH | None    | output folder                 |
| --bam_fn   | FILE | None    | input bam file                |
| --bed_fn   | FILE | None    | input target regions bed file |
| --ref_fn   | FILE | None    | input reference fasta         |

Other parameters:

[EXPERIMENTAL] are for testing and advance usage.

| Type                               | Parameters                    | Type  | Default | Description                                                                                                                                                                                                              |
|------------------------------------|-------------------------------|-------|---------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| General options                    | --threads                     | INT   | 48      | running threads, we recommend using 48 or above                                                                                                                                                                          |
| -                                  | --clair_ensemble_threads      | INT   | 16      | Clair-Ensemble threads, we recommend using 16                                                                                                                                                                            |
| -                                  | --subtype_parallel            | INT   | 3       | [EXPERIMENTAL] number of sutypes parallel run Clair                                                                                                                                                                      |
| initial filtering options          | --indel_l                     | INT   | 50      | filtering read with indel length > indel_l [50], set [0] to disable filtering                                                                                                                                            |
| clusterV options                   | --top_k                       | INT   | 25      | top k sutypes to output                                                                                                                                                                                                  |
| -                                  | --n_min_supports              | INT   | 50      | minimum read support for creating a subtype                                                                                                                                                                              |
| -                                  | --n_max_candidates            | INT   | 15      | [EXPERIMENTAL] number of selected candidates for clustering                                                                                                                                                              |
| -                                  | --min_af                      | FLOAT | 0.05    | [EXPERIMENTAL] minimum AF when cluastering                                                                                                                                                                               |
| -                                  | --n_max_coverage              | INT   | 10000   | [EXPERIMENTAL] max read for clustering                                                                                                                                                                                   |
| consensus and HIVDB report options | --hivdb_url                   | STR   | ""      | hivdb URL default query from the internet, for localizing the HIVDB, please check https://github.com/hivdb/sierra, and update this setting accordingly, e.g. using --hivdb_url http://localhost:8111/sierra/rest/graphql |
| -                                  | --number_of_read_for_consense | INT   | 1000    | [EXPERIMENTAL] number of original read for generating consense                                                                                                                                                           |
| -                                  | --flye_genome_size            | STR   | 5k      | [EXPERIMEANTAL], flye --genome-size for generating consensus, we recommand using 5k for HIV genome                                                                                                                       |
