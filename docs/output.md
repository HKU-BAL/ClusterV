

## Understand outputs file

### The directory structure for the clusterV results is as below:

    .
    ├── ${OUTPUT_DIR}
    │   ├── consensus                   # output of clusterV
    │   │   ├── all_report.tsv          # drug resistance report for found subtypes
    │   │   ├── all_info.tsv            # statistics for all found subtypes, with abundance information
    │   │   ├── id_1_report.tsv         # HIVDB report for subtype 1
    │   │   ├── id_1_report.json        # [raw HIVDB](https://hivdb.stanford.edu/page/release-notes/#json.output) report output for subtype 1
    │   │   ├── id_1.fasta              # consensus for subtype 1
    │   │   ├── id_1                    # data directory for subtype 1
    │   │   │   ├── id_1.bam            # all read alignment in bam file for subtype 1
    │   │   │   ├── id_1_cs.bam         # single consensus alignment for subtype 1
    │   │   │   ├── id_1_cs.vcf         # all variants in consensus for subtype 1
    │   │   │   ├── id_1_low_af.vcf     # low AF variants in subtypes but not in consensus for subtype 1
    │   │   │   └── ...
    │   │   ├── id_2                    # data directory for subtype 2
    │   │   │   └── ...
    │   │   └── ...
    │   ├── clustering                      # working files for clusterV algorithm
    │   │   ├── all_clusters_info.tsv       # clustering results, subtypes, abundance, and files location information 
    │   │   ├── 1                           # first iteration
    │   │   │   └── 1.tagged.bam  		# bam file with subtype labeling after 1st iteration                 
    │   │   │   └── 1.tagged.tag1.bam  	# subtype 1 after 1st iteration                 
    │   │   │   └── 1.tagged.tag2.bam  	# subtype 2 after 1st iteration                 
    │   │   │   └── ...                 
    │   │   ├── 1_1                     # subtype after first iteration 
    │   │   │   └── ...                 
    │   │   └── ...             
    │   ├── ${SAMPLE_ID}_f.bam          # BAM file after initalize filtering
    │   └── ...                 
    └── ...

### drug resistance report for found subtypes

all drug resistance report is available at `all_report.tsv`

| column          | description                                                                                                                           |
|-----------------|-------------------------------------------------------------------------------------------------------------------------------------|
| gene            | [Gene name](https://hivdb.stanford.edu/page/release-notes/#drm.classification) from HIVDB                                           |
| drugClass       | Drug class name from HIVDB                                                                                                          |
| drugName        | Drug name from HIVDB                                                                                                                |
| drugScore       | [Drug score](https://hivdb.stanford.edu/page/release-notes/#drm.penalty.scores.and.resistance.interpretation) from HIVDB            |
| resistanceLevel | [Drug resistance level](https://hivdb.stanford.edu/page/release-notes/#drm.penalty.scores.and.resistance.interpretation) from HIVDB |
| subtype_ori     | Subtype raw ID when running clustering                                                                                              |
| subtype         | Subtype ID                                                                                                                          |
| abundance       | Subtype's abundance                                                                                                                 |
| is_in_consensus | Whether the drug resistance mutation found in consensus (1 if found in subtype, 0 if found in low AF set)                           |
| VAF             | variants AF in subtype                                                                                                              |
| mutationPos     | mutation genome position in HIV reference                                                                                           |
| mutation        | mutation name                                                                                                                       |
| mutationScore   | [mutation sscore](https://hivdb.stanford.edu/page/release-notes/#drm.penalty.scores.and.resistance.interpretation) in HIVDB         |
| mutationType    | [mutation type](https://hivdb.stanford.edu/page/release-notes/#drm.classification) in HIVDB                                         |
| comments        | [mutation comments](https://hivdb.stanford.edu/page/release-notes/#comments) from HIVDB                                             |



### statistics for all found subtypes, with abundance information

all subtype information is available at `all_info.tsv`

| column                                               | description                                                                                                                                |
|------------------------------------------------------|--------------------------------------------------------------------------------------------------------------------------------------------|
| _idx                                                 | subtype ID                                                                                                                                 |
| _vcf                                                 | raw VCF for subtypes during clustering                                                                                                     |
| _bam                                                 | raw BAM for the subtype during clustering                                                                                                  |
| _out_dir                                             | working directory for the subtype during clustering                                                                                        |
| _sample_id                                           | sample ID                                                                                                                                  |
| _sample_idx                                          | raw index for the subtype                                                                                                                  |
| is_check_clustering                                  | checking status for subtypes, whether try splitting in this subtype                                                                        |
| percentage                                           | subtype's abundance                                                                                                                        |
| _v_info:coverage;snp_c;indel_c;median_af;af1;af2;af3 | subtype's coverage, # of SNP, # of INDEL, the median of AF, # of raw AF in [0, 0.3], # of raw AF in [0.3, 0.7], # of raw AF in [0.7, 1.0]  |

