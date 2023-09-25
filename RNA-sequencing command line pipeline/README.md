# RNA-Sequencing command line data processing pipeline
Included in this directory is a command line (shell script) RNA-sequencing data processing pipeline

The pipeline uses a Conda environment (YAML file) with version-controlled packages called within the script

The pipeline was implemented on Amazon Web Services EC2 to take advantage of faster computing/processing

To run the program on command line:

```sh
conda activate environment_name

./RNA-sequencing_analysis_pipeline.sh project_name organism_name path_to_fastq_directory fasta_reference seq_adapter gtf_ref
```

