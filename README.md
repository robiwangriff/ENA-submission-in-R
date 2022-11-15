# ENA-submission-in-R
Scripts for submission of amplicon sequence data to the ENA using R.

# Why submit to ENA using R?
Submitting sequences to the ENA is an involved process, and it can take time. There are a variety of ways to submit to the [ENA](https://ena-docs.readthedocs.io/en/latest/submit/general-guide.html), but I'm unaware of any routines for doing it in R. Since most who work on amplicon sequence data use R, it makes sense to integrate sequence submission into R workflows.
Personally, in the past I have tended to muddle through the submission process, 50:50 R/fiddling with spreadsheets. Formalising it all in R at least allows me to reproduce my mistakes (and hopefully improve the process). Contribs and feedback most welcome.

# The functions:
This repository contains a set of R functions to submit amplicon sequences, such as paired-end reads from Illumina MiSeq sequencing of taxonomic marker genes amplified from environmental DNA. A file containing the functions is [here](https://github.com/robiwangriff/ENA-submission-in-R/blob/main/ena_subm_fxns.R), and an example workflow is [here](https://github.com/robiwangriff/ENA-submission-in-R/blob/main/ena_subm_example.R). The example scripts can be worked through using a provided [environmental metadata file](https://github.com/robiwangriff/ENA-submission-in-R/blob/main/env_metadata.csv), and some small [exemplar fastq files](https://github.com/robiwangriff/ENA-submission-in-R/tree/main/fastq_files_from_sequencer) (16S rRNA gene and ITS assays performed on the same samples). All being well, these scripts should enable a successful submission of the example datasets to the ENA test server, all from within R. They can then be applied to personal datasets, though with large numbers of fastq files you may want to use a dedicated file transfer program to upload fastqs to ENA, and do go through the example first.

The functions use R to:

1. Transfer "fastqs relevant to a study" from a sequencing output dir (eg the MiSeq dump) to a local submission dir (eg a project/paper dir), based on samples present in an environmental metadata file.
2. Generate md5 checksums for the fastqs
3. Upload selected fastqs to the ENA 
4. Create a sample submission dataframe, based on the user environmental metadata file, and check/select relevant ENA MixS checklists
5. Link samples to fastq reads (multiple amplicon assays possible) by constructing ENA EXPERIMENT and RUN files
6. Convert dataframes to xml format required by the ENA for programatic submission
7. Submit xmls to the ENA 
8. Populate the environmental metadata file with ENA accession numbers, to enable reproducible coding of subsequent analyses.

# Usage
You will need an [ENA account](https://www.ebi.ac.uk/ena/submit/sra/#home) to submit sequences. 

Clone this repository if you are a github user, or click the green button to download the [zipped directory](https://github.com/robiwangriff/ENA-submission-in-R/archive/refs/heads/main.zip). You may need to modify paths in the example file depending on where you put the extracted dir. Then work through the [ena_subm_example.R](https://github.com/robiwangriff/ENA-submission-in-R/blob/main/ena_subm_example.R) file, which loads the functions contained in [ena_subm_fxns.R](https://github.com/robiwangriff/ENA-submission-in-R/blob/main/ena_subm_fxns.R) via a source() command. 

You should get a pretty speedy submission to the ENA dev server using the example files, but more and larger fastq files will of course take longer. Copying large files even locally takes time, as does the md5 generation - so be prepared to forfeit your R session for a bit when submitting "real" datasets. Some progress bars would be useful. These scripts have been used on submissions of ~2000 read files all from within R, but sometimes we have had failed R uploads of fastq files to the ENA, requiring a restart of the file transfer. Since there is nothing to really to be gained from using R for transfer - using a dedicated file transfer program (eg winscp) on your now neatly stored fastqs, may allow you to do something else useful in R whilst files are uploading.

With regard to submitting real datasets, since the user decides what metadata to upload, there is some manual work involved in selecting relevant metadata fields and populating the sample info file according to the MixS checklists. This is intended and necessary to ensure accurate and meaningful submissions.

After xml submission a receipt is returned indicating whether the submission is succesful. You may get failed test submissions because the project/samples/runs have already been registered (from your previous attempt)...even with the dev server. The dev server is supposed to be purged every day, and there doesn't appear to be a way of purging it yourself for testing purposes. On the plus side, this may mean your submission attempt to the dev server has kind of worked.

A final note: be careful with real submissions to the working submission server (dev is fine) as erroneous ones are hard to come back from, without renaming sample aliases and even your fastq files - hardly reproducible, and a pain. Always test on the dev server before submitting your samples. Dont try and submit anything to non-dev that may have been submitted before (eg the test files included here!).

