##Worked example for submitting sequences to the ENA


#The following worked example is for processing two amplicons from the same sample (eg 16S & ITS). 
#Will be simpler for one amplicon.

####################################################################################
#1. Load and define some variables
####################################################################################
	
	#########################################################
	#Define ENA webin login info
	#get account if needs https://www.ebi.ac.uk/ena/submit/sra/#registration
	
	#for security, enter YOUR details directly in terminal

	ena_user<-"Webin-xxxxx"
	ena_passwd<-"pxxssxxxd"

	#######################################################
	#set working directory (tut tut)
	#
	#If downloaded a zip file
	#setwd("/where/you/put/extractedzip/ENA-submission-in-R-main/")

	#If cloned github
	#setwd("/loca/github/dir/ENA-submission-in-R/") #note no "-main"
	
	#######################################################
	#load packages #If needs: install.packages("missing_package")
	
	library(fs)
	library(parallel)
	library(digest)
	library(RCurl)
	library(XML)
	library(httr)
	library(jsonlite)

	#####################################################
	#Load env metadata
	#Load environmental and sequence associated metadata (we call it an "env" file)
	#We'll update this as we go along and populate with various submission related fields
	#After doing this a few times it will become apparent that it is good practice
	#to pre-populate the env file with various fields which will facilitate ENA submission.
	
	env<-read.csv("env_metadata.csv",stringsAsFactors=F,check.names=F)

	####################################################
	#Load functions

	source("ena_subm_fxns.R") 
	
	#######################################################
	#Define directories which contain reads (eg path to your sequencer output dir).
	#May well contain more fastq files than you actually want to submit.
	#For worked example this file is provided:	

	seq_dir_16S<-"fastq_files_from_sequencer/16S_runs/"
	seq_dir_ITS<-"fastq_files_from_sequencer/ITS_runs/"

	########################################################
	# Create a local ENA submission directory where you can transfer reads to be uploaded
	# and any other outputs from the submission process

	dir.create("ENA_submission/")# here we are creating in our working dir. You may want to specify a more relevant path eg a project folder
	
	#Define folder paths where you'll transfer only the fastq files you want to submit - 
	#This directory will be created by the function - it just needs a path here

	upload_dir_16S<-"ENA_submission/16S_upload_fastqs/"
	upload_dir_ITS<-"ENA_submission/ITS_upload_fastqs/"


	######################################################
	#Define broad details about project
	#In the example data these are in conveniently found in the env file (or at least some fields which fit the bill), 
	#but you could also define them here: eg proj_name<-"Microbial communities of habitat X" 

	proj_name<-env$"Project name"[1]
	proj_title<-env$"sample_title"[1]
	proj_desc<-env$"sample_description"[1]

	###################################################
	#Define other key variables here?

	#The most important field which links everything is the sample_alias, ie the unique sample ID.
	#In the example env file this is (somewhat clumsily):

	sample_alias=env$"Sample ID for sequencing"

	#The other important columns are the fastq.gz filenames associated with each sample.
	#If your sequencer appends additional fields (ours does) and you dont actually know
	#the fastq filenames, the functions below (2a) should sort this. But you still need the
	#sample IDs as submitted to the sequencer. In our lab setup, these are not only the sample
	#alias, but also include details of location in 96 well PCR plates.

	#In the example env file, since we are dealing with both 16S and ITS amplicon assays
	#which were carried out on different runs, there are two columns
	#env$"Sequence ID 16S" and env$"Sequence ID ITS". Make sure such information is present in your env file
	#These not defined explicitly here but will be defined in the function in 2.

####################################################################################################################
#2. Transfer fastq files from seq output directory to local upload directory
####################################################################################################################
	
	###########################################################################
	#a. Retrieve fastq filenames which are to be submitted, based on sequence IDs in env file, and 
	# fastq filenames in sequencer output directory
	#(optional step: you may have full fastq names in env already - our sequencing protocol appends additional text 
	#separated by an underscore. So need to match the sequence ID in env file with output fastq filename)

	#function: seq_id_from_env
	
	f_r_reads_16S<-get_f_r_names(seq_id_from_env=env$"Sequence ID 16S",seq_dir=seq_dir_16S)
	f_r_reads_ITS<-get_f_r_names(seq_id_from_env=env$"Sequence ID ITS",seq_dir=seq_dir_ITS)

	#Recommend adding outputs to env file throughout- will create a full metadata record including ENA submission details
	#and routine saving (version tracked) will allow you to revisit the process should something go bad.

	env$"16Sf"<-f_r_reads_16S$f#add to env
	env$"16Sr"<-f_r_reads_16S$r

	env$"ITSf"<-f_r_reads_ITS$f
	env$"ITSr"<-f_r_reads_ITS$r


	######################################################################
	#b. Transfer selected files to local ENA upload directory

	#function: copy_seq_files
 	
	copy_seq_files(f_names=f_r_reads_16S[,1],r_names=f_r_reads_16S[,2],seq_dir=seq_dir_16S,ena_upload_dir=upload_dir_16S)
	copy_seq_files(f_names=f_r_reads_ITS[,1],r_names=f_r_reads_ITS[,2],seq_dir=seq_dir_ITS,ena_upload_dir=upload_dir_ITS)


	#These files can now be uploaded to the ena, using an external program (eg winscp) prior to running next
	#steps. If really want to do file transfer in R (not recomended) carry on to (3.) below. Otherwise go to (4.) 
	#whilst your sequences upload via winscp or somesuch.


#########################################################################################################
#3. Upload fastqs to ENA...in R
#########################################################################################################

#Note the below should work for transfer of low numbers of files. With over 400 files (eg 200 samples x 2 with F and R reads)
#the following fails for some reason and I havent investigated upload continuation,should it fail. Additionally, the R session
#will be out of use whilst it transfers files (which will be slow). For these reasons whilst the code is provided here, I 
#recomend using a different file transfer approach (eg winscp), and using your R session for generating the md5 codes (next step)


	r_upload(f_r_reads_16S[,1],f_r_reads_16S[,2],upload_dir_16S,ena_user,ena_passwd)
	r_upload(f_r_reads_ITS[,1],f_r_reads_ITS[,2],upload_dir_ITS,ena_user,ena_passwd)


#########################################################################################################
#4. Create md5 checksum codes for fastq files
#########################################################################################################

#An md5 checksum code is required in the RUN/EXP submission step for each fastq read to check file integrity after transfer
#The following generates them within R, using the f and r fastq name field to reference.
#You should then add these to your env file for completeness.
#Note md5 generating using the function below is slow (can be a few seconds per file, depending on size).
#Even with parallelisation in the function, expect this to hog your R session for a bit. Could do with a progress bar.

	md5s_16S<-md5_gen(f_names=f_r_reads_16S[,1],r_names=f_r_reads_16S[,2],upload_dir_16S)
	md5s_ITS<-md5_gen(f_names=f_r_reads_ITS[,1],r_names=f_r_reads_ITS[,2],upload_dir_ITS)

	env$md5_16Sf<-md5s_16S[,1]#add to env
	env$md5_16Sr<-md5s_16S[,2]

	env$md5_ITSf<-md5s_ITS[,1]
	env$md5_ITSr<-md5s_ITS[,2]

#########################################################################################################
#5. Create sample metadata files
#########################################################################################################

# There is no simple function to automate this, since the sample metadata file is populated according
# to user preference(though incorporating mandatory fields from the relevant ENA checklist). The following are 
# helper R functions/code designed to assist in choosing and using the ENA checklists, and then populating your sample
# metadata file (samp_info) according to mixS standards.

	############################################################################
	#a. Select a checklist appropriate for your samples
	require(jsonlite)
	checklists<-fromJSON("https://submission.ebi.ac.uk/api/checklists")
	data.frame(cl.id=checklists$"_embedded"$checklists$id,cl.desc=checklists$"_embedded"$checklists$displayName)

	select.cl<-"ERC000011" #not listed above (?) but the default checklist with no mandatory criteria (not recomended)
	select.cl<-"ERC000022" # a checklist for soil samples

	##################################################################
	#b. Inspect selected checklist(select.cl above) as an R data.frame. 
	#Includes fields, units, vocabularies, and whether mandatory or optional. 
	#Requires xml to df conversion...only tested for ERC000022 
	
	cl_022<-cl_as_df(select.cl="ERC000022")

	##################################################################
	#c. Build sample information dataframe (samp_info)
	#The following function will build a samp_info dataframe based on the user environmental
	#metadata file (env) according to the specified checklist. Its a biggie. 
	#Note: You dont need to call the output "samp_info", 
	#and your input metadata file doesnt need to be called "env"...but it will possibly help if they are (you'll see).

	#Alongside the checklist fields, you will need to identify a sample_alias - ie a unique name for each sample
	#and a "taxonomic id" (ENA lingo) and associated "scientific name" (here I default to tax_id= 256318; scientific_name="metagenome").
	#More on these IDs here:
	#https://ena-docs.readthedocs.io/en/latest/retrieval/programmatic-access/taxon-api.html
	#and you can look up from R with eg
	
	library(httr)
	content(GET("https://www.ebi.ac.uk/ena/taxonomy/rest/suggest-for-submission/soil*"))
	content(GET("https://www.ebi.ac.uk/ena/taxonomy/rest/scientific-name/metagenome"))

	#The function also allows you to specify whether to only include mandatory fields (all_fields=TRUE is the default),
	#and conveniently prints out some commands to help you populate your output samp_info dataframe. If a checklist FIELD
	#requires units to be specified, columns will be created using the the field name suffixed with "_UNITS".
	#You can also add any other fields from your env file when running the function, even if they not in the checklist.
	#Be careful though with (totally normal) missing data in checklist fields. "NA" (or other null designations) in any sample fields
	#could cause a submission failure. If in doubt, dont submit it, as your env file which goes with your paper/
	#to some other repo will have all the relevant metadata...even if one day the machines cant find it.

	samp_info<-build_sample_df(sample_alias=env$"Sample ID for sequencing",select.cl="ERC000022",all_fields=FALSE)

	#The function will print something like this, allowing you to copy it into your editor, and populate, then run in terminal.
	#more hints below

samp_info$"project name"<-env$"Project name"              #mandatory
samp_info$"sequencing method"<-"Illumina MiSeq"               #mandatory
samp_info$"investigation type"<-"mimarks-survey"             #mandatory
samp_info$"collection date"<-as.character(as.Date(as.character(env$sampling_date), format = "%d/%m/%Y") )             #mandatory
samp_info$"geographic location (country and/or sea)"<-"Austria"              #mandatory
samp_info$"geographic location (latitude)"<-env$lat               #mandatory
samp_info$"geographic location (longitude)"<-env$long               #mandatory
samp_info$"soil environmental package"<-"soil"            #mandatory
samp_info$"geographic location (depth)"<-"0.1"              #mandatory   #not in env :manually added
samp_info$"environment (biome)"<-"terrestrial biome"              #mandatory
samp_info$"environment (feature)"<-"soil"              #mandatory
samp_info$"environment (material)"<-"soil"             #mandatory
samp_info$"geographic location (elevation)"<-env$elev              #mandatory

#samp_info$pH<-env$Soil_pH # I tried to add this one but caused a failure on full dataset due to NA

samp_info$"UNIT_geographic location (latitude)"<- "DD"           #unit options:DD
samp_info$"UNIT_geographic location (longitude)"<- "DD"           #unit options:DD
samp_info$"UNIT_geographic location (depth)"<-     "m"       #unit options:m
samp_info$"UNIT_geographic location (elevation)"<- "m"           #unit options:m


	
	#Hints on the above: 
	#Select relevant fields and populate by reference to names(env) or enter appropriate text directly.
	#Google relevant env ontologies and envo terms. For soil, biome= terrestrial biome, feature=soil,material=soil *I think*
	#Date needs to be yyyy-mm-dd....or I think yyyy works. See trick above if date in wrong format.
	#Apparently use: "not collected","not provided", or "restricted access" if data not available..but not sure these work for numeric fields.


	#Note: You could save this samp_info file and submit via web webin uploader if dont want to use the R route
	#but will involve some non reproducible spreadsheet fiddling.
	#this may help:
	#select.cl<-"ERC000022"
	#comment<-paste("#checklist_accession",select.cl, sep="\t")
	#writeLines(comment,"samp_info.txt")
	#write.table(samp_info,"samp_info.txt",append=T,sep = "\t",row.names = FALSE,quote=FALSE)


	#################################################################################################################
	#d. Convert samp_info to xml file for R submission to ENA

	#function: create_sample_xml

	samples_xml<-create_sample_xml(samp_info=samp_info,center_name="UKCEH",ERC_checklist="ERC000022")
	saveXML(samples_xml, file="ENA_submission/samples.xml")
	

#################################################################################################################
#6. Prepare reads (RUN and EXPERIMENT) files for submission to ENA
#################################################################################################################
	
#The definition of RUNS and EXPERIMENT files is slightly ambiguous, but two xml files are needed nevertheless:
#an "EXPERIMENT" and a "RUN" xml. If you have multiple amplicons per sample they can be distinguished here by creating 
#unique EXPERIMENT and RUN aliases, with information on the amplicon used.
#These functions are a way to get it to work, it may not be what is actually required and meant
# within the ENA definitions of RUNS and EXPERIMENTS...

# A "reads_info" df will be created which references the sample alias to fastq reads, potentially for
#multiple amplicons. The RUN and EXP xmls are then created from this. Additional mandatory info is required
#Note the reads_info file can be uploaded via webin uploader, but goal here is for direct R submission. 

	##########################################################
	#a. Create reads_info dataframe

	#function: make_reads_info
	#It is barely a function, but conveniently fills in a few fields. Check them!

	reads_info_16S<-make_reads_info(library_name="16S",sample_alias=env$"Sample ID for sequencing",
	f_names=env$"16Sf",r_names=env$"16Sr",md5f=env$"md5_16Sf",md5r=env$"md5_16Sr")

	reads_info_ITS<-make_reads_info(library_name="ITS",sample_alias=env$"Sample ID for sequencing",
	f_names=env$"ITSf",r_names=env$"ITSr",md5f=env$"md5_ITSf",md5r=env$"md5_ITSr")

	#No need for "_16S" if only one amplicon

	#Important: Need to add exp_alias and run_alias fields now created in reads_info file
	# to the env file. This will later allow matching generated EBI accessions to samples 

	env$exp_alias_16S<-reads_info_16S$exp_alias
	env$run_alias_16S<-reads_info_16S$run_alias
	env$exp_alias_ITS<-reads_info_ITS$exp_alias
	env$run_alias_ITS<-reads_info_ITS$run_alias

	#Now rowbind the 16S and ITS files together into one reads_info file. No need for following if single amplicon assay.
	
	reads_info<-rbind(reads_info_16S,reads_info_ITS)

	#############################################################################
	#b. Create Experiment.xml from reads_info
	#function: create_exp_xml
	
	exp_xml<-create_exp_xml(reads_info=reads_info,proj_name=proj_name,title=proj_title)
	saveXML(exp_xml, file="ENA_submission/experiment.xml")


	#############################################################################
	#c. Create RUN.xml from reads_info
	#function: create_run_xml

	run_xml<-create_run_xml(reads_info=reads_info,center_name="UKCEH")
	saveXML(run_xml, file="ENA_submission/run.xml")


#########################################################################################################
#7. Submit xml files to ENA
#########################################################################################################

# All following scripts use the test server. This is best to test everything is working.
# remove "dev" from the url when ready to submit. Be aware some corrections are possible once submitted...
# but more xml...best get it right first time.


	#####################################################################################################
	# a. Register project

	##############################################################
	#Create Project xml with basic title and description of project

	proj_xml<-create_proj_xml(proj_name,proj_title,proj_desc)
	saveXML(proj_xml, file="ENA_submission/project.xml")

	##################################################################################################
	#Create submission xml for the project registration, to tell ENA to action the project submission

	# specify release date
	rel_date<-"2021-06-06"

	subm_proj_xml<-create_proj_subm_xml(rel_date)
	saveXML(subm_proj_xml, file="ENA_submission/project_submit.xml")


	######################################################################################
	#send project xml and submission xml to ena
	#The following pushes the xml to the ena, and a receipt is returned notifying of success, or not.
	#

	require(httr)
	#ena_user="******"
	#ena_passwd="******"
	project_xml="ENA_submission/project.xml"
	project_submit_xml="ENA_submission/project_submit.xml"
	url<-"https://wwwdev.ebi.ac.uk/ena/submit/drop-box/submit/"

	proj_receipt_xml <- POST(url,
	body=list( SUBMISSION=upload_file(project_submit_xml),
	PROJECT=upload_file( project_xml)),
	authenticate(ena_user,ena_passwd))

	status_code(proj_receipt_xml)
	http_status(proj_receipt_xml)
	headers(proj_receipt_xml)
	content(proj_receipt_xml)
	xmlParse(proj_receipt_xml) # check for errors
	
	saveXML(xmlParse(proj_receipt_xml), file="ENA_submission/RECEIPT_project.xml")



	######################################################################################
	##b. Register samples
	
	#send SAMPLE xmls to ena
	#make a submission xml to action the submission. This can also be used for actioning RUN and EXP submissions
	
	subm_SRE<-create_SRE_subm_xml()
	saveXML(subm_SRE, file="ENA_submission/SRE_submit.xml")

	#send sample xml and submission xml to ena

	samples_xml="ENA_submission/samples.xml"
	subm_SRE_xml="ENA_submission/SRE_submit.xml"

	samp_receipt_xml <- POST(url,
	body=list( SUBMISSION=upload_file(subm_SRE_xml),
	SAMPLE=upload_file(samples_xml)),
	authenticate(ena_user,ena_passwd))

	xmlParse(samp_receipt_xml ) # check for errors
      saveXML(xmlParse(samp_receipt_xml), file="ENA_submission/RECEIPT_samples.xml")


	######################################################################################
	##c. Register EXP and RUN xmls 

	run_xml="ENA_submission/run.xml"
	exp_xml="ENA_submission/experiment.xml"
	subm_SRE_xml="ENA_submission/SRE_submit.xml"

	RE_receipt_xml <- POST(url,
	body=list( SUBMISSION=upload_file( subm_SRE_xml),
	EXPERIMENT=upload_file( exp_xml),
	RUN=upload_file( run_xml)),
	authenticate(ena_user,ena_passwd))

	xmlParse(RE_receipt_xml) # check for errors
	saveXML(xmlParse(RE_receipt_xml), file="ENA_submission/RECEIPT_run_experiment.xml")


########Submit your sequences to ENA
#If all well, repeat 7 with 
#url<-"https://www.ebi.ac.uk/ena/submit/drop-box/submit/"
#NOT DONE - please dont actually submit these example files


#################
##DONE submission

#########################################################################################################
#8. Add accession details to env file
#########################################################################################################

#Now the sequence reads are submitted, we need to add the accession information into the metadata env file.
#This file can then be submitted to a metadata repository of your choosing (or article Supp info) to enable 
#any subsequent analyses to be fully reproducible.

#This can be done by parsing the xml receipts obtained above to dataframes and then adding accessions to env with simple lookups.

####add accessions to env file and export

#parse xml RUN and EXP ("RE") receipt to a dataframe to obtain run and experiment accession numbers

RE_df<-RE_receipt_to_df(RE_receipt_xml)

# Add accessions to env file by lookup against the run_alias and exp_alias created in 6a

env$run_acc_16S<-RE_df[match(env$run_alias_16S,RE_df$run_alias),]$run_acc 
env$exp_acc_16S<-RE_df[match(env$exp_alias_16S,RE_df$exp_alias),]$exp_acc 

env$run_acc_ITS<-RE_df[match(env$run_alias_ITS,RE_df$run_alias),]$run_acc
env$exp_acc_ITS<-RE_df[match(env$exp_alias_ITS,RE_df$exp_alias),]$exp_acc 


#write.csv(env,"ENA_submission/env_with_seq_accessions.csv")

#########################################################################################################
#9. Do some reproducible analyses
#########################################################################################################

#Get fastq files from ENA https://www.ebi.ac.uk/ena/browser/view/PRJN.....

#env<-read.csv(ENA_submission/env_with_seq_accessions.csv)

library(DADA2)
library(vegan)
...






