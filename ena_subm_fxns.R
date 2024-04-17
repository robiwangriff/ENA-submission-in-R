#################################################################
######ENA Upload Functions

#############################################################################
#Helper function: get_f_r_names

#' We need two columns in the metadata env file containing the f and r filenames.
#' This step would ideally not be required if they were included in the original env metadata
#' file - this should be possible at the "lab stage" before samples are submitted to the sequencer.
#' IF you have them skip this step, if not the following
#' works for how our sequencer is set up to name the output files.
#'
#' In our case we give the sequencer a "-" separated unique ID for a sample:
#' eg ALPINE-3-C-F-04-Sol-8-3-16S-1 (should be in env file)
#' 
#' then the sequencer outputs:
#'
#' ALPINE-3-C-F-04-Sol-8-3-16S-1_S240_L001_R2_001.fastq.gz
#' 
#' ...so we need to get a list of all the fastqs orginally sequenced, and match these to the sequence 
#' names in the env file. Like this:

# The following function takes a vector of sequence names (eg seq_id_16S or seq_id_ITS in envfile) 
# and outputs a two column dataframe of f and r read names. You may want to append these to your envfile

get_f_r_names<-function(seq_id_from_env,seq_dir){

#get sequence names of all files in the original sequence file directory
seqnames<-data.frame(nm=list.files(seq_dir))

#recreate the orginal sample name provided by the user. This could be specific to our usage situation
seqnames$ID<-sub('_.*', '', seqnames$nm)

#add columns with R1 and R2 file names (EBI wants separate columns denoting these)
seqnames$name_read1<-ifelse(grepl("_R1_",seqnames$nm),paste(seqnames$nm),"")
seqnames$name_read2<-ifelse(grepl("_R2_",seqnames$nm),paste(seqnames$nm),"")

#now remove duplicates (due to having f & r reads)
tmp1<-seqnames[seqnames$name_read1!="",]
tmp2<-seqnames[seqnames$name_read2!="",]

#output two columns containing names of f & r reads
out<-tmp1
out$name_read2<-tmp2$name_read2#adds read 2 names back in
out$nm<-NULL

#now subset the output to only those seq_ids in the env file
f<-out[match(seq_id_from_env,out$ID),]$name_read1
r<-out[match(seq_id_from_env,out$ID),]$name_read2

return(data.frame(f,r,stringsAsFactors=FALSE))
}


#example usage:
#f_r_reads<-get_f_r_names(seq_id_from_env=env$"Sequence ID 16S",seq_dir=seq_dir)
#you might want to add to your env file to generate a final env with all sequence submission associated inf
#env$"16Sf"<-f_r_reads$f
#env$"16Sr"<-f_r_reads$r


####################################################################
#function: copy_seq_files

#copy the selected fastq files (ie those files associated with samples in env file) from your sequencer output directory
#to your local ENA upload directory

copy_seq_files<-function(f_names,r_names,seq_dir,ena_upload_dir){
myfiles<-c(c(f_names,r_names))#create vector of the required f and r fastqs

#remove NAs and output notification

if (any(is.na(myfiles))) {print("Alert:The dataframe of sequence names contains NAs. Check this is intended.")}
myfiles<-na.omit(myfiles)

myfiles1<-paste(seq_dir,myfiles,sep="")#add in path to directory where these files are
dir.create(file.path(ena_upload_dir), showWarnings = FALSE)
require(fs)
file_copy(myfiles1, ena_upload_dir)#copy selected files to ENA upload directory
}


#example usage: be careful you are about to copy many files. Can take time, will fail if sequence files already exist in target dir.
#copy_seq_files(f_names=f_r_reads[,1],r_names=f_r_reads[,2],seq_dir=seq_dir,ena_upload_dir=ena_upload_dir)


#####################################################################
#function: md5_gen

# MD5 checksum codes are required by the ENA to check files have transferred correctly. 
# Here we will create a 32 character string MD5 unigue to each fastq file, so ENA can validate after upload 
#
# Create function using the forward and reverse sequence names to generate a two column data.frame
# containing the md5 codes. Generating md5 codes is slow (several seconds per file) so this is parallelised.


md5_gen<-function(f_names,r_names,ena_upload_dir){
require(parallel)
require(digest)

n.cores <- detectCores()
cl <- makeCluster(n.cores)
md5f<-parSapply(cl,paste0(ena_upload_dir,f_names), digest,file=TRUE, algo="md5")
stopCluster(cl)

cl <- makeCluster(n.cores)
md5r<-parSapply(cl,paste0(ena_upload_dir,r_names), digest,file=TRUE, algo="md5")
stopCluster(cl)

return(data.frame(md5f,md5r,stringsAsFactors=FALSE))}


#example use
#md5s<-md5_gen(f_names=f_r_reads[,1],r_names=f_r_reads[,2],ena_upload_dir)
#you may also want to add to your env
#env$md5_16Sf<-md5s[,1]
#env$md5_16Sr<-md5s[,2]


###############################################################################
#function: r_upload

#function to upload the fastq files from local ENA submission directory
#to ENA upload area associated with the user webin account.

# Important note: uploading from R works, but is slow and can fail after a few hundred files have been uploaded.
# This could perhaps be resolved with a better function, but its simple enough to use a non R 
# alternative such as winscp.


r_upload<-function(f_names,r_names,ena_upload_dir,ena_user,ena_passwd){
require(RCurl)

#make a vector containing all the files you want to upload
files2upload_A<-c(f_names,r_names)
files2upload_B<-paste0(ena_upload_dir,files2upload_A)

files2upload=data.frame(files2upload_A, files2upload_B, stringsAsFactors = F)

for(i in 1:nrow(files2upload)){
ftpUpload(files2upload[i,2],
paste0("ftp://",ena_user,":",ena_passwd,"@webin.ebi.ac.uk/",files2upload[i,1]))
  }

}

#example use
#r_upload(f_r_reads[,1],r_names=f_r_reads[,2], ena_upload_dir)

#############################################################################################
##Helper function: cl_as_df

# Inspect a selected checklist as an R data frame. Includes the fields, units, vocabularies, and whether mandatory or not. 
#Requires xml to df conversion...only tested for ERC000022 

cl_as_df<-function(select.cl){
require("XML")
require(httr)
fileURL <- paste0("https://www.ebi.ac.uk/ena/browser/api/xml/",select.cl)
cl_xml <- GET(fileURL)
cl_xml_r <-xmlParse(cl_xml)
#cl_df<-xmlToDataFrame(nodes = getNodeSet(cl_xml_r, "//FIELD"))


fields <- getNodeSet(cl_xml_r,"//FIELD")
labels<-sapply(fields, function(x) xpathSApply(x, ".//LABEL", xmlValue))
description<-sapply(fields, function(x) xpathSApply(x, ".//DESCRIPTION", xmlValue))
mandatory<-sapply(fields, function(x) xpathSApply(x, ".//MANDATORY", xmlValue))
units<-sapply(fields, function(x) xpathSApply(x, ".//UNIT", xmlValue))
units<-unlist(lapply(units,toString))
value<-sapply(fields, function(x) xpathSApply(x, ".//VALUE", xmlValue))
value<-unlist(lapply(value,toString))

df <- data.frame(LABEL=labels,DESCRIPTION=description,MANDATORY=mandatory,UNIT=units,VALUE=value,stringsAsFactors=F)
return(df)
}

#example use
#cl_df_022<-cl_as_df(select.cl="ERC000022")


##############################################################################################
#function: build_sample_df

#creates a sample_info dataframe, which will need populating manually based on user env file

build_sample_df<-function(sample_alias,tax_id="256318", scientific_name="metagenome",select.cl,all_fields=TRUE,
sample_title=sample_alias,...){

df_core<-data.frame("sample_alias"=sample_alias,"tax_id"=tax_id,"scientific_name"=scientific_name,
"sample_title"=sample_title,stringsAsFactors=FALSE,...)

df_cl<-cl_as_df(select.cl)
if(all_fields){
df_cl_use<-df_cl
}
else{
df_cl_use<-df_cl[df_cl$MANDATORY=="mandatory",]
}

#add checklist fields to the core dataframe
df_core[,df_cl_use$LABEL]=""

#make units fields
units.f<-ifelse(df_cl_use$UNIT!="",paste0("UNIT_",df_cl_use$LABEL),"")
df_core[, units.f[units.f!=""]]=""

message(cat(""))
message(cat("Copy this into editor to help populate the samp_info data.frame:"))
message(cat(""))
message(cat(paste0("samp_info","$",'"',df_cl_use$LABEL,'"',"<-","env$","               #",df_cl_use$MANDATORY),sep = "\n"))
message(cat(paste0("samp_info","$",'"',units.f[units.f!=""],'"',"<-","            #unit options:",df_cl_use[units.f!="",]$UNIT),sep = "\n"))
return(df_core)


}

#example use:
#samp_info<-build_sample_df(sample_alias=env$"Sample ID for sequencing",select.cl="ERC000022",all_fields=FALSE)


################################################################################################
#function: create_sample_xml

#converts populated samp_info file to xml for programmatic ena submission

create_sample_xml<-function(samp_info,center_name,ERC_checklist){

subm_samples = newXMLDoc()
root = newXMLNode("SAMPLE_SET", doc = subm_samples)
# WRITE XML NODES AND DATA

#simple to do nodes first
for (i in 1:nrow(samp_info)){
samp = newXMLNode("SAMPLE",attrs =c(alias = samp_info$sample_alias[i],center_name=center_name),parent = root)
tit = newXMLNode("TITLE",samp_info$sample_title[i],parent = samp)  
samp_name = newXMLNode("SAMPLE_NAME",parent = samp)  
tax_id = newXMLNode("TAXON_ID",samp_info$tax_id[i],parent = samp_name) 
sci_nm = newXMLNode("SCIENTIFIC_NAME",samp_info$scientific_name[i],parent = samp_name) 
com_nm = newXMLNode("COMMON_NAME",samp_info$common_name[i],parent = samp_name) 

#now add other attributes
#units are an issue. We added "UNITS_" as columns to those fields which have specified units.
#we dont want these fields in the output xml (but we do want the units!)

#which columns are these
selcols<-names(samp_info[ , grepl( "UNIT_" , names( samp_info ) ) ])

#make df referencing fields that have units, matched to fields :
unit.df<-data.frame(unit_name=selcols,var=gsub(pattern="UNIT_", replacement = '', selcols))

#make df without these "UNITS_" column fields
samp_info_no_units<-samp_info[,!names(samp_info)%in%selcols]

#create attributes node
samp_attrs = newXMLNode("SAMPLE_ATTRIBUTES",parent = samp) 
#loop over and populate new attributes for each column
for (j in 6:ncol(samp_info_no_units)){
samp_attr<-newXMLNode("SAMPLE_ATTRIBUTE",parent = samp_attrs)
tag<-newXMLNode("TAG",names(samp_info_no_units)[j],parent = samp_attr) 
val<-newXMLNode("VALUE",samp_info_no_units[i,j],parent = samp_attr) 
if(names(samp_info_no_units)[j]%in%unit.df$var){
unit<-newXMLNode("UNITS",samp_info[i,paste0("UNIT_",names(samp_info_no_units)[j])]
,parent = samp_attr) }


}
samp_attr<-newXMLNode("SAMPLE_ATTRIBUTE",parent = samp_attrs)
tag<-newXMLNode("TAG","ENA-CHECKLIST",parent = samp_attr) 
val<-newXMLNode("VALUE",ERC_checklist,parent = samp_attr) 
}

return(subm_samples)
}

#example use
#samp_info<-cl_as_df(select.cl="ERC000022")



#########################################################################################
#function: make_reads_info

make_reads_info<-function(sample_alias,f_names,r_names,md5f,md5r,instrument_model="Illumina MiSeq",
library_source="METAGENOMIC",library_selection="PCR",library_strategy="AMPLICON",
insert_size="300",library_name="16S"){

dfs<-data.frame("sample_alias"=sample_alias,"instrument_model"=instrument_model,
"library_source"=library_source,"library_selection"=library_selection,
"library_strategy"=library_strategy,"insert_size"=insert_size,"forward_file_name"=f_names,
"forward_file_md5"=md5f,"reverse_file_name"=r_names,"reverse_file_md5"=md5r,"library_name"=library_name,stringsAsFactors=FALSE )

dfs$exp_alias<-paste(sample_alias,library_name,"exp",sep="_")
dfs$run_alias<-paste(sample_alias,library_name,"run",sep="_")

return(dfs)
}

############################################################################
#function: create_exp_xml
create_exp_xml<-function(reads_info,proj_name,title){

subm_exp = newXMLDoc()
root = newXMLNode("EXPERIMENT_SET", doc = subm_exp)
# WRITE XML NODES AND DATA

for (i in 1:nrow(reads_info)){

exp = newXMLNode("EXPERIMENT",attrs =c(alias = reads_info$exp_alias[i]),parent = root)
tit = newXMLNode("TITLE",title,parent = exp)  
proj = newXMLNode("STUDY_REF",attrs =c(refname  = proj_name),parent = exp)  
des = newXMLNode("DESIGN",parent = exp) 
des_desc<-newXMLNode("DESIGN_DESCRIPTION",parent = des)
samp_desc<-newXMLNode("SAMPLE_DESCRIPTOR",attrs =c(refname = reads_info$sample_alias[i]),parent = des)
lib_desc <- newXMLNode("LIBRARY_DESCRIPTOR",parent = des) 
lib_name<-newXMLNode("LIBRARY_NAME",reads_info$library_name[i],parent = lib_desc)
lib_strat<-newXMLNode("LIBRARY_STRATEGY",reads_info$library_strategy[i],parent = lib_desc)
lib_sour<-newXMLNode("LIBRARY_SOURCE",reads_info$library_source[i],parent = lib_desc)
lib_sel<-newXMLNode("LIBRARY_SELECTION",reads_info$library_selection[i],parent = lib_desc)
lib_lay<-newXMLNode("LIBRARY_LAYOUT",parent = lib_desc)
lib_typ<-newXMLNode("PAIRED",attrs =c(NOMINAL_LENGTH=reads_info$insert_size[i],NOMINAL_SDEV="30"),parent = lib_lay)
lib_pro<-newXMLNode("LIBRARY_CONSTRUCTION_PROTOCOL",parent = lib_desc)
plat = newXMLNode("PLATFORM",parent = exp) 
ill = newXMLNode("ILLUMINA",parent = plat)
mod<-newXMLNode("INSTRUMENT_MODEL",reads_info$instrument_model[i],parent=ill)
#not going to add any other experiment attributes, but if wanted to:
#exp_atts<-newXMLNode("EXPERIMENT_ATTRIBUTES",parent=exp)
#then add any additional fields with tag and values as per the samp_info xml


}

return(subm_exp)
}


##############################################################################
#function: create_run_xml

create_run_xml<-function(reads_info, center_name){


subm_run = newXMLDoc()
root = newXMLNode("RUN_SET", doc = subm_run)
# WRITE XML NODES AND DATA

for (i in 1:nrow(reads_info)){

run = newXMLNode("RUN",attrs =c(alias = reads_info$run_alias[i],center_name=center_name),parent = root)
exp_ref = newXMLNode("EXPERIMENT_REF",attrs =c(refname = reads_info$exp_alias[i]),parent = run)  
db = newXMLNode("DATA_BLOCK",parent = run)  
files = newXMLNode("FILES",parent = db) 
file<-newXMLNode("FILE", attrs =c(filename = reads_info$forward_file_name[i],
filetype="fastq",checksum_method="MD5",checksum=reads_info$forward_file_md5[i]),parent = files)
file<-newXMLNode("FILE", attrs =c(filename = reads_info$reverse_file_name[i],
filetype="fastq",checksum_method="MD5",checksum=reads_info$reverse_file_md5[i]),parent = files)
}


return(subm_run)
}


#########################################################################################
#function:create_proj_xml

create_proj_xml<-function(proj_name,proj_title,proj_desc){
require(XML)
project = newXMLDoc()
root = newXMLNode("PROJECT_SET", doc = project)
# WRITE XML NODES AND DATA
proj = newXMLNode("PROJECT",attrs =c(alias = proj_name),parent = root)
tit = newXMLNode("TITLE",proj_title,parent = proj)
desc=newXMLNode("DESCRIPTION",proj_desc,parent = proj)
sub_proj<-newXMLNode("SUBMISSION_PROJECT",parent = proj)
seq_proj<-newXMLNode("SEQUENCING_PROJECT",parent = sub_proj)
return(project)
}

###########################################################################################
#function:create_proj_subm_xml

create_proj_subm_xml<-function(rel_date){
	subm_proj = newXMLDoc()
	root = newXMLNode("SUBMISSION", doc = subm_proj)
	# WRITE XML NODES AND DATA
	actions = newXMLNode("ACTIONS",parent = root)
	action = newXMLNode("ACTION",parent = actions)
	add=newXMLNode("ADD",parent = action)
	action_1 = newXMLNode("ACTION",parent = actions)
	hold<-newXMLNode("HOLD",attrs =c(HoldUntilDate= rel_date),parent = action_1)
return(subm_proj)
}

########################################################################################
#function: create_SRE_subm_xml
	
#create submision xml #actions the submission, can be used to register SAMPLES,RUNS, and EXPERIMENTS

create_SRE_subm_xml<-function(){
     subm_samp = newXMLDoc()
	root = newXMLNode("SUBMISSION", doc = subm_samp)
	# WRITE XML NODES AND DATA
	actions = newXMLNode("ACTIONS",parent = root)
	action = newXMLNode("ACTION",parent = actions)
	add=newXMLNode("ADD",parent = action)
return(subm_samp)
	}

#########################################################################################
#function: RE_receipt_to_df

#parse xml RUN and EXP ("RE") receipt to obtain run and experiment accession numbers

RE_receipt_to_df<-function(RE_receipt_xml){

doc_xml <- xmlParse(RE_receipt_xml)
run_nodes<- getNodeSet(doc_xml, "//RUN")
exp_nodes<- getNodeSet(doc_xml, "//EXPERIMENT")

run_acc <- lapply(run_nodes, xmlGetAttr, "accession")
exp_acc <- lapply(exp_nodes, xmlGetAttr, "accession")

run_alias<-lapply(run_nodes, xmlGetAttr, "alias")
exp_alias<-lapply(exp_nodes, xmlGetAttr, "alias")

df<-data.frame( run_acc = unlist(run_acc),run_alias = unlist(run_alias),exp_acc=unlist(exp_acc),exp_alias=unlist(exp_alias))
return(df)
}


#########################################################################################
#function: samp_receipt_to_df

#parse xml sample ("samp") receipt to obtain sample accession numbers

samp_receipt_to_df<-function(samp_receipt_xml){

doc_xml <- xmlParse(samp_receipt_xml)
samp_nodes<- getNodeSet(doc_xml, "//SAMPLE")
ext_id_nodes<- getNodeSet(doc_xml, "//EXT_ID")

samp_acc <- lapply(samp_nodes, xmlGetAttr, "accession")
ext_id_acc <- lapply(ext_id_nodes, xmlGetAttr, "accession")

samp_alias<-lapply(samp_nodes, xmlGetAttr, "alias")


df<-data.frame( samp_acc = unlist(samp_acc),samp_alias = unlist(samp_alias),ext_id_acc=unlist(ext_id_acc))
return(df)
}
