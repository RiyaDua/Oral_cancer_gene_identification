# using trimmomatic tool
java -jar /home/pallavis/pooja/trimommatic/Trimmomatic-0.36/trimmomatic-0.36.jar PE /home/pallavis/pooja/sra/sratoolkit../bin/SRR21496995_1.fastq /home/pallavis/ pooja/sra/sratoolkit../bin/SRR21496995_2.fastq SRR21496995_p_1.fastq unpaired_cut_1_SRR21496995.fastq SRR21496995_p_2.fastq unpaired_cut_2_SRR21496995.fastq HEADCROP:3 MINLEN:10 


# using star aligner
./STAR --runMode alignReads --outFileNamePrefix output83 --genomeDir sta-index --readFilesIn /home/pallavis/pooja/sra/sratoolkit.3.0.0-ubuntu64/bin/SRR11080783_p_1.fastq /home/pallavis/pooja/sra/sratoolkit.3.0.0-ubuntu64/bin/SRR11080783_p_2.fastq --outSAMtype BAM SortedByCoordinate --runThreadN 4 

##Locating bamfile
file <- list.files(pattern = "\\.bam$")

####Defining gene models
gtffile1 <- file.path("Homo_sapiens.GRCh38.gtf")

###Counting with featureCounts
library ("Rsubread”)
fc <- featureCounts(files=file,
                    annot.ext=gtffile1,
                    isGTFAnnotationFile=TRUE,
                    isPairedEnd=TRUE) 
