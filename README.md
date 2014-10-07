
#########################################
                                  
      Compressor for aligned reads 
       (described in a SAM file)              
                                  
#########################################

#########################################
		Contact               
#########################################
Please contact the authors for any bug/comment:

Idoia Ochoa:	iochoa@stanford.edu
Mikel Hernaez:	mhernaez@stanford.edu

#########################################
		Purpose
#########################################

This is a program for compression and decompression of aligned reads presented in a SAM file.
Note that the purpose of this algorithm is to compress the necessary information to reconstruct the reads contained in the SAM file, i.e.,
we do not compress the identifiers and the quality values.
The decoder produces a file with only the reconstructed reads (i.e., we do not reconstruct the SAM file). 

Assumptions made by the program:
1) The SAM file is ordered by position.
If the sam file is not ordered, please convert it first to a bam file, order it, and then convert back to sam file (one can do this with SamTools, for example).
2) The SAM file contains ONLY mapped reads, i.e., those reads that did not map should not be included in the SAM file.
Recall that the purpose of the proposed method is to compress the aligned reads only.
3) The SAM file contains the AUX field MD.

#########################################
		Usage
#########################################

1) Compress the reads contained in the SAM file
2) Reconstruct the reads using the output of the compressor and the reference used for the alignment (to generate the SAM file)


#########################################
		Download
#########################################

1) Download the software from sourceForge. 
Choose between Ubuntu 64 bits and MAC distributions.
2) Unzip it: you will find the following folders and files

- ProposedMethod
	- README
	- simulations
		- main.c
		- arith.inc
		- input_file.inc 
		- os_stream.inc
		- sam_stream.inc
		- stats.inc
		- ceReference.txt 
	- SAMfiles
		- SRR065390_1_bowtie2_c_elegans_ws235.mapped.example.sam
	- c_elegans
		- CHROMOSOME_I.fa
	- compressedFiles
		- 
	- reconstructedReads
		-

#########################################
              Installation
#########################################

Go to the folder simulations, and run the following command:
gcc -o main.run main.c

Alternatively, one can use
main_MAC.run (if using MAC)
main_Ubuntu64.run (if using Ubuntu 64 bits)

#########################################
		Usage
#########################################

For an example on how to run the program and generate the necessary files, please see example below.

1) Compress the reads contained in the SAM file
Go to the folder simulations
Run the following command

./main.run <samFile> <prefixOfOutput>

The <samFile> is a SAM file ordered by position, that contains only the aligned reads, i.e., those reads that failed to align to the reference should not be included in the SAM file. We further assume the AUX field MD is present.
The <prefixOfOutput> is the path/prefix that the compressed files will have. In this case the program will generate the following files:
prefixOfOutput_F.ido
prefixOfOutput_I.ido
prefixOfOutput_M.ido
prefixOfOutput_P.ido
prefixOfOutput_P_A.ido
prefixOfOutput_S.ido
prefixOfOutput_char.ido
prefixOfOutput_pos.ido

3) Decompress the reads
Go to the folder simulations
Run the following command

./main.run <reconstructedReads> <fileWithPathToRefChromosomes> <prefixOfOutput>

The <reconstructedReads> is the file where the program is going to write the reconstructed reads.
The <fileWithPathToRefChromosomes> is a file where line X contains the path to the chromosome X (see example below)
The <prefixOutput> should be the same as the one specified in step 1) (for compression)

#########################################
		Example
#########################################

In the following example, we compress the SAM file SRR065390_1_bowtie2_c_elegans_ws235.mapped.example.sam, generated from the FASTQ file SRR065390_1.fastq with bowtie2 with c_elegans_ws235 as reference. For the purpose of the example, the sam file includes only around 3 million reads that mapped to the chromosome1. 

We then decompress it and reconstruct the reads.

1) Compress the sam file
Go to the folder simulations
Run the following command:

./main.run ../SAMfiles/SRR065390_1_bowtie2_c_elegans_ws235.mapped.example.sam ../compressedFiles/test

The command generates the files 
../compressedFiles/test_F.ido
../compressedFiles/test_I.ido
../compressedFiles/test_M.ido
../compressedFiles/test_P.ido
../compressedFiles/test_P_A.ido
../compressedFiles/test_S.ido
../compressedFiles/test_char.ido
../compressedFiles/test_pos.ido

2) Decompress the reads
Go to the folder simulations
Run the following command:

./main.run ../reconstructedReads/example.reads ceReference.txt ../compressedFiles/test

The file ceReference.txt contains in the first line the path to the first chromosome. If there were more chromosomes, they should be specified one per line. In this case we are just using the first chromosome.
The reconstructed reads will be generated in reconstructedReads/example.reads
The compressedFiles/test is the prefix of the compressed files generated in step 1).

#########################################
    Thanks for using our program!!               
#########################################

