
#########################################
                                  
      Compressor for aligned reads 
       (described in a SAM file)              

#########################################
		Contact               
		
Please contact the authors for any bug/comment:

Idoia Ochoa:	iochoa@stanford.edu
Mikel Hernaez:	mhernaez@stanford.edu

#########################################
		Purpose

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

1) Compress the reads contained in the SAM file
2) Reconstruct the reads using the output of the compressor and the reference used for the alignment (to generate the SAM file)

#########################################
		Download

1) Download the software from github. (https://github.com/mikelhernaez/cbc) 
2) Unzip it.

#########################################
              Installation

Go to the main folder and run the following command:

$ make

This will create a folder named "bin" where you will find the binary to run the program.

#########################################
		Usage

For an example on how to run the program and generate the necessary files, please see example below.

1) Compress the reads contained in the SAM file:

$ PATH_TO_BINARY/cbc -c <samFile> <outputFile> <ReferenceFile>

The <samFile> is a SAM file ordered by position, that contains only the aligned reads, i.e., those reads that failed to align to the reference should not be included in the SAM file. We further assume the AUX field MD is present.
The <outputFile> is the name of the compressed file will have.
The <ReferenceFile> is the FASTA file used as reference for the creation of the SAM file.

3) Decompress the reads:

$ PATH_TO_BINARY/cbc -d <compressedFile> <outputFile> <ReferenceFile>

The <outputFile> is the file where the program is going to write the reconstructed reads.
The <compressedFile> is the file already compressed by cbc.
The <ReferenceFile> is the FASTA file used as reference for the creation of the SAM file.

#########################################
    Thanks for using our program!!               
#########################################
