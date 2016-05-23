#test file
#!/bin/bash -x 

#************---Pray it runs---***********
#                 (' - ')
#                 |  0 0__
#               .'     ___O      
#              /   .-.' \\
#             J    |`.\\  \\
#             | |_.|  | | |
#              \\__.'`.|-' /
#              L      `--'\\
#              |           | 
#              J           \\                              
#              \\         /  \\
#               \\      .'`.  `.                                 .'
#              ___) /\ (___`.  `-._____________________________.'/
#            _///__/__\___\\\\_`-.____________________________.-'
#
#************---Pray it runs---***********

##Final project - Chip-seq
#Geno Villafano

Base_run="off"
	#***********---Define---***********
		#Path to raw data files (fastq)
		inpath="/home/gvillafano/../../archive/MCB5429/Final_data/ChIPseq/"
			#Makes a list of files within the inpath directory but skips directories
			fastqs=$(ls -p $inpath | grep -v '/$')
			#If only using certain files in list indicated with $1 $2 $3 to choose file 1,2,3 or $fastqs for all files
			set -- $fastqs
			fastqs=$(echo $1 $2 $3 $4 $5 $6) 
		#Path to FOLDER of outpath-highgest branch up to which output files will  be
		outpath="/home/gvillafano/Final_project/ChIPseq/outfiles/"
		#Pattern used to subset the alignment to certain regions/chromsomes
		pattern="chr19"
		#Adaptor sequence to remove
		adaptor="GATCGGAAGAGCTCGTATGCCGTCTTCTGCTTGAAA"
		#Genome index to be aligned to
		genome="/home/gvillafano/../../archive/MCB5429/genomes/hg19/bowtieIndex/hg19"
		#Chromosome size file
		chrsizes="/home/gvillafano/Final_project/helper_files/hg19_chromInfo.txt"
	
	#***********---What to run---***********
	QC_and_trimming="off"
		QC_and_trimming_fastqc="on"
		QC_and_trimming_fastx_clipper="on"
		QC_and_trimming_fastx_trimmer="on"

	Align_to_genome="off"
		Align_to_genome_alignment="on"
		Align_to_genome_extract="on"
		Align_to_genome_SamtoBed="on"
		Align_to_genome_sort="on"
		Align_to_genome_Bedgraph="on"
			Bedgraph_separate_strands="off"
			Bedgraph_add_tracklines="on"
				Tracklines_colour_treated="on"
	
Analysis="on"
	#***********---Define---***********
		#Looks within the outpath from the previous set of functions for the processed data
		an_inpath=${outpath}
		#A name specifier as to which are controls/untreated (Needed for MACS)
		controls="untr"
		treated="treat"
		#Input here the size of the chromosome or region being extracted
		regionsize=59128983
		#Name of region or chromosome name	
		regionName="${pattern}" 
		#Set the area around summits to capture window when looking for motifs
		#I.e. windows=25 will give a +/- 50 bp window around each peak location
		windows="50"
		#After getting summits with MACS, subset out top_peaks=<integer> number of peaks
		top_peaks="200"
		#fasta formatted chromosomes
		faChr="/home/gvillafano/../../archive/MCB5429/genomes/hg19/fasta/"

	#***********---What to run---***********
	MACS_peaks="on"
		MACS_run="on"
			MACS_summit_slop="on"
				MACS_subset_slopped_summits="on"
			MACS_subset_summits="off"
				MACS_MEME_run="on"
					MACS_MEME_MAST="on"
		
#***********---Show defined paths and files---***********	
echo "inpath being used: ${inpath}"
echo "outpath being used: ${outpath}"
echo -e "inpath to be used:\n ${inpath}"
echo -e "All the files within inpath dir\n:::: ${fastqs}" | tee -a log_file.txt

#$$$***********---BASE RUN---***********$$$

if [ $Base_run == "on" ]
	then
			
	for file in ${fastqs}
		do
	
	
		#***********---Setup file structure---***********
		#Sets variable as inpath + fastq file currently being used within loop
		inPATH=${inpath}${file} 
		echo "File currently being worked on: ${file}" | tee -a log_file.txt

		#Sets variable to be used as base file name for saving outputs without extension
		outPATH=${outpath}${file%.*} 
		echo "outPATH to be used:"
		echo "${outPATH}"
		outFILE=${outPATH}/${file%.*}

		#Makes separate output directory for each file in loop
		mkdir ${outPATH}/
	
	
		#***********---Quality checks and trimming---***********
		if [ $QC_and_trimming == "on" ]
			then
				#Runs QC analysis on data. Generates two new files for each fastq being checked
				if [ $QC_and_trimming_fastqc == "on" ]
					then	
				
					echo "${file} - Fastqcin' it"  | tee -a log_file.txt
					fastqc -o ${outPATH} ${inPATH} | tee -a log_file.txt
				fi
		
				#Removes adaptor sequence from reads as defined in $adaptor variable, generates a file called {outPATH}_c.tmp		
				if [ $QC_and_trimming_fastx_clipper	== "on" ]
					then
				
					echo "${file} - Clippin' it" 
					echo -e "Adaptor sequence: \n ${adaptor} " | tee -a log_file.txt
					fastx_clipper -v -Q33 -i ${inPATH} -o ${outFILE}_c.tmp -a ${adaptor} | tee -a log_file.txt
				fi

				#Takes _c.tmp file and trims bases bellow a quality score of 20, should be mostly ends, makes {outPATH}_c1.temp
				if [ $QC_and_trimming_fastx_trimmer == "on" ]
					then
				
					echo "${file} - Trimmin' it"  | tee -a log_file.txt
					fastq_quality_trimmer -v -i ${outFILE}_c.tmp -Q33 -t 20 -l 5 -o ${outFILE}_c1.tmp | tee -a log_file.txt
			
					rm ${outFILE}_c.tmp  #Removes the clipped file, as only to carry on with trimmed file
				fi
		fi		
	
	
		#***********---Align to genome(s)---***********
		if [ $Align_to_genome == "on" ]
			then
				#Aligns trimmed file to genome, makes 3 files: noalign, yesalign, _align.log (bowtie stats) 
				if [ $Align_to_genome_alignment == "on" ]
					then
				
					echo "Alignin' ${file} to the genome: ${genome}" | tee -a log_file.txt
					bowtie -p8 -v2 -m1 -q --sam --un ${outFILE}_noalign ${genome} ${outFILE}_c1.tmp ${outFILE}_yesalign.sam 2>&1 | tee ${outFILE}_align.log
				
					rm ${outFILE}_noalign #Removes reads that didn't align to save space
				fi
						
				#Extracting only certain chromosome reads
				if [ $Align_to_genome_extract == "on" ]
					then
				
					echo "Extractin' reads just for ${pattern}"  | tee -a log_file.txt
					samtools view ${outFILE}_yesalign.sam | grep ${pattern} > ${outFILE}_yesalign_${pattern}.sam | tee -a log_file.txt
					#Adds header to the grepped out sam file, overwrites the header-less file afterwards
					samtools view -H ${outFILE}_yesalign.sam | cat - ${outFILE}_yesalign_${pattern}.sam > temp && mv temp ${outFILE}_yesalign_${pattern}.sam  | tee -a log_file.txt
				
					rm ${outFILE}_yesalign.sam #Removes original sam file to save space
				
				fi
			
				#Convertin' SAM to BAM to BED
				if [ $Align_to_genome_SamtoBed == "on" ]
					then

					echo "${file} - Convertin' SAM to BAM to BED" | tee -a log_file.txt
					samtools view ${outFILE}_yesalign_${pattern}.sam -Sb  >  ${outFILE}_yesalign_${pattern}.bam | tee -a log_file.txt
					bamToBed -i ${outFILE}_yesalign_${pattern}.bam > ${outFILE}_yesalign_${pattern}.bed | tee -a log_file.txt
				
					rm ${outFILE}_yesalign_${pattern}.sam #Removing original sam file to save space
				
				fi				

				#Sortin' the BED
				if [ $Align_to_genome_sort == "on" ]
					then
				
					echo "${file} - Sortin' the BED" | tee -a log_file.txt
					sort -k 1,1 -k2,2n ${outFILE}_yesalign_${pattern}.bed > ${outFILE}_yesalign_${pattern}_sorted.bed | tee -a log_file.txt
				
					rm ${outFILE}_yesalign_${pattern}.bed #Removing unsorted bed to save space
				fi
	
	
			#***********---Generating Bedgraphs and tracks---***********
			if [ $Align_to_genome_Bedgraph == "on" ]
				then
				
				echo "Creatin' Bedgraphs for ${outFILE}_yesalign_${pattern}_sorted.bed"  | tee -a log_file.txt
				for chromosome in ${chrsizes}
					do
		
				#Separate strands or not
				if [ $Bedgraph_separate_strands == "on" ]
					then
						
					echo "Separatin' the strands of ${outFILE} to two separate bedgraphs"  | tee -a log_file.txt
					awk '$6 == "+"' ${outFILE}_yesalign_${pattern}_sorted.bed | genomeCoverageBed -5 -i /dev/stdin -bg -g ${chromosome} > ${outFILE}_yesalign_${pattern}_sorted_plus.bedgraph 
					awk '$6 == "-"' ${outFILE}_yesalign_${pattern}_sorted.bed | genomeCoverageBed -5 -i /dev/stdin -bg -g ${chromosome} > ${outFILE}_yesalign_${pattern}_sorted_m.bedgraph
					awk '{ $4=$4*-1; print }' ${outFILE}_yesalign_${pattern}_sorted_m.bedgraph > ${outFILE}_yesalign_${pattern}_sorted_minus.bedgraph #adds minus to read properly			

					rm ${outFILE}_yesalign_${pattern}_sorted_m.bedgraph
				
					else
				
					echo "Makin' ${outFILE} into a standard bedgraph"  | tee -a log_file.txt
					genomeCoverageBed -5 -i ${outFILE}_yesalign_${pattern}_sorted.bed -bg -g ${chromosome} > ${outFILE}_yesalign_${pattern}_sorted.bedgraph 
				fi 
				done
		
				#Adding Tracklines
				if [ $Bedgraph_add_tracklines == "on" ]
					then
			
					if [ $Tracklines_colour_treated == "on" ]
						then
			
						if [[ ${outFILE} == *"treat"*  ]] #searching path of input file (${outFILE}) for the word treat and changes the colour if found
							then
					
							echo "Addin' some different colour tracklines to ${file}"  | tee -a log_file.txt
							awk 'BEGIN {  print "browser position ${pattern}:5,289,521-5,291,937"
										  print "track type=bedGraph name=\"$\" description=\"${outFILE}\"  visibility=full autoScale=on alwaysZero=on color=250,20,0 windowingFunction=maximum"} 
									   {  print $0}' ${outFILE}_yesalign_${pattern}_sorted.bedgraph  > ${outFILE}_yesalign_${pattern}_sorted_header.bedgraph

							else
				
							echo "Addin' tracklines to ${file}"  | tee -a log_file.txt
							awk 'BEGIN {  print "browser position ${pattern}:5,289,521-5,291,937"
									  print "track type=bedGraph name=\"$\" description=\"${outFILE}\"  visibility=full autoScale=on alwaysZero=on color=0,253,0 windowingFunction=maximum"} 
								   {  print $0}' ${outFILE}_yesalign_${pattern}_sorted.bedgraph  > ${outFILE}_yesalign_${pattern}_sorted_header.bedgraph
						fi					
			
						else
				
						echo "Addin' tracklines to ${file}"  | tee -a log_file.txt
						awk 'BEGIN {  print "browser position ${pattern}:5,289,521-5,291,937"
									  print "track type=bedGraph name=\"$\" description=\"${outFILE}\"  visibility=full autoScale=on alwaysZero=on color=0,253,0 windowingFunction=maximum"} 
								   {  print $0}' ${outFILE}_yesalign_${pattern}_sorted.bedgraph  > ${outFILE}_yesalign_${pattern}_sorted_header.bedgraph
					fi
				fi			
			fi
		

		fi		
	done
fi


#$$$***********---ANALYSIS---***********$$$

if [ $Analysis == "on" ]
	then	
	
	
	echo -e "\n\n\n\n\nStarting Analysis of data\n\n\n\n\n" | tee -a ${an_inpath}${an_inpath}${an_inpath}analysis_log_file.txt
				
	#***********---Calling MACS---***********
	if [ $MACS_peaks == "on" ]
		then
		
		echo -e "Running MACS-based analyses\n"  | tee -a ${an_inpath}analysis_log_file.txt
		#Creates new folder for MACS outputs
		an_outpath="${an_inpath}MACS_output"
		mkdir ${an_outpath}
		
		if [ $MACS_run == "on" ]
			then
				
			echo -e "Calling MACS on path: $an_inpath \n output directed to ${an_outpath}"	| tee -a ${an_inpath}analysis_log_file.txt
			#Creates array of paths to treated samples
			an_paths_treat=($(find ${an_inpath} -maxdepth 2 -name "*${treated}*" -name "*${regionName}*" -name "*.bam*" -print | sort))
			echo -e "Treated samples: ${an_paths_treat[*]}\n\n"  | tee -a ${an_inpath}analysis_log_file.txt
			#Creates array of paths to controls 
			an_paths_con=($(find ${an_inpath} -maxdepth 2 -name "*${controls}*" -name "*${regionName}*" -name "*.bam*"  -print | sort))
			echo -e "Control samples: ${an_paths_con[*]}\n\n"  | tee -a ${an_inpath}analysis_log_file.txt
			
			echo -e "Using MACS to call peaks on ${an_paths_treat[*]}\n\n" | tee -a ${an_inpath}analysis_log_file.txt
			for i in $(seq 0 $(expr ${#an_paths_treat[*]} - 1)) #for each of the treated files...
				do

				#Stores dir from which script was run
				current_wd=$(pwd)
				#Sets output directory because I can't figure out how to redirect MACS
				cd ${an_outpath}
			
				#Getting file name without extension
				the_file="$(basename ${an_paths_treat[ $i ]})"
				file_name="${the_file%.*}"
				echo "Calling peaks for ${file_name} with control data ${an_paths_con[ $i ]}\n"	 | tee -a ${an_inpath}analysis_log_file.txt
				macs14 -t ${an_paths_treat[ $i ]} -c ${an_paths_con[ $i ]} -n ${file_name}  -g ${regionsize} 2>&1| tee -a ${an_inpath}analysis_log_file.txt
			
				#Change back to starting dir
				cd ${current_wd}
			done
		fi
			
		if [ $MACS_summit_slop == "on" ]
			then 
			
			echo -e "Sloppin' around the summits file +/- ${windows}\n"  | tee -a ${an_inpath}analysis_log_file.txt
			#Creating list of summit files
			summits=($(ls ${an_outpath}/*_summits.bed))
			echo -e "Found these summits files:\n ${summits[*]}\n"  | tee -a ${an_inpath}analysis_log_file.txt
			echo -e "Capturing +/- ${windows} window\n"	 | tee -a ${an_inpath}analysis_log_file.txt
		
			for summitf in ${summits[*]}
				do
			
				#Uses bedtools slop to create window
				bedtools slop -i ${summitf} -g ${chrsizes} -b ${windows} -header > ${summitf%.*}_slop.bed 2>&1| tee -a ${an_inpath}analysis_log_file.txt
				
				#Take only the top x-amount of peaks as defined by user
				if [ $MACS_subset_slopped_summits == "on" ]
					then
					
					#Creates new files with top # called peaks after sorting by column 5 (peak intensity)
					echo -e "Removing all but top ${top_peaks} peaks\n"  | tee -a ${an_inpath}analysis_log_file.txt
					sort -n -r -k5 ${summitf%.*}_slop.bed | head -${top_peaks} > ${summitf%.*}_slop_top${top_peaks}.bed
					rm ${summitf%.*}_slop.bed
				fi
				
				echo -e "Done sloppin' ${summitf}\n"  | tee -a ${an_inpath}analysis_log_file.txt
				done	
		fi	
		
		if [ $MACS_subset_summits == "on" ]
			then
			
			#Creating list of summit files
			summits=($(ls ${an_outpath}/*_summits.bed))
			echo -e "Found these summits files:\n ${summits[*]}\n"  | tee -a ${an_inpath}analysis_log_file.txt
			
			for summitf in ${summits[*]}
				do
				
				#Creates new files with top # called peaks after sorting by column 5 (peak intensity)
				echo -e "Removing all but top ${top_peaks} peaks\n"  | tee -a ${an_inpath}analysis_log_file.txt
				sort -n -r -k5 ${summitf} | head -${top_peaks} > ${summitf%.*}_slop_top${top_peaks}.bed 
				rm ${summitf}
				
			done
		fi	

		#***********---Calling MEME---***********
		if [ $MACS_MEME_run == "on" ]
			then
		
			#Makes new dir and saves current for saving between
			MEME_outpath="${an_outpath}/MEME_output"
			mkdir ${MEME_outpath}
			current_wd=$(pwd)

			#Getting fasta sequence for slopped summit files 
			echo -e "Checking for slopped summits files"  | tee -a ${an_inpath}analysis_log_file.txt
			slopped_summits=($(ls ${an_outpath}/*_summits_slop*.bed))
			echo -e "Found these slopped-summit files:\n ${slopped_summits[*]}\n"  | tee -a ${an_inpath}analysis_log_file.txt

			for slopped in ${slopped_summits[*]}
				do
				
				#gets just file name from end of path
				slopped_base=$(basename ${slopped})
				slopped_base_noext=${slopped_base%.*}
				#Moves into MEME output directory and makes new dir for each slopped summit file
				mkdir ${MEME_outpath}/${slopped_base%.*}
				cd ${MEME_outpath}/${slopped_base%.*}
				temp_out_dir=$(pwd)
			
				#Retrieves fasta sequences for summit regions
				echo -e "Gettin' fasta sequences for ${slopped_base}\n"  | tee -a ${an_inpath}analysis_log_file.txt
				bedtools getfasta -name -fi ${faChr}${pattern}.fa -bed ${slopped} -fo ${slopped_base_noext}.fa  | tee -a ${an_inpath}analysis_log_file.txt
			
				#Runs MEME based on fasta sequences and find motif
				echo -e "Runnin' MEME on ${slopped_base}\n"
				meme ${slopped_base_noext}.fa -oc MEME_out -dna -nmotifs 2 -minw 5 -maxw 10 -revcomp -mod zoops  2>&1| tee -a ${an_inpath}analysis_log_file.txt
			
				if [ $MACS_MEME_MAST == "on" ]
					then
				
					#Runs MAST to search if motif from MEME is within fasta regions
					echo -e "Running MAST for ${slopped_base}\n"  | tee -a ${an_inpath}analysis_log_file.txt
					mast "$(pwd)/MEME_out/meme.txt" ${slopped_base_noext}.fa -oc MAST_out  2>&1| tee -a ${an_inpath}analysis_log_file.txt
				fi
				cd ${current_wd}
			done
			
		fi
	fi
fi

#   And if it don't
# 
#    ____..(\,;,/)
#   '=====' (o o)\//,
#  ..//'_____\ /     \,
#  '.~-------`+'(  (   \    )
#               //  \   |_./
#               '~' '~----'     