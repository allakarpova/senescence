#January 20 2022
#Reyka Jayasinghe
#reyka@wustl.edu


#Usage: sc_isoform_mapping.py -b sample.bam -c coordinates -r barcodes
#-r barcodes is not needed to run - currently doesnt have functionality

#Coordinates file 
#Example file - two exons junctions of interest to scan bam data to identify 
#chr	start	stop	exon	isoformtype
#chr1	198692373	198703298	exonA	PTPRC-RO
#chr1	198692373	198699564	exonA	PTPRC-RB+RBC
#chr1	198692373	198696712	exonA	PTPRC
#chr1	198699704	198702387	exonA	PTPRC-RBC
#chr1	198699704	198703298	exonA	PTPRC-RB
#exonA is a placeholder for functionality in the future

import re
import sys
import getopt
import argparse
from subprocess import call
from subprocess import check_output
from subprocess import check_call
import subprocess
from collections import defaultdict
import os
import shlex


def main(argv):
   bamfile = ''
   coordinates = ''
   barcodes = ''
   try:
      opts, args = getopt.getopt(argv,"hc:b:r:",["bfile=","coords=","bcode="])
   except getopt.GetoptError:
      print('python sc_isoform_mapping.py -b <bamfile> -c <coordinates> -r <barcodes>')
      print('-r barcodes is not needed to run - currently doesnt have functionality')
      sys.exit(3)
   for opt, arg in opts:
      if opt == '-h':
         print('python sc_isoform_mapping.py -b <bamfile> -c <coordinates> -r <barcodes>')
         sys.exit()
      elif opt in ("-b", "--bfile"):
         bamfile = arg
      elif opt in ("-c", "--coords"):
         coordinates = arg
      elif opt in ("-r","--bcode"):
      	 barcodes = arg
      #return(bamfile,coordinates,barcodes)

if __name__ == "__main__":
   main(sys.argv[1:])

position={}
exon={}
junction_bc={}
bc={}
bc_dict={}
gene_coord=open(sys.argv[4])
bamfile=sys.argv[2]
#Import coordinates data
#Be sure that your chromosome coordinates for bam are exactly the same as how they show up in the bam file
#to check chromosome annotation: samtools view -H bam 
#check "SN: " options include chr1,1,GRCh38_chr1,mm10_chr1 
for cline in gene_coord:
	#remove header from coordinates file
		(chr,start,stop,e,isoform)=cline.strip().split('\t')
		key=chr+"__"+start+"__"+stop
		position[key]=isoform
		exon[key]=e	
		junction_bc[key]=[] #create empty list to store barcodes
gene_coord.close()

for k in position:
    (c,start,stop)=k.strip().split('__')
    search_junction=c+":"+start+"-"+stop
    #Search bam file for any junction reads that overlap coordinates for exon/intron boundaries
    samtools_cmd="samtools view "+bamfile+" "+search_junction
    #print(samtools_cmd)
    #process = check_call(shlex.split(samtools_cmd),shell=False) #works but pritns output to terminal
    process=subprocess.check_output(shlex.split(samtools_cmd),shell=False).decode(sys.stdout.encoding).strip()
    reads=process.split("\n")
    junction_count=0
    total_reads=0
    for read in reads:
    	total_reads+=1
    	read_details=read.strip().split("\t")
    	if not re.match("\tCB:",read): #only keep reads that have cell barcode details
    		#Determine if read spans exon junction
    		#Reference: https://www.biostars.org/p/89581/
    		#read_details[5] 68M53005N27M3S
    		if re.match("^(\d+)M(\d+)N(\d+)M$",read_details[5]):
    			#if read_details[0] == "A00585:97:HKWTJDSXX:4:1124:6379:3959":
	    			#print(read)
	    		m=re.match("^(\d+)M(\d+)N(\d+)M$",read_details[5])
	    		exon1_match=m.group(1)
	    		intron_m=m.group(2)
	    		exon2_match=m.group(3)
	    			#Calculate exon/intron boundary of aligned read
	    			###APPLIES FOR BOTH POSITIVE AND NEG STRAND DATA
	    		eib1=(int(read_details[3])+int(exon1_match))-int(1)
	    		eib2=(int(eib1)+int(intron_m))+int(1)
	    		chr_suffix=search_junction.split(':')
	    		check_junction=chr_suffix[0]+"__"+str(eib1)+"__"+str(eib2)
	    		if check_junction in position:
	    			#cellbarcode=read_details[24].split(":")[2]
	    			for element in read_details:
	    				if re.match("^CB:",element):
	    					cellbarcode=element.split(":")[2]
	    					junction_bc[check_junction].append(cellbarcode)
	    					isoformtype=position[check_junction]
	    					if cellbarcode in bc_dict:
	    						bc_dict[cellbarcode].append(isoformtype)
	    					else:
	    						bc_dict[cellbarcode]=[]
	    						bc_dict[cellbarcode].append(isoformtype)
	    					#bc_list[cellbarcode].append(isoformtype)
	    				if re.match("^RG:",element):
	    				#Extract sample info from bam
	    				#RG:Z:MC200_NT-6:0:1:HKWTJDSXX:4
	    				#sampleinfo=read_details[15].split(":")[2]
	    					sampleinfo=element.split(":")[2]

outfile=sampleinfo+"_junctions_bc.txt"
fh = open(outfile, "w")
#header="Junction"+"\t"+"Sample"+"\t"+"Barcode"+"\t"+"Exon"+"\t"+"Isoform"+"\t"
header="Sample"+"\t"+"Barcode"+"\t"+"Isoform"
fh.write(header+"\n")

for i in bc_dict:
	itype=bc_dict[i]
	unique_isoforms=set(itype)
	info=sampleinfo+"\t"+i+"\t"+";".join(unique_isoforms)+"\n"
	fh.write(info)

#for i in junction_bc:
#	isoform_type=position[i] #determine gene annotation
#	exon_number=exon[i] #determine gene annotation
#	list_bc=junction_bc[i]
#	for bc_i in set(list_bc): #some barcodes have multiple reads - take unique of list
#		info=i+"\t"+sampleinfo+"\t"+bc_i+"\t"+exon_number+"\t"+isoform_type+"\n"
#		fh.write(info)