#!/usr/bin/python2
#
#Reformats UCSC knownGene table dump
#
#INPUT
#   stdin: UCSC knownGene table in BED format
#	argv[1]: id list to filter

import csv
import sys
		
def main():

    uid = []
    with open(sys.argv[1], 'r') as idFile:
        idTable = csv.reader(idFile, delimiter='\t')
        for record in idTable:
            uid.append(record[0])
    uid = frozenset(uid)
        
    reader = csv.reader(sys.stdin, delimiter='\t')

    for record in reader:

        if record[0].startswith('#'):
            continue
        
        #UCSC BED coordinates are 0-based
        chrom = record[0]
        txStart = int(record[1])
        txEnd = int(record[2])
        name = record[3]
        strand = record[5]
        cdsStart = int(record[6])
        cdsEnd = int(record[7])
        exonCount = int(record[9])
       	exonSizes = (record[10].split(','))
    	exonStarts = (record[11].split(','))
    	exonSizes = map(int, exonSizes[0:(len(exonSizes)-1)])
    	exonStarts = map(int, exonStarts[0:(len(exonStarts)-1)])
    	         
        #Skip if not in ID list
        if not (name in uid):
            continue
            
        #Output null for records that have 0 CDS length?
#        if cdsStart==cdsEnd:
#            sys.stdout.write('\t'.join([name, 'NULL', '\n']))        
#            continue

        cdsStart = cdsStart -txStart
        cdsEnd = cdsEnd - txStart
            
        splicedLength = sum(exonSizes)
        
        #Filter out exons that come after cds start/end positions
        #Exon fields should be sorted
        LeftUTRExonStarts = [exonStart for exonStart in exonStarts if exonStart <= cdsStart]
        RightUTRExonStarts = [exonStart for exonStart in exonStarts if exonStart <= cdsEnd]
        
        #Add up all the exon lengths up to cdsStart/end
        leftUTR = sum(exonSizes[0:(len(LeftUTRExonStarts)-1)]) + cdsStart - exonStarts[len(LeftUTRExonStarts)-1]
        rightUTR = sum(exonSizes[0:(len(RightUTRExonStarts)-1)]) + cdsEnd - exonStarts[len(RightUTRExonStarts)-1]
        cdsLength = rightUTR - leftUTR
   
        #OUTPUT

        if strand=="-":
            leftUTR = splicedLength - leftUTR
            rightUTR = splicedLength - rightUTR
#                if (rightUTR == 0) or (leftUTR == splicedLength):
#                    continue
            if cdsLength==0:
                sys.stdout.write('\t'.join([name,'0',str(splicedLength), "NOCDS", '0', '+', '\n']))
            else:
                sys.stdout.write('\t'.join([name,'0',str(splicedLength), "TX", '0', '+', '\n']))
                if rightUTR!=0:
                    sys.stdout.write('\t'.join([name, '0', str(rightUTR), "5UTR", '0', '+', '\n']))
                if leftUTR!=splicedLength:
                    sys.stdout.write('\t'.join([name, str(leftUTR), str(splicedLength), "3UTR", '0', '+', '\n']))
                sys.stdout.write('\t'.join([name, str(rightUTR), str(leftUTR), "ORF", '0', '+', '\n']))

                if cdsLength>60:
                    sys.stdout.write('\t'.join([name, str(rightUTR+45), str(leftUTR-15),"ORF_COUNT",'0','+','\n']))

        else:
#                if (leftUTR == 0) or (rightUTR == splicedLength):
#                    continue
            if cdsLength==0:
                sys.stdout.write('\t'.join([name,'0',str(splicedLength), "NOCDS", '0', '+', '\n']))
            else:
                sys.stdout.write('\t'.join([name,'0',str(splicedLength), "TX", '0', '+', '\n']))
                if leftUTR!=0:
                    sys.stdout.write('\t'.join([name, '0', str(leftUTR), "5UTR", '0', '+', '\n']))
                if rightUTR!=splicedLength:
                    sys.stdout.write('\t'.join([name, str(rightUTR), str(splicedLength), "3UTR", '0', '+', '\n']))
                sys.stdout.write('\t'.join([name, str(leftUTR), str(rightUTR), "ORF", '0', '+', '\n']))

				#START+45nt, END-15nt for effective ORF
                if cdsLength>60:
                    sys.stdout.write('\t'.join([name, str(leftUTR+45), str(rightUTR-15),"ORF_COUNT",'0','+','\n']))


if __name__ == '__main__':
    main()
