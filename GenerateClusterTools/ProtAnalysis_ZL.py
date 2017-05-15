from Bio.Seq import Seq, MutableSeq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature,FeatureLocation
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.SeqUtils import GC
from Bio import SeqIO, Entrez, Alphabet
from Bio.Graphics import GenomeDiagram
from Bio.Graphics.GenomeDiagram import CrossLink
from reportlab.lib.colors import *
from reportlab.lib import colors
import urllib2
import os
import csv
import numpy as np

#from user_plasmids import *
#from expanding import *
#from copy import deepcopy

###########################################################


### Tools for downloading and retrieving sequences

def generateRecords_ZL(input="pcHitTable.csv"):
    #We will use a hit table for each result
    inhandle = open(input,"rU")
    csv_file = csv.reader(inhandle,delimiter=',',dialect=csv.excel)
    #List of gi nums for the Entrez retrieval
    gi_num = []
    for row in csv_file:
        if len(row)>0: #Prevents error on empty rows
            gi_num.append(row[1].split("|")[3]) #This just adds the gi number
    Entrez.email = "ryanbotts@pointloma.edu"
    print "Fetching nucleotide sequences..."
    handle=Entrez.efetch(db='nucleotide',id=gi_num,rettype='gb')
    print "Retrieval complete."
    records = []
    print("Parsing records...")
    count = 0
    for Record in SeqIO.parse(handle,"genbank"):
        records.append(Record)
        count+=1
        print "Appending ("+str(count)+")..."
    print "{0}"+str(records[0:3])
    return records
        

### General utilities for handling lists of sequences, retrieving records from a list, extracting plasmid names, downloading records from a list, etc.
def get_Target_Names(infile="ISCR1extAlign.gb"):
    # get all of the record id's for sequences in a multi-genbank file
    # built to allow a user to manually remove records from the analysis
    in_handle = open(infile,"rU")
    inseqs = SeqIO.parse(in_handle,"genbank")
    print "Using name stored in description."
    # Full name stored in the record because it was too long to store as the id or name.
    return [record.description.split(' ')[0] for record in inseqs]


def extract_plnames_from_file(infile = "Plasmids.gb"):
    # Given a file of genbank records, try to strip the plasmid names
    # outputs into a csv that has the name, accession number and can be modified to have the Inc group
    # Automates finding the plasmid name in the genbank record
    # In places where the name cannot be called well, the csv output can be modified to manually call the correct name
    # Output:
    #   csv file with the following columns:  AccNum, Plasmid name, genus, species
    
    outfname = infile.strip(".gb")+".csv"
    print "Writing names to %s." % outfname
    print "Records that do not have the term plasmid in their name will be dropped."
    csv_file = csv.writer(open(outfname,"w"),delimiter=',',dialect=csv.excel)

    count = 0  # a counter for the number of plasmids in the file
    total = 0  # count the total number of sequences in the file
    print 'Extracting names from file, requires that the name of the plasmid follows the word plasmid in the genbank description'

    with open(infile,"rU") as inhandle:
        for sq in SeqIO.parse(inhandle,"genbank"):
            print sq.description
            total+=1
            Genus=sq.description.split(" ")[0]
            Species=sq.description.split(" ")[1]
            AccNum=sq.name
            if sq.description.find('plasmid')>0:
                # logic takes the name after the word plasmid if it exists, otherwise takes the name before it.
                #prior = sq.description.split("plasmid")[0].split(" ")
                
                #name = sq.description.split("plasmid")[1].split(" ")[1].strip(",").strip(".")
                try:
                    if len(sq.description.split('plasmid ')[1].split(' ')) > 0:
                        print "^^"+str(sq.description.split("plasmid "))
                        name = sq.description.split('plasmid ')[1].split(' ')[0].strip(",") # take the word after plasmid as the name
                    elif len(sq.description.split('plasmid ')[1].split(' ')) == 0:
                        name = sq.description.split('plasmid')[0].split(' ')[0] # take the word after plasmid as the name
                except IndexError:
                    if sq.description.find('plasmid.'):
                        name = "unnamed"
                    else:
                        if len(sq.description.split('plasmid,')[0].split(' ')) > 0:
                            stuff = sq.description.split('plasmid,')[1].split(' ')
                            name = stuff[-1].strip(",") # take the word after plasmid as the name
                print "%%%%%"+name
                
                count+=1
                #print count
            else:
                continue
            row = [AccNum,name,Genus,Species]
            csv_file.writerow(row)
            total+=1
    print "Found %i seqs in expected list of %i." % (count, total)

def StripCDSAddInc(ingb,incsv):
    # routine opens a file with plasmid names and inc groups, if known, along with gb file with seqs
    # strips all features from each genbank record found in the csv file and saves the amino acid sequence, named with inc group
    # open all genbank records in ingb and reads all names and groups in the csv
    # naming convention of genes is:
    #   "gene_name"-"number"_:_"IncGroup"_:_"PlasmidName"
    # Inc Group is set as IncUnk if it is not known
    # Inputs:
    #   ingb-input file of all of the genbank records for analysis
    #   incsv-input file with the plasmid name, inc group, and accession number, used to name the features appropriately
    #
    # Outputs:
    #   outfa-fasta output of amino acid sequences of all of the genes, named accordingly, file is named the same as the input
    #           ingb with the appendix changed to AA.fa.
    inseqs = SeqIO.parse(open(ingb,"rU"),"genbank")
    
    seqdict = {x[0]:x for x in GetSeqsAndGroupsFromCSV(incsv)}
    outfa = ingb.strip(".gb")+"AA.fa"
    outhandle = open(outfa,"w")
    print "Total records in csv file %i", len(seqdict)
    print "Stripping all CDS's from seqs"
    
    k=0 #count the number of records
    m=0 # count the number of features
    for record in inseqs:
        #k+=1
        if record.name in seqdict:
            #print seqdict[record.name]
            plname=seqdict[record.name][1] # Save the plasmid name
            k+=1
            j = 0 # count the number of orfs on each plasmid
            for feature in record.features:
                try:
                    if feature.type == "CDS":
                        m+=1
                        j+=1
                        # name the feature based on the gene or product name.  If neither is available, set to orf something
                        if 'product' in feature.qualifiers:
                            feat_name=feature.qualifiers['product'][0].replace(" ","")
                        elif 'gene' in feature.qualifiers:
                            feat_name=feature.qualifiers['gene'][0].replace(" ","")
                        else:
                            feat_name=str(['orfX'+str(j)]).replace(" ","")
#                        print feat_name + " found on "+ record.id
                        try:
                            sq=feature.extract(record)
                            # Using description to store name because the naems are too long with the positions in them
                            sq.id = feat_name+"_"+plname
                            #sq.id = feat_name+"-"+str(j)+"_:_"+group+"_:_"+plname
                            # only right features that are nonempty (not sure why some are)
                            if len(sq)>0:
                                sq.seq = sq.seq.translate(table="Bacterial", to_stop=True)
                                SeqIO.write(sq,outhandle,"fasta")
                        except ValueError:
                            print "Referenced another sequence..."
            
                except KeyError:
                    print record.id + ' has no gene features, skipping sequence'
        #else:
        #    print record.name+" not found in csv"
    outhandle.close()

    print "%s sequences in analysis" % k
    print "%s features for analysis" % m


def extract_seqs_from_file(infile = "ISCR1Keeps.csv",targetfile = "ISCR1Align.gb", outseqfile = "ISCR1AlignSubset"):
    # Utility for selecting a subset of sequences from a file with multiple records
    # Given a file of seqs and a csv file with names, extract all sequences in the file with those names
    # csv file must have the seqs listed in the first column
    # Inputs:
    #   infile-csv file with all of the names of the desired sequences, searches by accession number
    #   targetfile-fasta or genbank file with all records
    #   outseqfile-name of the file for outputting only the extracted sequences
    # Outputs:
    #   outseqfile-all of the records retrieved from the larger set of records
    inhandle = open(infile,"rU")
    csv_file = csv.reader(inhandle,delimiter=',',dialect=csv.excel)
    seqnames =  [row[0] for row in csv_file]
    inhandle.close()

    seqhandle = open(targetfile,"rU")
    
    # identify the type of reads we are using
    if targetfile.find("gb") > 0:
        type = "genbank"
        outname = outseqfile + '.gb'
    elif targetfile.find("fasta") > 0:
        type = "fasta"
        outname =  outseqfile + '.fasta'

    count = 0
    print 'Extracting features from file.  Sequence name must be a subset of names in'
    outseqhandle = open(outname, 'w')
    
    for sqname in seqnames:
        print "Searching for %s" % sqname
        for sq in SeqIO.parse(open(targetfile,"rU"),type):
        
            if type == "genbank":
                sqnm = sq.description.split(' ')[0]
            
            else:
                sqnm = sq.id
            
            if sqnm in sqname:
                print "record id %s" % sq.id
    
    #for sq in SeqIO.parse(seqhandle,type):
    #    print sq.id
    #    if any(sq.id in sqname for sqname in seqnames):
                SeqIO.write(sq, outseqhandle, type)
                count+=1
    outseqhandle.close()

    print "Found %i seqs in expected list of %i." % (count, len(seqnames))


#========================================================================
# Note that this method had to be modified to work on a MAC, approx 2000 rows of the table were not read in
def GetSeqsAndGroupsFromCSV(infile):
    ### Reads csv file with Acc number, name and Inc Group
    # Assumes table format AccNum, Name, Genus, Species, Group(optional)
    # Keeps track of Inc Group if known, otherwise assigned IncUKN
    # Returns list, where each entry has this info
    print "Reading sequences and groups from CSV file, check that file contains updated Inc groups"
    with open(infile,"rU") as inhandle:
        csv_file = csv.reader(inhandle,delimiter=',',dialect=csv.excel)
        #Reads csv file and adds unknown group if that column is unknown
        try:
            out = [row if row[4]!= "" else row[0:4] +["IncUNK"] for row in csv_file]
        except IndexError:
            csv_file = csv.reader(inhandle,delimiter=',',dialect=csv.excel)
            out = [row for row in csv_file]

        # for some reason this line causes rows to be skipped and I can't find the pattern
        #out = out[1:][::2] # Slicing: to eliminate blank lines on Windows OS

        return out




### Tools for handling the output of USEARCH clustering
def get_clusters(clusterfilename):
    ### finds clusters from default tab delimited output from USEARCH
    # creates multiple lists each containing all items in a cluster
    # input:
    #    clusterfilename - tab delimited default output from USEARCH
    #        center of cluster begins with a S in the first column and name in the 9th column
    # output:
    #     groups - list containing all elements of each cluster
    print "May not correctly extract clusters of singletons"
    groups = []  # empty dictionary to store each cluster
    
    with open(clusterfilename, 'rU') as f:
        reader = csv.reader(f, delimiter='\t')
        notfirst = False # identifies if we are starting at the first group or not
        group=[]
        for row in reader:
            if row[0].strip() == "S" and notfirst:
                # add previous group to dictionary
                groups.append(group)
                # start a new group
                group=[]
            # the rows marked with a C indicate the center of a cluster and may be omitted.
            elif row[0].strip() == "C":
                groups.append(group)
                break
            notfirst = True
            name = row[8].strip() # extract the name from the row
            group.append(name)
        # append the last cluster to the list
        groups.append(group)
    return groups

def get_clust_name(pro_list):
    ### Attempts to determine an appropriate name for all proteins in a cluster
    # logic needs improvement,
    # names the group with the following priorities
    #   1.  a protein in the list has one of the preferred names, then that is used first
    #   2.  the first 4 letter protein name in the list
    #   3.  the first 6 (for 4 letter names with a number) letter named protein in the list
    #   4.  the first protein name in the list that is not hypothetical or orf
   
    #print 'Creating protein family names, attempting to name the cluster'
    #print 'Uses the following names if possible:  sul1, qacEdelta, blaCTX-M,aac, aad4'
    
    # keep track of how long the name is and whether it is the best we have or not.
    min_len_name = 0
    found_best = False
    
    pref_names = ['sul1','qacEdelta','blaCTX-M','aac','aad4','ISCR', 'tnpR','tnpA','ins1','mer', 'TEM']
    
    # default name is the first name in the group of nothing better is found
    group_name = pro_list[0]

    # take the next item satsfying the following criteria, if one is not found, procede to the next statement
    # if none of these criteria works, use default name

    group_name = next((x for x in pro_list for y in pref_names if y in x.split('_:_')[0]),\
        next((x for x in pro_list if len(x.split('_:_')[0].split('-')[0]) == 4),\
        next((x for x in pro_list if len(x.split('_:_')[0].split('-')[0]) == 6),\
            group_name)
        ))

    return group_name







##########################################################


###########################################################

if __name__=="__main__":

    ### tools for downloading ISCR sequences
    #Find_Sequence("M1-ISCR1.fasta")
    proj_name = 'Plasmids2016'
    #1. Done through BLAST, download all sequences in BLAST output
    
   # retrieve_from_blast_table(proj_name+'HitTable.csv', proj_name+'Complete.gb', proj_name+'Align.fasta',proj_name+'ExtAlign.gb', 1015, proj_name, True)
    #retrieve_from_gb_and_blast_table(infile = proj_name+'HitTable.csv', ingb = proj_name+'Complete.gb', outfile = proj_name+'Align.fasta', extoutfile = proj_name+'ExtAlign.gb', targetlen = 513, featname = proj_name, down_only = True)
    ##extract_seqs_from_file(targetfile = proj_name+'Align.gb', outseqfile = proj_name+'AlignSubset')
    #infile = proj_name+'Keeps.csv',

    # note that if usearch fails for empty file it is likely due to having written a duplicate sequence in the AA file, remove sequence
    #2. Extract all proteins from all records
   # StripCDS(infile=proj_name+'.gb', outfile = proj_name+'FeaturesAA.fasta',nucoutfile = proj_name+'FeaturesNN.fasta')

#==============================================================================
# #2. Will need to be adapted for either Windows or Mac
#==============================================================================
    #2. Cluster proteins into protein families
  #  os.system("usearch.exe -cluster_fast "+proj_name+"FeaturesAA.fasta -id 0.7 -target_cov 0.7 -centroids clust.fasta -uc "+proj_name+"Clusters.tab -userfields query+target+id+ql+tl+alnlen+qcov+tcov")

    #3. identigy protein families on each sequence and export into matrix
   # write_matrix_file(proj_name+'FeatureMatrixShortNames.csv',targetnamefile = proj_name+'ExtAlign.gb',clusterfile = proj_name+'Clusters.tab',writefullnames=False)

    #4. Find the patterns in the gene combinations in R
    #5. Map the groups to visualize patters in the protein combinations
 #   map_seq_groups(outfilename=proj_name+'MapGroup', plasmidftfile = proj_name+'FeatureMatrixShortNames.csv', plasmidclustfile = proj_name+'1Groups.csv', plasmidseqs=proj_name+'ExtAlign.gb')
    #        seqclusters = proj_name+'Clusters.tab', pathtoplasmids='./', featname = proj_name)

    # prior to running matrices through heat map, removed all proteins found in 2 or fewer sequences.

    proj = "Plasmids20-200kb-6-9-2016"
    #extract_plnames_from_file(proj+".gb")
#==============================================================================
# Had to modify code to work on a MAC
#  there should be ~4349 plasmids and ~328435 sequences, but without the update there were only ~2300 plasmids
#==============================================================================
    #StripCDSAddInc(proj+".gb",proj+".csv")
    os.system('./usearch -cluster_fast '+proj+'AA.fa -id 0.5 -target_cov 0.7 -centroids clust_'+proj+'.fasta -uc '+proj+'_Clusters.tab -userfields query+target+id+ql+tl+alnlen+qcov+tcov')
    #os.system('./usearch -cluster_fast clust_'+proj+'.fasta -id 0.5 -target_cov 0.5 -centroids '+proj+'_2nd.fasta -uc '+proj+'Clusters_2nd.tab -userfields query+target+id+ql+tl+alnlen+qcov+tcov')

