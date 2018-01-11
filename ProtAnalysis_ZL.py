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


def extract_plnames_from_file(infile = "Plasmids.gb"):
    # Given a file of genbank records, try to strip the plasmid names
    # outputs into a csv that has the name, accession number and can be modified to have the Inc group
    # Automates finding the plasmid name in the genbank record
    # In places where the name cannot be called well, the csv output can be modified to manually call the correct name
    # Output:
    #   csv file with the following columns:  AccNum, Plasmid name, genus, species
    
    outfname = infile.strip(".gb")+".csv"
    #print "Writing names to %s." % outfname
    #print "Records that do not have the term plasmid in their name will be dropped."
    csv_file = csv.writer(open(outfname,"wb"),delimiter=',')

    count = 0  # a counter for the number of plasmids in the file
    total = 0  # count the total number of sequences in the file
    #print 'Extracting names from file, requires that the name of the plasmid follows the word plasmid in the genbank description'
	
    # This will be used when adding inc groups -ZL
    incGroupDict = make_incGroup_dict()

    with open(infile,"rU") as inhandle:
        for sq in SeqIO.parse(inhandle,"genbank"):
            #print sq.description
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
                        #print "^^"+str(sq.description.split("plasmid "))
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
                #print "%%%%%"+name
                
                count+=1
                ##print count
            else:
                continue
				
	# Get the inc-group for this plasmid using the list made above -ZL
            group = get_incGroup_from_plasmid(name,incGroupDict)
            row = [AccNum,name,Genus,Species,group]
            csv_file.writerow(row)
            total+=1
    #print "Found %i seqs in expected list of %i." % (count, total)

def make_incGroup_dict(file="IncGroups.csv"):
	with open(file,"rU") as inhandle:
		csv_file = csv.reader(inhandle,delimiter=',')
		groupDict = {row[0].lower() : row[1].lower() for row in csv_file}
		return groupDict
		
def get_incGroup_from_plasmid(plasmidName, incGroupDict):
	incGroup = ""
	try:
		incGroup = incGroupDict[plasmidName.lower()]
	except KeyError:
		incGroup = "Inc?"
	return incGroup

#========================================================================
# Fixed for Windows 10 OS (5/18/2017) -ZL
#========================================================================
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
    
    seqdict = []
    try:
        seqdict = {x[0]:x for x in GetSeqsAndGroupsFromCSV(incsv)}
        print "SUCCESS"
    except IndexError:
        print seqdict
        print GetSeqsAndGroupsFromCSV(incsv)[0:4]
        print "indErrorInStrip"
    outfa = ingb.strip(".gb")+"AA.fa"
    outhandle = open(outfa,"w")
    #print "Total records in csv file %i", len(seqdict)
    #print "Stripping all CDS's from seqs"
    
    k=0 #count the number of records
    m=0 # count the number of features
    storedSeqs = []
    storedSeqs_ws = []
    well_studied = ["pNDM-1_Dok01", "F", "RP4", "R751", "pRA3", "R7K", "pSK41"]
    bad = ["trans","recomb","integ","resist","tnp","tetracy","lactam","ins"]
    
    # "m" variable can give a unique identifier to each seq! -ZL (2017)
    
    for record in inseqs:
        #k+=1
        if record.name in seqdict:
            ##print seqdict[record.name]
            accNum=seqdict[record.name][0]
            plname=seqdict[record.name][1] # Save the plasmid name
            try:
                group = seqdict[record.name][4] # Save the group name
            except IndexError:
                group = ""
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
                            #print feat_name + " found on "+ record.id
                        if any_in_str(bad,feat_name.lower()):
                            continue
                        try:
                            sq=feature.extract(record)
                            # Using description to store name because the naems are too long with the positions in them
                            #sq.id = feat_name+"_"+plname
                            sq.id = feat_name+"_"+group+"_"+plname+"_"+str(m)
                            #sq.id = feat_name+"-"+str(j)+"_:_"+group+"_:_"+plname
                            # only right features that are nonempty (not sure why some are)
                            if len(sq)>0:
                                sq.seq = sq.seq.translate(table="Bacterial", to_stop=True)
                                if any_eq_str(well_studied, plname):
                                    storedSeqs_ws.append(sq)
                                else:
                                    storedSeqs.append(sq)
                                #SeqIO.write(sq,outhandle,"fasta")
                        except ValueError:
                            print "Referenced another sequence..."
                except KeyError:
                        print record.id + ' has no gene features, skipping sequence'
                        
    #Sort ALL sequences??
    storedSeqs.sort(key=lambda x: len(x.seq))
    storedSeqs_ws.sort(key=lambda x: len(x.seq))
    
    for storedSeq in reversed(storedSeqs + storedSeqs_ws):
        #print(len(storedSeq.seq))
        #ignore sequences less than 50 by request of Celeste -ZL
        if len(storedSeq.seq) > 50:
            SeqIO.write(storedSeq,outhandle,"fasta")
                    
        #else:
        #    #print record.name+" not found in csv"
    outhandle.close()

    #print "%s sequences in analysis" % k
    #print "%s features for analysis" % m

#========================================================================
# Fixed for Windows 10 OS (5/18/2017) - ZL
#========================================================================

def GetSeqsAndGroupsFromCSV(infile):
    ### Reads csv file with Acc number, name and Inc Group
    # Assumes table format AccNum, Name, Genus, Species, Group(optional)
    # Keeps track of Inc Group if known, otherwise assigned IncUKN
    # Returns list, where each entry has this info
    #print "Reading sequences and groups from CSV file, check that file contains updated Inc groups"
    with open(infile,"rU") as inhandle:
        csv_file = csv.reader(inhandle,delimiter=',')
        out = [row for row in csv_file]
        # for some reason this line causes rows to be skipped and I can't find the pattern
        #out = out[1:][::2] # Slicing: to eliminate blank lines on Windows OS
        print out

        return out

# Returns true if any item in the list is found IN the string
def any_in_str(listChecks, str):
    for item in listChecks:
        if item in str:
            return True
    return False
    
# Returns true if any item in the list EQUALS the string
def any_eq_str(listChecks, str):
    for item in listChecks:
        if item == str:
            return True
    return False

# Just like any_eq_str, but returns the first matching STR
def which_eq_str(listChecks, str):
    for item in listChecks:
        if item == str:
            return item
    return "NONE"

#==============================================================================
# Get clusters into a suitable format for R
# June 2017 edition
# Zac Lindsey
#==============================================================================
def gather_clusters_case_ZL(proteinFile = "Backbones_4.csv", 
                      clusterFile = "Plasmids20-200kb-6-9-2016_Clusters.tab",
                      outfile = "ClusterGroups.csv",
                      heatmap = False):
    import csv

    cRows = [] # Rows which contain a cluster summary
    sRows = [] # Rows which contain a centroid or hit 
    
    with open(clusterFile,"r") as clusterfile:
        readTab = csv.reader(clusterfile, delimiter="\t")
        for clustRow in readTab:
            if clustRow[0] == "C":
                cRows.append(clustRow)
            else:
                sRows.append(clustRow)
    
    titles = []
    Names = []
    ClusterIDs = []
    ProtInfo = {}
    index = -1
    bad = ["trans","recomb","integ","resist","tnp","tetracy","lactam","ins"]
    well_studied = ["pNDM-1_Dok01", "F", "RP4", "R751", "pRA3", "R7K", "pSK41"]
    
    with open(proteinFile,"r") as backbones:
        backboneFile = csv.reader(backbones)
        for backboneRow in backboneFile:
            titles.append(backboneRow[0])
            Names.append({})
            ClusterIDs.append(set())
            index += 1
            for protein in backboneRow: 
                for srow in sRows:
                    currProt = srow[8]
                    cid = srow[1]
                    ProtInfo[currProt] = [srow[2],srow[3],0]
                    if protein in currProt:
                        ClusterIDs[index].add(cid)
                    ProtInfo[currProt][2] = which_eq_str(well_studied,currProt.split("_")[2])
                    if cid in ClusterIDs[index]:
                        try:
                            Names[index][cid].add(currProt)
                        except KeyError:
                            Names[index][cid] = set()
                            Names[index][cid].add(currProt)   
    
    CG = []
    with open(outfile, "w") as out:
        Writer = csv.writer(out)
        for i in range(len(Names)):
            for j in Names[i]:
                rows = []
                for k in Names[i][j]:
                    row = []
                    row.append(titles[i])   # Append 4-letter name
                    row.append(j)   # Append Cluster ID
                    crPrInfo = k.split("_")
                    protName = crPrInfo[0]
                    incGroup = crPrInfo[1]
                    plasName = crPrInfo[2]
                    uniqueID = crPrInfo[-1]
                    row.append(protName)
                    row.append(incGroup)
                    row.append(plasName)
                    row.append(uniqueID)
                    row.append(ProtInfo[k][0])  # Append Seq Length
                    row.append(ProtInfo[k][1])  # Append Identity
                    row.append(ProtInfo[k][2])  # Append exemplar identifier
                    rows.append(row)
                rows.sort(key=lambda x: int(x[6]))
                for row in reversed(rows):
                    Writer.writerow(row)
                    CG.append(row)
    return CG

#========================================================================
# Method for generating heat map in R (1/11/2018) - ZL
# CG is essentially a list of rows. each row has a certain format.           

#========================================================================
def Pre_heat_map_table(CG, outfile):
     
    Labels = set([row[0] for row in CG])
    plasNames = set([row[4] for row in CG])
    
    # w, h = len(set(plasNames)), len(set(Labels));
    # Matrix = {[0 for x in range(w)] for y in range(h)}
    Matrix = {x:{y:0 for y in Labels} for x in plasNames}
    
    for plasName in plasNames:
        LabelsByPlasmid = []
        for row in CG:
            if row[4] == plasName:
                LabelsByPlasmid.append(row[0])
        for Label in set(LabelsByPlasmid):
            Matrix[plasName][Label] = 1
      
    with open(outfile, "w") as out:
        Writer = csv.writer(out)
        Writer.writerow([key for key in Matrix['1'].keys()])
        for plasName in plasNames:
            Writer.writerow([plasName]+[Matrix[plasName][key] for key in Matrix[plasName].keys()])
            
#    LabelsByPlasmid = {plasmid:[] for plasmid in plasNames}
#    for plasName in plasNames:
#        for row in CG:
#            if row[4] == plasName:
#                LabelsByPlasmid[plasName].append(row[0])
#                
#    with open(outfile, "w") as out:
#        for key in LabelsByPlasmid.keys():
#            Writer = csv.writer(out)
#            Writer.writerow( [key]+ list(set( LabelsByPlasmid[key] )) )
            
##########################################################

#c("4LetterName","cid","protName","IncGroup","Plasmid","uid","SeqLength","Identity","Ex")

###########################################################

if __name__=="__main__":

    proj = "Plasmids20-200kb-6-9-2016"
    #extract_plnames_from_file(proj+".gb")
#==============================================================================
# Had to modify code to work on a MAC
#  there should be ~4349 plasmids and ~328435 sequences, but without the update there were only ~2300 plasmids
#==============================================================================
    #StripCDSAddInc(proj+".gb",proj+".csv")
    
    #os.system('usearch.exe -cluster_fast '+proj+'AA.fa -id 0.7 -target_cov 0.7 -centroids clust_'+proj+'.fasta -uc '+proj+'_Clusters.tab -userfields query+target+id+ql+tl+alnlen+qcov+tcov')
    CG = gather_clusters_case_ZL(proteinFile = "Backbones_6-13-2017.csv",
                      clusterFile = "Plasmids20-200kb-6-9-2016_Clusters.tab",
                      outfile = "ClusterGroups_1-10-2018_1.csv")
    Pre_heat_map_table(CG, "matrixData_1-11-2018_1.csv")