import sys, getopt, os, random, math
from random import randint
from decimal import *
from collections import deque
from Bio import Phylo, AlignIO
from Bio.Phylo.Applications import PhymlCommandline
from Bio.File import as_handle
import sqlite3
conn = sqlite3.connect('gene')

class Organism:

    genome = ''
    dnaGenome = ''
    score = 0
    genelist = []
    isParent = False

    #initialize constructor
    def __init__(self, name):
        self.name = name
        self.genome = self.createGenome()
        #fix switch for abr and meta
        self.score = self.getScore(False, False)
        self.genelist = self.getGenes()
        self.dnaGenome = self.converttoDNA()

    #refresh variables when new generation is produced
    def refreshGenes(self, genome):
        self.genome = genome
        self.score = self.getScore(False, False)
        self.genelist = self.getGenes()
        self.dnaGenome = self.converttoDNA()

    #generate genome for initial populations
    def createGenome(self):
        c = conn.cursor()
        sequence = ''

        c.execute('SELECT BinaryEncoding FROM gene WHERE GeneName="clr" ORDER BY RANDOM() LIMIT 1')
        gene1 = str(c.fetchone())
        gene1 = gene1.replace(' ', '')[:-3].upper()
        gene1 = gene1.replace(' ', '')[3:].upper()
        sequence += gene1

        c.execute('SELECT BinaryEncoding FROM gene WHERE GeneName="clr" ORDER BY RANDOM() LIMIT 1')
        gene2 = str(c.fetchone())
        gene2 = gene2.replace(' ', '')[:-3].upper()
        gene2 = gene2.replace(' ', '')[3:].upper()
        sequence += gene2

        c.execute('SELECT BinaryEncoding FROM gene WHERE GeneName="face" ORDER BY RANDOM() LIMIT 1')
        gene3 = str(c.fetchone())
        gene3 = gene3.replace(' ', '')[:-3].upper()
        gene3 = gene3.replace(' ', '')[3:].upper()
        sequence += gene3

        c.execute('SELECT BinaryEncoding FROM gene WHERE GeneName="rand" ORDER BY RANDOM() LIMIT 1')
        gene4 = str(c.fetchone())
        gene4 = gene4.replace(' ', '')[:-3].upper()
        gene4 = gene4.replace(' ', '')[3:].upper()
        sequence += gene4

        c.execute('SELECT BinaryEncoding FROM gene WHERE GeneName="rand" ORDER BY RANDOM() LIMIT 1')
        gene5 = str(c.fetchone())
        gene5 = gene5.replace(' ', '')[:-3].upper()
        gene5 = gene5.replace(' ', '')[3:].upper()
        sequence += gene5

        c.execute('SELECT BinaryEncoding FROM gene WHERE GeneName="rand" ORDER BY RANDOM() LIMIT 1')
        gene6 = str(c.fetchone())
        gene6 = gene6.replace(' ', '')[:-3].upper()
        gene6 = gene6.replace(' ', '')[3:].upper()
        sequence += gene6

        c.execute('SELECT BinaryEncoding FROM gene WHERE GeneName="meta" ORDER BY RANDOM() LIMIT 1')
        gene7 = str(c.fetchone())
        gene7 = gene7.replace(' ', '')[:-3].upper()
        gene7 = gene7.replace(' ', '')[3:].upper()
        sequence += gene7

        c.execute('SELECT BinaryEncoding FROM gene WHERE GeneName="abr" ORDER BY RANDOM() LIMIT 1')
        gene8 = str(c.fetchone())
        gene8 = gene8.replace(' ', '')[:-3].upper()
        gene8 = gene8.replace(' ', '')[3:].upper()
        sequence += gene8

        c.execute('SELECT BinaryEncoding FROM gene WHERE GeneName="essn" ORDER BY RANDOM() LIMIT 1')
        gene9 = str(c.fetchone())
        gene9 = gene9.replace(' ', '')[:-3].upper()
        gene9 = gene9.replace(' ', '')[3:].upper()
        sequence += gene9

        c.execute('SELECT BinaryEncoding FROM gene WHERE GeneName="usls" ORDER BY RANDOM() LIMIT 1')
        gene10 = str(c.fetchone())
        gene10 = gene10.replace(' ', '')[:-3].upper()
        gene10 = gene10.replace(' ', '')[3:].upper()
        sequence += gene10

        return sequence


    def converttoDNA(self):
        test = ''
        temp1 = self.genome[0:2]
        temp2 = self.genome[2:4]
        temp3 = self.genome[4:6]
        temp4 = self.genome[6:8]
        temp5 = self.genome[8:10]
        temp6 = self.genome[10:12]
        temp7 = self.genome[12:14]
        temp8 = self.genome[14:16]
        temp9 = self.genome[16:18]
        temp10 = self.genome[18:20]
        temp11 = self.genome[20:22]
        temp12 = self.genome[22:24]
        temp13 = self.genome[24:26]
        temp14 = self.genome[26:28]
        temp15 = self.genome[28:30]
        temp16 = self.genome[30:32]
        temp17 = self.genome[32:34]
        temp18 = self.genome[34:36]
        temp19 = self.genome[36:38]
        temp20 = self.genome[38:40]
        temp21 = self.genome[40:42]
        temp22 = self.genome[42:44]
        temp23 = self.genome[44:46]
        temp24 = self.genome[46:48]
        temp25 = self.genome[48:50]
        temp26 = self.genome[50:52]
        temp27 = self.genome[52:54]
        temp28 = self.genome[54:56]
        temp29 = self.genome[56:58]
        temp30 = self.genome[58:60]

        if(temp1 == '00'):
            test += 'A'
        elif(temp1 == '01'):
            test += 'C'
        elif(temp1 == '10'):
            test += 'G'
        elif(temp1 == '11'):
            test += 'T'

        if(temp2 == '00'):
            test += 'A'
        elif(temp2 == '01'):
            test += 'C'
        elif(temp2 == '10'):
            test += 'G'
        elif(temp2 == '11'):
            test += 'T'

        if(temp3 == '00'):
            test += 'A'
        elif(temp3 == '01'):
            test += 'C'
        elif(temp3 == '10'):
            test += 'G'
        elif(temp3 == '11'):
            test += 'T'

        if(temp4 == '00'):
            test += 'A'
        elif(temp4 == '01'):
            test += 'C'
        elif(temp4 == '10'):
            test += 'G'
        elif(temp4 == '11'):
            test += 'T'

        if(temp5 == '00'):
            test += 'A'
        elif(temp5 == '01'):
            test += 'C'
        elif(temp5 == '10'):
            test += 'G'
        elif(temp5 == '11'):
            test += 'T'

        if(temp6 == '00'):
            test += 'A'
        elif(temp6 == '01'):
            test += 'C'
        elif(temp6 == '10'):
            test += 'G'
        elif(temp6 == '11'):
            test += 'T'

        if(temp7 == '00'):
            test += 'A'
        elif(temp7 == '01'):
            test += 'C'
        elif(temp7 == '10'):
            test += 'G'
        elif(temp7 == '11'):
            test += 'T'

        if(temp8 == '00'):
            test += 'A'
        elif(temp8 == '01'):
            test += 'C'
        elif(temp8 == '10'):
            test += 'G'
        elif(temp8 == '11'):
            test += 'T'

        if(temp9 == '00'):
            test += 'A'
        elif(temp9 == '01'):
            test += 'C'
        elif(temp9 == '10'):
            test += 'G'
        elif(temp9 == '11'):
            test += 'T'

        if(temp10 == '00'):
            test += 'A'
        elif(temp10 == '01'):
            test += 'C'
        elif(temp10 == '10'):
            test += 'G'
        elif(temp10 == '11'):
            test += 'T'

        if(temp11 == '00'):
            test += 'A'
        elif(temp11 == '01'):
            test += 'C'
        elif(temp11 == '10'):
            test += 'G'
        elif(temp11 == '11'):
            test += 'T'

        if(temp12 == '00'):
            test += 'A'
        elif(temp12 == '01'):
            test += 'C'
        elif(temp12 == '10'):
            test += 'G'
        elif(temp12 == '11'):
            test += 'T'

        if(temp13 == '00'):
            test += 'A'
        elif(temp13 == '01'):
            test += 'C'
        elif(temp13 == '10'):
            test += 'G'
        elif(temp13 == '11'):
            test += 'T'

        if(temp14 == '00'):
            test += 'A'
        elif(temp14 == '01'):
            test += 'C'
        elif(temp14 == '10'):
            test += 'G'
        elif(temp14 == '11'):
            test += 'T'

        if(temp15 == '00'):
            test += 'A'
        elif(temp15 == '01'):
            test += 'C'
        elif(temp15 == '10'):
            test += 'G'
        elif(temp15 == '11'):
            test += 'T'

        if(temp16 == '00'):
            test += 'A'
        elif(temp16 == '01'):
            test += 'C'
        elif(temp16 == '10'):
            test += 'G'
        elif(temp16 == '11'):
            test += 'T'

        if(temp17 == '00'):
            test += 'A'
        elif(temp17 == '01'):
            test += 'C'
        elif(temp17 == '10'):
            test += 'G'
        elif(temp17 == '11'):
            test += 'T'

        if(temp18 == '00'):
            test += 'A'
        elif(temp18 == '01'):
            test += 'C'
        elif(temp18 == '10'):
            test += 'G'
        elif(temp18 == '11'):
            test += 'T'

        if(temp19 == '00'):
            test += 'A'
        elif(temp19 == '01'):
            test += 'C'
        elif(temp19 == '10'):
            test += 'G'
        elif(temp19 == '11'):
            test += 'T'

        if(temp20 == '00'):
            test += 'A'
        elif(temp20 == '01'):
            test += 'C'
        elif(temp20 == '10'):
            test += 'G'
        elif(temp20 == '11'):
            test += 'T'

        if(temp21 == '00'):
            test += 'A'
        elif(temp21 == '01'):
            test += 'C'
        elif(temp21 == '10'):
            test += 'G'
        elif(temp21 == '11'):
            test += 'T'

        if(temp22 == '00'):
            test += 'A'
        elif(temp22 == '01'):
            test += 'C'
        elif(temp22 == '10'):
            test += 'G'
        elif(temp22 == '11'):
            test += 'T'

        if(temp23 == '00'):
            test += 'A'
        elif(temp23 == '01'):
            test += 'C'
        elif(temp23 == '10'):
            test += 'G'
        elif(temp23 == '11'):
            test += 'T'

        if(temp24 == '00'):
            test += 'A'
        elif(temp24 == '01'):
            test += 'C'
        elif(temp24 == '10'):
            test += 'G'
        elif(temp24 == '11'):
            test += 'T'

        if(temp25 == '00'):
            test += 'A'
        elif(temp25 == '01'):
            test += 'C'
        elif(temp25 == '10'):
            test += 'G'
        elif(temp25 == '11'):
            test += 'T'

        if(temp26 == '00'):
            test += 'A'
        elif(temp26 == '01'):
            test += 'C'
        elif(temp26 == '10'):
            test += 'G'
        elif(temp26 == '11'):
            test += 'T'

        if(temp27 == '00'):
            test += 'A'
        elif(temp27 == '01'):
            test += 'C'
        elif(temp27 == '10'):
            test += 'G'
        elif(temp27 == '11'):
            test += 'T'

        if(temp28 == '00'):
            test += 'A'
        elif(temp28 == '01'):
            test += 'C'
        elif(temp28 == '10'):
            test += 'G'
        elif(temp28 == '11'):
            test += 'T'

        if(temp29 == '00'):
            test += 'A'
        elif(temp29 == '01'):
            test += 'C'
        elif(temp29 == '10'):
            test += 'G'
        elif(temp29 == '11'):
            test += 'T'

        if(temp30 == '00'):
            test += 'A'
        elif(temp30 == '01'):
            test += 'C'
        elif(temp30 == '10'):
            test += 'G'
        elif(temp30 == '11'):
            test += 'T'

        return test

    #given a genome and environmental variables, determine an organisms fitness score
    def getScore(self, switch1, switch2):
        temp_score = 0
        c = conn.cursor()
        genescore1 = self.genome[0:6]
        genescore2 = self.genome[6:12]
        genescore3 = self.genome[12:18]
        genescore4 = self.genome[18:24]
        genescore5 = self.genome[24:30]
        genescore6 = self.genome[30:36]
        genescore7 = self.genome[36:42]
        genescore8 = self.genome[42:48]
        genescore9 = self.genome[48:54]
        genescore10 = self.genome[54:60]

        c.execute('SELECT score1 FROM gene WHERE BinaryEncoding =?', (genescore1,))
        temp = str(c.fetchone())
        temp = temp.replace(' ', '')[:-2].upper()
        temp = temp.replace(' ', '')[1:].upper()
        temp_score += int(temp)

        c.execute('SELECT score1 FROM gene WHERE BinaryEncoding =?', (genescore2,))
        temp = str(c.fetchone())
        temp = temp.replace(' ', '')[:-2].upper()
        temp = temp.replace(' ', '')[1:].upper()
        temp_score += int(temp)

        c.execute('SELECT score1 FROM gene WHERE BinaryEncoding =?', (genescore3,))
        temp = str(c.fetchone())
        temp = temp.replace(' ', '')[:-2].upper()
        temp = temp.replace(' ', '')[1:].upper()
        temp_score += int(temp)

        c.execute('SELECT score1 FROM gene WHERE BinaryEncoding =?', (genescore4,))
        temp = str(c.fetchone())
        temp = temp.replace(' ', '')[:-2].upper()
        temp = temp.replace(' ', '')[1:].upper()
        temp_score += int(temp)

        c.execute('SELECT score1 FROM gene WHERE BinaryEncoding =?', (genescore5,))
        temp = str(c.fetchone())
        temp = temp.replace(' ', '')[:-2].upper()
        temp = temp.replace(' ', '')[1:].upper()
        temp_score += int(temp)

        c.execute('SELECT score1 FROM gene WHERE BinaryEncoding =?', (genescore6,))
        temp = str(c.fetchone())
        temp = temp.replace(' ', '')[:-2].upper()
        temp = temp.replace(' ', '')[1:].upper()
        temp_score += int(temp)

        if switch1 == True:
            c.execute('SELECT score2 FROM gene WHERE BinaryEncoding =?', (genescore7,))
            temp = str(c.fetchone())
            temp = temp.replace(' ', '')[:-2].upper()
            temp = temp.replace(' ', '')[1:].upper()
            temp_score += int(temp)

        elif switch1 == False:
            c.execute('SELECT score1 FROM gene WHERE BinaryEncoding =?', (genescore7,))
            temp = str(c.fetchone())
            temp = temp.replace(' ', '')[:-2].upper()
            temp = temp.replace(' ', '')[1:].upper()
            temp_score += int(temp)

        else:
            sys.exit()

        if switch2 == True:
            c.execute('SELECT score2 FROM gene WHERE BinaryEncoding =?', (genescore8,))
            temp = str(c.fetchone())
            temp = temp.replace(' ', '')[:-2].upper()
            temp = temp.replace(' ', '')[1:].upper()
            temp_score += int(temp)

        elif switch2 == False:
            c.execute('SELECT score1 FROM gene WHERE BinaryEncoding =?', (genescore8,))
            temp = str(c.fetchone())
            temp = temp.replace(' ', '')[:-2].upper()
            temp = temp.replace(' ', '')[1:].upper()
            temp_score += int(temp)

        else:
            sys.exit()

        c.execute('SELECT score1 FROM gene WHERE BinaryEncoding =?', (genescore9,))
        temp = str(c.fetchone())
        temp = temp.replace(' ', '')[:-2].upper()
        temp = temp.replace(' ', '')[1:].upper()
        temp_score += int(temp)

        c.execute('SELECT score1 FROM gene WHERE BinaryEncoding =?', (genescore10,))
        temp = str(c.fetchone())
        temp = temp.replace(' ', '')[:-2].upper()
        temp = temp.replace(' ', '')[1:].upper()
        temp_score += int(temp)

        return temp_score

    #given a genome, obtain a list of genes present in an organism
    #add condiion so that essential genes are not counted twice
    def getGenes(self):
        gene_list= []
        c = conn.cursor()
        genename1 = self.genome[0:6]
        genename2 = self.genome[6:12]
        genename3 = self.genome[12:18]
        genename4 = self.genome[18:24]
        genename5 = self.genome[24:30]
        genename6 = self.genome[30:36]
        genename7 = self.genome[36:42]
        genename8 = self.genome[42:48]
        genename9 = self.genome[48:54]
        genename10 = self.genome[54:60]

        c.execute('SELECT GeneName FROM gene WHERE BinaryEncoding =?', (genename1,))
        temp = str(c.fetchone())
        temp = temp.replace(' ', '')[:-3].lower()
        temp = temp.replace(' ', '')[3:].lower()
        c.execute('SELECT Allele FROM gene WHERE BinaryEncoding =?', (genename1,))
        temp2 = str(c.fetchone())
        temp2 = temp2.replace(' ', '')[:-2].lower()
        temp2 = temp2.replace(' ', '')[1:].lower()
        temp_list1 = (temp + temp2)

        c.execute('SELECT GeneName FROM gene WHERE BinaryEncoding =?', (genename2,))
        temp = str(c.fetchone())
        temp = temp.replace(' ', '')[:-3].lower()
        temp = temp.replace(' ', '')[3:].lower()
        c.execute('SELECT Allele FROM gene WHERE BinaryEncoding =?', (genename2,))
        temp2 = str(c.fetchone())
        temp2 = temp2.replace(' ', '')[:-2].lower()
        temp2 = temp2.replace(' ', '')[1:].lower()
        temp_list2 = (temp + temp2)

        c.execute('SELECT GeneName FROM gene WHERE BinaryEncoding =?', (genename3,))
        temp = str(c.fetchone())
        temp = temp.replace(' ', '')[:-3].lower()
        temp = temp.replace(' ', '')[3:].lower()
        c.execute('SELECT Allele FROM gene WHERE BinaryEncoding =?', (genename3,))
        temp2 = str(c.fetchone())
        temp2 = temp2.replace(' ', '')[:-2].lower()
        temp2 = temp2.replace(' ', '')[1:].lower()
        temp_list3 = (temp + temp2)

        c.execute('SELECT GeneName FROM gene WHERE BinaryEncoding =?', (genename4,))
        temp = str(c.fetchone())
        temp = temp.replace(' ', '')[:-3].lower()
        temp = temp.replace(' ', '')[3:].lower()
        c.execute('SELECT Allele FROM gene WHERE BinaryEncoding =?', (genename4,))
        temp2 = str(c.fetchone())
        temp2 = temp2.replace(' ', '')[:-2].lower()
        temp2 = temp2.replace(' ', '')[1:].lower()
        temp_list4 = (temp + temp2)

        c.execute('SELECT GeneName FROM gene WHERE BinaryEncoding =?', (genename5,))
        temp = str(c.fetchone())
        temp = temp.replace(' ', '')[:-3].lower()
        temp = temp.replace(' ', '')[3:].lower()
        c.execute('SELECT Allele FROM gene WHERE BinaryEncoding =?', (genename5,))
        temp2 = str(c.fetchone())
        temp2 = temp2.replace(' ', '')[:-2].lower()
        temp2 = temp2.replace(' ', '')[1:].lower()
        temp_list5 = (temp + temp2)

        c.execute('SELECT GeneName FROM gene WHERE BinaryEncoding =?', (genename6,))
        temp = str(c.fetchone())
        temp = temp.replace(' ', '')[:-3].lower()
        temp = temp.replace(' ', '')[3:].lower()
        c.execute('SELECT Allele FROM gene WHERE BinaryEncoding =?', (genename6,))
        temp2 = str(c.fetchone())
        temp2 = temp2.replace(' ', '')[:-2].lower()
        temp2 = temp2.replace(' ', '')[1:].lower()
        temp_list6 = (temp + temp2)

        c.execute('SELECT GeneName FROM gene WHERE BinaryEncoding =?', (genename7,))
        temp = str(c.fetchone())
        temp = temp.replace(' ', '')[:-3].lower()
        temp = temp.replace(' ', '')[3:].lower()
        c.execute('SELECT Allele FROM gene WHERE BinaryEncoding =?', (genename7,))
        temp2 = str(c.fetchone())
        temp2 = temp2.replace(' ', '')[:-2].lower()
        temp2 = temp2.replace(' ', '')[1:].lower()
        temp_list7 = (temp + temp2)

        c.execute('SELECT GeneName FROM gene WHERE BinaryEncoding =?', (genename8,))
        temp = str(c.fetchone())
        temp = temp.replace(' ', '')[:-3].lower()
        temp = temp.replace(' ', '')[3:].lower()
        c.execute('SELECT Allele FROM gene WHERE BinaryEncoding =?', (genename8,))
        temp2 = str(c.fetchone())
        temp2 = temp2.replace(' ', '')[:-2].lower()
        temp2 = temp2.replace(' ', '')[1:].lower()
        temp_list8 = (temp + temp2)

        c.execute('SELECT GeneName FROM gene WHERE BinaryEncoding =?', (genename9,))
        temp = str(c.fetchone())
        temp = temp.replace(' ', '')[:-3].lower()
        temp = temp.replace(' ', '')[3:].lower()
        c.execute('SELECT Allele FROM gene WHERE BinaryEncoding =?', (genename9,))
        temp2 = str(c.fetchone())
        temp2 = temp2.replace(' ', '')[:-2].lower()
        temp2 = temp2.replace(' ', '')[1:].lower()
        temp_list9 = (temp + temp2)

        c.execute('SELECT GeneName FROM gene WHERE BinaryEncoding =?', (genename10,))
        temp = str(c.fetchone())
        temp = temp.replace(' ', '')[:-3].lower()
        temp = temp.replace(' ', '')[3:].lower()
        c.execute('SELECT Allele FROM gene WHERE BinaryEncoding =?', (genename10,))
        temp2 = str(c.fetchone())
        temp2 = temp2.replace(' ', '')[:-2].lower()
        temp2 = temp2.replace(' ', '')[1:].lower()
        temp_list10 = (temp + temp2)

        gene_list = [temp_list1, temp_list2, temp_list3, temp_list4, temp_list5, temp_list6, temp_list7, temp_list8, temp_list9, temp_list10]

        return gene_list

    #generate mutation throughout an organims genome (called during new generation)
    def generateMutation(self):
        temp_genome = ''
        j = 0
        while j < len(self.genome):
             temp_letter = self.genome[j]
             lucky_number = randint(1,500)
             #picks number between 1 and 60. If number = 36, flips 0 to 1 or vice versa.
             if(lucky_number == 45):
                if(self.genome[j] == '0'):
                    temp_genome += '1'
                else:
                    temp_genome += '0'
             else: temp_genome += self.genome[j]

             j += 1
        self.genome = ''
        self.genome = temp_genome

    #updates sum of populations fitness
    def updateTotalScore(self):
        #if indivdual score < 600, it lacks essential genes to survive.
        if (self.score < 600):
            tempp_score = 0
        else:
            tempp_score = self.score
        return tempp_score

def createFasta(inputList):
    file = open("phylogenetics.fasta", 'w')
    for s in range(100):
        s = inputList[s]
        file.write(">" + s.name + '\n')
        file.write(s.dnaGenome + '\n')
    file.close()
    #egfr_tree = Phylo.read("phylogenetics.fasta", "newick")
    #Phylo.draw_ascii(egfr_tree)
    align = AlignIO.read("phylogenetics.fasta", "fasta")
    print(align)
    cmdline = PhymlCommandline(input='phylogenetics.fasta', datatype='nt', model='WAG', alpha='e', bootstrap=100)
    out_log, err_log = cmdline()


#selects genes in current generation to reproduce for next generation
def selectNextGeneration(inputOrgList):
    #initialize list
    list = []
    #creates a list of size = sum of all organisms finess scores.
    for r in range(100):
        r = inputOrgList[r]
        k = 0
        #if indivdual score < 600, it lacks essential genes to survive.
        if(r.score > 600):
            #higher score = more entries in list = better probability of getting selected
            while k < r.score:
                list.append(r)
                k += 1
    complete_list = random.sample(list, 20)
    return complete_list

#uniform crossover ratio of 0.5 for 2 randomly selected organisms
def crossOver(inputOrgList, nextGenList):

    crossNextGen = []
    #parents need to be paired up
    while(len(nextGenList) > 0):
        entry1 = nextGenList.pop()
        entry2 = nextGenList.pop()

        child = 1
        for i in range(1,11):
            rand1 = random.randint(0,1)
            rand2 = random.randint(0,1)
            rand3 = random.randint(0,1)
            rand4 = random.randint(0,1)
            rand5 = random.randint(0,1)
            rand6 = random.randint(0,1)
            rand7 = random.randint(0,1)
            rand8 = random.randint(0,1)
            rand9 = random.randint(0,1)
            rand10 = random.randint(0,1)
            #0: gets gene from dad. 1: gets gene from mom.
            if(rand1 == 0):
                genename1 = entry1.genome[0:6]
            elif(rand1 == 1):
                genename1 = entry2.genome[0:6]
            if(rand2 == 0):
                genename2 = entry1.genome[6:12]
            elif(rand2 == 1):
                genename2 = entry2.genome[6:12]
            if(rand3 == 0):
                genename3 = entry1.genome[12:18]
            elif(rand3 == 1):
                genename3 = entry2.genome[12:18]
            if(rand4 == 0):
                genename4 = entry1.genome[18:24]
            elif(rand4 == 1):
                genename4 = entry2.genome[18:24]
            if(rand5 == 0):
                genename5 = entry1.genome[24:30]
            elif(rand5 == 1):
                genename5 = entry2.genome[24:30]
            if(rand6 == 0):
                genename6 = entry1.genome[30:36]
            elif(rand6 == 1):
                genename6 = entry2.genome[30:36]
            if(rand7 == 0):
                genename7 = entry1.genome[36:42]
            elif(rand7 == 1):
                genename7 = entry2.genome[36:42]
            if(rand8 == 0):
                genename8 = entry1.genome[42:48]
            elif(rand8 == 1):
                genename8 = entry2.genome[42:48]
            if(rand9 == 0):
                genename9 = entry1.genome[48:54]
            elif(rand9 == 1):
                genename9 = entry2.genome[48:54]
            if(rand10 == 0):
                genename10 = entry1.genome[54:60]
            elif(rand10 == 1):
                genename10 = entry2.genome[54:60]

            newSequence = (genename1 + genename2 + genename3 + genename4 + genename5 + genename6 + genename7 + genename8 + genename9 + genename10)

            crossNextGen.append(newSequence)

    return crossNextGen

#replaces the 100 previous organism's genomes with tose produced from crossover(). Also mutates.
def insertnewGenes(inputOrgList, crossNextGen):
    outputNewGene = []

    global_fitness_score = 0
    for r in range(100):
        p = r

        r = inputOrgList[r]
        p = crossNextGen[p]
        #deletes old genome and replaces it wit crossover value
        r.genome = ''
        r.genome = p
        #introduce mutations
        r.generateMutation()
        #gets new score and geneList values
        r.refreshGenes(r.genome)
        outputNewGene.append(r)
        #calculate total fitness score & append to end of list
        global_fitness_score += r.score
    outputNewGene.append(global_fitness_score)
    return outputNewGene

#lists organisms and their properties
def showOrganisms(newOrganisms, k):
    print('Generation ' + str(k) + ' complete!')
    print(newOrganisms[100])


def initialGeneration():
    global_fitness_score = 0
    #list of organism objects
    organismlist = []

    for i in range(100):
        i = Organism('Organism_' + str(i+1))
        global_fitness_score += i.updateTotalScore()
        organismlist.append(i)
        #print(i.name, i.genome, i.genelist, i.score)
    #total fitness level appended to the end of the list (position 100)
    organismlist.append(global_fitness_score)
    print(global_fitness_score)
    return organismlist

def newGeneration(inputOrgList, counter):
    global_fitness_score = 0
    #select organisms for reproduction
    nextGenList = selectNextGeneration(inputOrgList)
    #crossover and get new organisms
    newGeneList = crossOver(inputOrgList, nextGenList)
    #insert new genes into old organisms & mutate
    newOrganisms = insertnewGenes(inputOrgList, newGeneList)
    #refresh genes and get new scores
    showOrganisms(newOrganisms, counter)

    return newOrganisms



def main():
    #generate initial organism list
    initialOrgList = initialGeneration()
    createFasta(initialOrgList)
    test = selectNextGeneration(initialOrgList)
    for element in test:
        print element.name
    counter = 1
    iterate = newGeneration(initialOrgList, counter)
    temp = iterate
    counter = 2
    for i in range(5):
        temp2 = newGeneration(temp, counter)
        temp = temp2
        counter += 1
    for element in temp:
        try:
            print(element.dnaGenome)
        except:
            print(element)

if __name__ == "__main__": main()
