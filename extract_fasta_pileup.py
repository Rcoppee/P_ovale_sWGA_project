#!/usr/bin/env python

import sys, argparse

parser = argparse.ArgumentParser()
parser.add_argument('--file', '-f', type = str, help = "mpileup file", required = True)
parser.add_argument('--output', '-o', type = str, help = "nom de l'output", required = True) 
args = parser.parse_args()
fichier_entree = args.file
fichier_sortie = args.output

cpt = 1
tab={}

FOUT = open(fichier_sortie+".fasta", 'w')

t = 0
j = 0

with open(fichier_entree) as Fic:	
	for nb_line, line in enumerate(Fic):
		data = line.split("\t")
		chro = data[0]
		nucleotide = data[2].strip().split("\n")[0]
		pos = data[1].strip()

		nb_pos = nb_line

		if chro in tab: #if chromosome in table
			cpt = int(tab[chro][0])

			if (nucleotide == ""): #if the nucleotide is not present, then add '-'
				j += 1
				nucleotide = "-"
			if (nucleotide == "*"): #if the nucleotide is '*' (indel), update to '-'
				j += 1
				nucleotide = "-"	

			if (int(pos) != int(cpt)): #we check the presence of lacking positions
				pos_manquante = int(pos) - int(cpt)
				cpt2 = 1 
				while (int(cpt2) <= int(pos_manquante)): #attribute '-' to lacking positions
					sequence = sequence + "-"
					cpt2 += 1 
					t += 1 
				cpt = int(pos)
				
			sequence = sequence + nucleotide.upper()

			cpt += 1 
			tab[chro] = (cpt, pos, sequence)

		#if chromosome not in table
		else:
			if (pos == "1"): #if the starting position is 1
				sequence = "" #initialisation
				cpt = int(pos)+1 
				if (nucleotide == ""): #if the nucleotide is not present, then add '-'
					nucleotide = "-"
				sequence = sequence+nucleotide.upper() 
				tab[chro]=(cpt, pos, sequence) 
			else: #if the starting position is not 1, complete lacking positions
				cpt = int(pos)+1 
				cpt2 = 1 
				sequence= "" 
				while (int(cpt2) < int(pos)): 
					sequence = sequence + "-" 
					cpt2 += 1 
					t += 1
				if (nucleotide == ""): 
					nucleotide = "-" 
				sequence = sequence+nucleotide.upper() 
				tab[chro]=(cpt, pos, sequence)


taille_chro=0
for x in tab:
	FOUT.write(">"+str(x)+"\n") #writing header
	FOUT.write(tab[x][2]) #writing sequence
	FOUT.write("\n")
	taille_chro = taille_chro + int(tab[x][1])
print("The file contains: "+str(nb_pos+1)+" lines")
print("The length of all chromosomes is: "+str(taille_chro))
print("The size of the final sequence is: "+str(len(tab[x][2].strip("\n"))))
print("The number of lacking positions is: "+str(t))
print("The number of lacking nucleotides is: "+str(j))
FOUT.close()

