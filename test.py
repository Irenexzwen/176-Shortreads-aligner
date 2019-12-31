from project import *
from shared import *
from evaluation import *

def import_transcriptome(file):
	genes = dict()
	isoforms = dict()
	exons = dict()

	unknown_genes = dict()
	unknown_isoforms = dict()
	unknown_exons = dict()

	with open(file) as f:
		lines = f.read().splitlines()

	i=0
	for line in lines:
		split = line.split()

		object_type = split[0]
		object_id = split[1]

		if object_type == "gene":
			genes[object_id] = Gene(object_id, [])

			isoform_ids = split[2].split(";")

			for isoform_id in isoform_ids:
				isoforms[isoform_id] = Isoform(isoform_id, [])

				genes[object_id].isoforms.append(isoforms[isoform_id])

		elif object_type == "isoform":
			exon_ids = split[2].split(";")

			for exon_id in exon_ids:
				i+=1
				exons[exon_id] = Exon(exon_id, 0, 0)

				isoforms[object_id].exons.append(exons[exon_id])

		elif object_type == "exon":
			exons[object_id].start = int(split[2])
			exons[object_id].end = int(split[3])

		elif object_type == "unknown_gene":
			unknown_genes[object_id] = Gene(object_id, [])

			isoform_ids = split[2].split(";")

			for isoform_id in isoform_ids:
				unknown_isoforms[isoform_id] = Isoform(isoform_id, [])

				unknown_genes[object_id].isoforms.append(unknown_isoforms[isoform_id])

		elif object_type == "isoform":
			exon_ids = split[2].split(";")

			for exon_id in exon_ids:
				i+=1
				unknown_exons[exon_id] = Exon(exon_id, 0, 0)

				unknown_isoforms[object_id].exons.append(unknown_exons[exon_id])

		elif object_type == "exon":
			unknown_exons[object_id].start = int(split[2])
			unknown_exons[object_id].end = int(split[3])

	return genes.values(), unknown_genes.values()

def import_genome(file):
	with open(file) as f:
		lines = f.read().splitlines()

	return lines[1]

def import_reads(file):
	with open(file) as f:
		lines = f.read().splitlines()

	return dict([(lines[2*i], lines[2*i+1]) for i in range(len(lines)//2)])