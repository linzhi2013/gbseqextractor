#!/usr/bin/python3
"""
Copyright (c) 2017-2018 Guanliang Meng <mengguanliang@foxmail.com>.

This file is part of MitoZ.

MitoZ is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MitoZ is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MitoZ.  If not, see <http://www.gnu.org/licenses/>.

"""

import argparse
from Bio import SeqIO
import re
import os
import sys

def get_para():
	description = """
Extract any CDS or rNRA or tRNA DNA sequences of genes from Genbank file.

Seqid will be the value of '/gene=' or '/product=', if they both were not
present, the gene will not be output!

version 20201128:
    Now we can handle compounlocation (feature location with "join")!


Warning:
    Each sequence in the result files corresponds to ONE feature
in the GenBank file, I will NOT combine multiple CDS of the same gene into ONE!


Please cite:
Guanliang Meng, Yiyuan Li, Chentao Yang, Shanlin Liu,
MitoZ: a toolkit for animal mitochondrial genome assembly, annotation
and visualization, Nucleic Acids Research, https://doi.org/10.1093/nar/gkz173

	"""

	parser = argparse.ArgumentParser(
		description=description,
		formatter_class=argparse.RawDescriptionHelpFormatter)

	parser.add_argument("-f", required=True, metavar="<STR>", help="Genbank file")


	parser.add_argument("-prefix", metavar="<STR>", required=True,
		help="prefix of output file. required.")

	parser.add_argument("-seqPrefix", metavar="<STR>", default="",
		help="prefix of each seq id. default: None")

	parser.add_argument("-types", nargs="+", default="CDS",
		choices=["CDS", "rRNA", "tRNA", "wholeseq", "gene"],
		help="what kind of genes you want to extract? wholeseq for whole fasta seq.[%(default)s]")

	parser.add_argument("-gi", default=False, action='store_true',
		help='''use gi number as sequence ID instead of accession number when "
			gi number is present. (default: accession number)''')

	parser.add_argument("-p", default=False, action='store_true',
		help="output the position information on the ID line. Warning: the position on ID line is 0 left-most! [%(default)s]")

	parser.add_argument("-t", default=False, action='store_true',
		help="output the taxonomy lineage on ID line [%(default)s]")

	parser.add_argument("-s", default=False, action='store_true',
		help="output the species name on the ID line [%(default)s]")

	parser.add_argument("-l", default=False, action='store_true',
		help="output the seq length on the ID line [%(default)s]")

	parser.add_argument("-rv", default=True, action='store_true',
		help="reverse and complement the sequences if " + \
			"the gene is on minus strand. Always True! ")

	parser.add_argument("-F", default=False, action='store_true',
		help="only output full length genes,i.e., exclude the genes with '>' or '<' in their location [%(default)s]")


	if len(sys.argv) == 1:
		parser.print_help()
		sys.exit()
	else:
		args = parser.parse_args()

	return args


def main():
	args = get_para()

	fh_gene = fh_cds = fh_rrna = fh_trna = fh_wholeseq = ""
	if "gene" in args.types:
		fh_gene = open(args.prefix+".gene", 'w')

	if "CDS" in args.types:
		fh_cds = open(args.prefix+".cds", 'w')

	if "rRNA" in args.types:
		fh_rrna = open(args.prefix+".rrna", 'w')

	if "tRNA" in args.types:
		fh_trna = open(args.prefix+".trna", 'w')

	if "wholeseq" in args.types:
		fh_wholeseq = open(args.prefix + ".fasta", 'w')

	for rec in SeqIO.parse(args.f, 'gb'):
		ass_num = rec.id
		#print(rec.id)
		#print(rec.name)
		#print(rec.description)
		#print(rec.letter_annotations)
		#print(rec.annotations)

		# wholeseq
		if fh_wholeseq:
			if args.seqPrefix:
				wholeseq_idline = ">" + args.seqPrefix
			else:
				wholeseq_idline = ">"

			if args.gi:
				try:
					wholeseq_idline += rec.annotations['gi']
				except KeyError:
					wholeseq_idline += rec.id
			else:
				wholeseq_idline += rec.id

			if args.l:
				wholeseq_idline += ";len="+str(len(rec))

			if args.s:
				species = str(rec.annotations['organism'])
				wholeseq_idline += ";" + species

			if args.t:
				taxonomy = str(rec.annotations['taxonomy'])
				taxonomy = taxonomy.replace("'", "")
				wholeseq_idline += ";"+ taxonomy

			print(wholeseq_idline, file=fh_wholeseq)
			print(rec.seq, file=fh_wholeseq)

		# CDS, tRNA, rRNA, gene
		for fea in rec.features:
			if fea.type in args.types:
				if args.F:
					if '>' in str(fea.location)  or '<' in str(fea.location):
						continue

				start = fea.location.start
				end = fea.location.end

				gene = ""
				product = ""

				if 'gene' in fea.qualifiers:
					gene = fea.qualifiers['gene'][0]
				elif 'product' in fea.qualifiers:
					product = fea.qualifiers['product'][0]
					gene = product
				else:
					print(ass_num, "Warning: NO gene or product tag! this gene is not output!\n", fea.location, file=sys.stderr)
					continue

				if args.seqPrefix:
					idline = ">" + args.seqPrefix
				else:
					idline = ">"

				if args.gi:
					try:
						idline += "%s;%s" % (rec.annotations['gi'], gene)
					except KeyError:
						idline += "%s;%s" % (rec.id, gene)
				else:
					idline += "%s;%s" % (rec.id, gene)

				if args.l:
					idline += ";len=" + str(len(fea.location))

				if args.p:
					idline += ";" + str(fea.location)


				if args.s:
					species = str(rec.annotations['organism'])
					idline += ";" + species

				if args.t:
					taxonomy = str(rec.annotations['taxonomy'])
					taxonomy = taxonomy.replace("'", "")
					idline += ";"+ taxonomy

				# iter each part
				if len(fea.location.parts) > 1:
					gene_seq = ''
					for each_part in fea.location.parts:
						strand = each_part.strand
						each_part_seq = each_part.extract(rec.seq)
						if strand == -1:
							each_part_seq = each_part_seq.reverse_complement()
						gene_seq += each_part_seq

				else:
					gene_seq = fea.location.extract(rec.seq)
					strand = fea.strand
					if strand == -1:
						gene_seq = gene_seq.reverse_complement()


				if fea.type == "gene":
					print(idline, file=fh_gene)
					print(gene_seq, file=fh_gene)
				elif fea.type == "CDS":
					print(idline, file=fh_cds)
					print(gene_seq, file=fh_cds)
				elif fea.type == "rRNA":
					print(idline, file=fh_rrna)
					print(gene_seq, file=fh_rrna)
				elif fea.type == "tRNA":
					print(idline, file=fh_trna)
					print(gene_seq, file=fh_trna)

	if fh_gene:
		fh_gene.close()

	if fh_cds:
		fh_cds.close()

	if fh_rrna:
		fh_rrna.close()

	if fh_trna:
		fh_trna.close()

	if fh_wholeseq:
		fh_wholeseq.close()


if __name__ == '__main__':
	main()
