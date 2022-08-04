


import sys, os, datetime
from string import ascii_lowercase as lowercase
import pickle
import argparse
import re

from Bio import SeqIO
from Bio import AlignIO
from Bio.Align.Applications import ClustalOmegaCommandline
from tqdm import tqdm

aa_code = ["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y","-"]

_nsre = re.compile('([0-9]+)')
def natural_sort_key(s):
	return [int(text) if text.isdigit() else text.lower() for text in re.split(_nsre, s)]

## Arguments
parser = argparse.ArgumentParser()
mainGroup = parser.add_argument_group('Main Arguments')
mainGroup.add_argument('--action', dest='action', choices=['align','analyze','query'], required=True, help="Action to Perform")
mainGroup.add_argument('--anchor', type=str, help="Anchor sequence")
mainGroup.add_argument('--ali_dbfile', type=str, help="Sequence database (FASTA format)")
mainGroup.add_argument('--query', type=str, nargs='?', help="Query string")
mainGroup.add_argument('--report', type=str, nargs='?', help="Report string")
mainGroup.add_argument('--o', type=str, help="Output filename for query (omit extension)")


advancedGroup = parser.add_argument_group('Advanced Arguments')
advancedGroup.add_argument('--db_root', default=os.getcwd(), type=str, help="DB root (default is $PWD)")
advancedGroup.add_argument('--db_pickle', type=str, help="DB pickle file name (default is $PWD/db.pickle)")
advancedGroup.add_argument('--ali_root', type=str, help="Alignment root (default is <db_root>/alignments)")
advancedGroup.add_argument('--ali_infile', type=str, help="Temp file for ClustalO alignment (default is <db_root>/inFile.fasta)")


## Do pairwise alignments
def align_db(hxb2, ali_root, db_content):

	if not os.path.isdir(ali_root):
		os.mkdir(ali_root)

	hxb2_seq =  next(SeqIO.parse(hxb2, "fasta"))
	print("HXB2 sequence loaded (%s, %sbp)" % (hxb2_seq.id, len(hxb2_seq.seq)))

	db_content_seqs = {}
	records = SeqIO.parse(db_content, "fasta")
	for record in records:
		db_content_seqs[record.id] = record

	print("Sequence database loaded (%s records)" % len(db_content_seqs))

	print("Aligning Sequences to HXB2 Reference")
	total = len(db_content_seqs)
	count = 0
	for record in tqdm(db_content_seqs.keys()):
		records = [hxb2_seq, db_content_seqs[record]]
		ali_outfile = ali_root + db_content_seqs[record].id
		
		with open(ali_infile, 'w') as handle:
			SeqIO.write(records, handle, "fasta")

		clustalomega_cline = ClustalOmegaCommandline(infile=ali_infile, outfile=ali_outfile, verbose=True, auto=True, force=True, threads=16)
		clustalomega_cline()
		count += 1

	print("Done, pairwise alignments written out to %s" % ali_root)

## Analyze the alignments and build the db.
def analyze_db(ali_root, db_pickle):

	L1 = {} #Position Keyed Dictionary

	print("Analyzing Alignments")
	for alignment in tqdm(os.listdir(ali_root)):

		ali_raw = AlignIO.read(open(os.path.join(ali_root + alignment)),"fasta")
		hxb2_seq = ali_raw[0].seq
		sub_seq = ali_raw[1].seq

		hxb2_count = 1
		hxb2_alpha = 0
		real_count = 0

		for hxb2_res, sub_res in zip(hxb2_seq, sub_seq):
			if real_count == 0 and hxb2_res == '-':
				real_count += 1
				continue
			if hxb2_res != '-':
				if str(hxb2_count) not in L1:
					L1[str(hxb2_count)] = {}
				if hxb2_res not in L1[str(hxb2_count)]:
					L1[str(hxb2_count)][hxb2_res] = set()
				L1[str(hxb2_count)][hxb2_res].add("Anchor")
				if sub_res not in L1[str(hxb2_count)]:
					L1[str(hxb2_count)][sub_res] = set()
				L1[str(hxb2_count)][sub_res].add(alignment)

				real_count += 1
				hxb2_count += 1
				hxb2_alpha = 0
			else:
				if hxb2_alpha > 25:
					alpha = lowercase[25] + lowercase[hxb2_alpha-26]
				else:
					alpha = lowercase[hxb2_alpha]

				if str(str(hxb2_count) + alpha) not in L1:
					L1[str(str(hxb2_count) + alpha)] = {}
				if hxb2_res not in L1[str(str(hxb2_count) + alpha)]:
					L1[str(str(hxb2_count) + alpha)][hxb2_res] = set()
				L1[str(str(hxb2_count) + alpha)][hxb2_res].add("Anchor")
				if sub_res not in L1[str(str(hxb2_count) + alpha)]:
					L1[str(str(hxb2_count) + alpha)][sub_res] = set()
				L1[str(str(hxb2_count) + alpha)][sub_res].add(alignment)

				real_count += 1
				hxb2_alpha += 1

	pickle.dump(L1, open(db_pickle, 'wb'))
	print("Done, written out to %s" % db_pickle)

def query_db(anchor_res, report_res, db_pickle):

	if len(sys.argv) < 2:
		
		sys.exit()

	print("Performing Anchor Database Search...")
	print("Anchor Points: %s" % anchor_res)
	if report_res:
		print("Report Points: %s" % report_res)
	else:
		print("Report Points: ALL RES")
	print("loading db...")
	L1 = pickle.load(open(db_pickle, 'rb'))

	anchor_chain = 'init'
	for anchor in anchor_res:
		anchor = [x.strip() for x in anchor.split(":")]
		if "/" in anchor[1]:
			if anchor_chain == 'init':
				anchor_chain = set()
				for res in anchor[1].split('/'):
					anchor_chain = anchor_chain | L1[anchor[0]][res]
			else:
				temp_chain = set()
				for res in anchor[1].split('/'):
					temp_chain = temp_chain | (anchor_chain & L1[anchor[0]][res])
				anchor_chain = anchor_chain & temp_chain
		else:
			if anchor_chain == 'init':
				anchor_chain = L1[anchor[0]][anchor[1]]
			else:
				anchor_chain = anchor_chain & L1[anchor[0]][anchor[1]]
		print("queried %s, %s total entries match" % (anchor, len(anchor_chain)))
	print("writing report...")
	if not report_res:
		report_res = L1.keys()
	report_res = sorted(report_res, key=natural_sort_key)

	with open("result_%s.csv" % outTime, 'w') as handle:
		handle.write("POS," + ','.join(aa_code) + '\n')
		for res in report_res:
			write_string = res
			for aa in aa_code:
				if aa in L1[res].keys():
					write_string += ",%s" % str(len(anchor_chain & L1[res][aa]))
				else:
					write_string += ",0"
			handle.write(write_string + '\n')

if __name__ == '__main__':
	args = parser.parse_args()

	action = args.action
	db_root = args.db_root
	anchor = args.anchor
	db_content = args.ali_dbfile

	if args.o:
		_outfile = args.o + ".csv"
		_outfileall = args.o + "_all.csv"
	else:
		outTime = datetime.datetime.now().strftime("%Y%m%d-%H%M")
		_outfile = "result_%s.csv" % outTime
		_outfileall = "result_%s_all.csv" % outTime

	if args.ali_root:
		ali_root = args.ali_root
	else:
		ali_root = os.path.join(db_root, "alignments/")

	if args.ali_infile:
		ali_infile = args.ali_infile
	else:
		ali_infile = os.path.join(db_root, "inFile.fasta")

	if args.db_pickle:
		db_pickle = args.db_pickle
	else:
		db_pickle = os.path.join(db_root, "db.pickle")

	if action == 'align':
		align_db(anchor, ali_root, db_content)
	elif action == 'analyze':
		analyze_db(ali_root, db_pickle)
	elif action == 'query':

		try:
			if args.query:
				anchor_res = args.query.split(',')
				if args.report:
					report_res = [x.strip() for x in args.report.split(',')]
				else:
					report_res = False
			else:
				raise
		except:
			print("You must at minimum define anchor residues as comma separated k:v pairs EX: --query \"240:T,232:T/K\"")
			print("You can set positions to report (based on hxb2 seq) as comma separated values EX: --report \"240, 241, 242\"")
			sys.exit()
		
		query_db(anchor_res, report_res, db_pickle)

	print(args.query)
