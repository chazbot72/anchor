db_root = "/mnt/c/Users/bowman/Desktop/anchor/"
ali_root = "/mnt/c/Users/bowman/Desktop/anchor/alignments/"
ali_infile = "/mnt/c/Users/bowman/Desktop/anchor/inFile.fasta"
hxb2 = db_root + "HXB2.fasta"
db_content = db_root + "HIV1_SFL_2018_env_PRO.fasta"
aa_code = ["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y","-"]
line_break = "--------------------------------------------------------------------------------"

align = False
analyze = False
query = True

import sys, os, datetime
from string import ascii_lowercase as lowercase
import pickle

import re

_nsre = re.compile('([0-9]+)')
def natural_sort_key(s):
    return [int(text) if text.isdigit() else text.lower()
            for text in re.split(_nsre, s)]

from Bio import SeqIO
from Bio import AlignIO
from Bio.Align.Applications import ClustalOmegaCommandline

def print_progress(iteration, total, prefix = '', suffix = '', decimals = 2, barLength = 100):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : number of decimals in percent complete (Int) 
        barLength   - Optional  : character length of bar (Int) 
    """
    filledLength    = int(round(barLength * iteration / float(total)))
    percents        = round(100.00 * (iteration / float(total)), decimals)
    bar             = '#' * filledLength + '-' * (barLength - filledLength)
    sys.stdout.write('%s [%s] %s%s %s\r' % (prefix, bar, percents, '%', suffix)),
    sys.stdout.flush()
    if iteration == total:
        print("\n")

if align:
    hxb2_seq =  next(SeqIO.parse(hxb2, "fasta"))
    print(line_break)
    print("HXB2 sequence loaded (%s, %sbp)" % (hxb2_seq.id, len(hxb2_seq.seq)))

    db_content_seqs = {}
    records = SeqIO.parse(db_content, "fasta")
    for record in records:
    	db_content_seqs[record.id] = record

    print("Sequence database loaded (%s records)" % len(db_content_seqs))
    print(line_break)

    print("Aligning Sequences to HXB2 Reference")
    total = len(db_content_seqs)
    count = 0
    for record in db_content_seqs.keys():
    	print_progress(count, total, prefix = 'Aligning...', suffix = '', barLength = 50)
    	records = [hxb2_seq, db_content_seqs[record]]
    	ali_outfile = ali_root + db_content_seqs[record].id
    	
    	with open(ali_infile, 'w') as handle:
    		SeqIO.write(records, handle, "fasta")

    	clustalomega_cline = ClustalOmegaCommandline(infile=ali_infile, outfile=ali_outfile, verbose=True, auto=True, force=True, threads=16)
    	clustalomega_cline()
    	count += 1

if analyze:

    L1 = {} #Position Keyed Dictionary

    for alignment in os.listdir(ali_root):

        ali_raw = AlignIO.read(open(ali_root + alignment),"fasta")
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
                L1[str(hxb2_count)][hxb2_res].add("HXB2")
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
                L1[str(str(hxb2_count) + alpha)][hxb2_res].add("HXB2")
                if sub_res not in L1[str(str(hxb2_count) + alpha)]:
                    L1[str(str(hxb2_count) + alpha)][sub_res] = set()
                L1[str(str(hxb2_count) + alpha)][sub_res].add(alignment)

                real_count += 1
                hxb2_alpha += 1

    pickle.dump(L1, open('db.pickle', 'wb'))

if query:

    if len(sys.argv) < 2:
        print("You must at minimum define anchor residues as comma separated k:v pairs EX: \"240:T,232:T/K\"")
        print("You can set positions to report (based on hxb2 seq) as comma separated values EX: \"240, 241, 242\"")
        sys.exit()

    anchor_res = sys.argv[1].split(',')
    try:
        report_res = [x.strip() for x in sys.argv[2].split(',')]
    except:
        report_res = False

    print(line_break)
    print("Performing Anchor Database Search...")
    print("Anchor Points: %s" % anchor_res)
    if report_res:
        print("Report Points: %s" % report_res)
    else:
        print("Report Points: ALL RES")
    print(line_break)
    print("loading db...")
    L1 = pickle.load(open('db.pickle', 'rb'))

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
    report_res.sort(key=natural_sort_key)

    outTime = datetime.datetime.now().strftime("%Y%m%d-%H%M")

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

    with open("result_%s_all.csv" % outTime, 'w') as handle:
        handle.write("POS," + ','.join(aa_code) + '\n')
        for res in report_res:
            write_string = res
            for aa in aa_code:
                if aa in L1[res].keys():
                    write_string += ",%s" % str(len(L1[res][aa]))
                else:
                    write_string += ",0"
            handle.write(write_string + '\n')





