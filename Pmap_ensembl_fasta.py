import pandas as pd
import Bio
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os

"""
example of ENSEMBL FASTA
    >ENSP00000474133.1 pep:known chromosome:GRCh37:15:20211158:20211188:-1 gene:ENSG00000270783.1 transcript:ENST00000604950.1 gene_biotype:IG_D_gene transcript_biotype:IG_D_gene
"""

def fastaToDF(filename):
    # parse sequence fasta file
    identifiers = [seq_record.id for seq_record in SeqIO.parse(filename, "fasta")]
    descr = [seq_record.description for seq_record in SeqIO.parse(filename, "fasta")]
    lengths = [len(seq_record.seq) for seq_record in SeqIO.parse(filename, "fasta")]
    proSeq = [seq_record.seq for seq_record in SeqIO.parse(filename, "fasta")]
    ensp = []
    enspv = []
    enst = []
    enstv = []
    ensg = []
    ensgv = []
    proseq = []
    for row in descr:
        splitrow = row.split(" ")
        for val in splitrow:
            if 'ENSG' in val:
                gene = val.split(":")[1]
                stablegene = gene.split(".")[0]
                ensgv.append(gene)
                ensg.append(stablegene)
            if 'ENST' in val:
                tx = val.split(":")[1]
                stabletx = tx.split(".")[0]
                enstv.append(tx)
                enst.append(stabletx)
    for id in identifiers:
        splitID = id.split('.')
        stable = splitID[0]
        ensp.append(stable)
        enspv.append(id)
    for ps in proSeq:
        s = str(ps)
        proseq.append(s)
    s1 = pd.Series(enspv, name='ENSPv')
    s2 = pd.Series(ensp, name= 'ENSP')
    s3 = pd.Series(enstv, name='ENSTv')
    s4 = pd.Series(enst, name= 'ENST')
    s5 = pd.Series(ensgv, name='ENSGv')
    s6 = pd.Series(ensg, name= 'ENSG')
    s7 = pd.Series(lengths, name='Length')
    s8 = pd.Series(proseq, name='proSequence')
    series = [s1, s2, s3, s4, s5, s6, s7, s8]
    df = pd.concat(series, axis=1)
    return(df)



