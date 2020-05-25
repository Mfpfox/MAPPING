

import pandas as pd
import Bio
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os

os.chdir("/Users/mariapalafox/Box Sync/CODE_DATA/dir_MAPpaper/TSV_UNIPROT_xref/")

"""
from processUniprotFasta.py

0    >ENSP00000488240.1
1    pep
2    chromosome:GRCh38:CHR_HSCHR7_2_CTG6:142847306:142847317:1
3    gene:ENSG00000282253.1
4    transcript:ENST00000631435.1
5    gene_biotype:TR_D_gene
6    transcript_biotype:TR_D_gene gene_symbol:TRBD1 description:T cell 
        receptor beta
7    diversity
8    1
9    [Source:HGNC Symbol;Acc:HGNC:12158]
    GGGGGGGGGGGGG


    >ENSP00000474133.1 pep:known chromosome:GRCh37:15:20211158:20211188:-1 gene:ENSG00000270783.1 transcript:ENST00000604950.1 gene_biotype:IG_D_gene transcript_biotype:IG_D_gene
XYYDFWTGYYT

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

    #converting lists to pandas Series
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


def main():
    pep85 = "v85Homo_sapiens.GRCh37.pep.all.fa"
    pep92 = "v92Homo_sapiens.GRCh38.pep.all.fa"
    pep94 = "v94Homo_sapiens.GRCh38.pep.all.fa"
    pep96 = "v96Homo_sapiens.GRCh38.pep.all.fa"
    pep97 = "v97Homo_sapiens.GRCh38.pep.all.fa"
    v85 = fastaToDF(pep85)
    print("v85 shape: ", v85.shape)
    v92 = fastaToDF(pep92)
    print("v92 shape: ", v92.shape)
    v94 = fastaToDF(pep94)
    print("v94 shape: ", v94.shape)
    v96 = fastaToDF(pep96)
    print("v96 shape: ", v96.shape)
    v97 = fastaToDF(pep97)
    print("v97 shape: ", v97.shape)
    v85.to_csv("v85Homo_sapiens.GRCh37.pep.all.csv", index=False)
    v92.to_csv("v92Homo_sapiens.GRCh38.pep.all.csv", index=False)
    v94.to_csv("v94Homo_sapiens.GRCh38.pep.all.csv", index=False)
    v96.to_csv("v96Homo_sapiens.GRCh38.pep.all.csv", index=False)
    v97.to_csv("v97Homo_sapiens.GRCh38.pep.all.csv", index=False)

    # v85 shape:  (104763, 8)
    # v92 shape:  (107844, 8)
    # v94 shape:  (109095, 8)
    # v96 shape:  (109914, 8)
    # v97 shape:  (110048, 8)

    # test = "fastaTest.fa"
    # vtest = fastaToDF(test)
    # print(vtest)

main()



