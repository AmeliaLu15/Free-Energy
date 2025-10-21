from seqfold import fold, dg, dg_cache, dot_bracket
import pandas

meta = pandas.read_csv("../data/Free energy masterfile.csv")
with open("../data/rCRS.fasta") as f:
    lines = f.read()
lines = lines.split("\n")
refseq = "".join(lines[1:len(lines)])

cellid = 22
fname = "../data/"+str(cellid)+".fasta"
with open(fname) as f:
        lines = f.read()
mtDNA_lines = lines.split("\n")
mtDNA = mtDNA_lines[1]


def alignedFreeEnergy(fname, break1, break2, windowsize):
    with open(fname) as f:
        lines = f.read()

    mtDNA_lines = lines.split("\n")
    mtDNA = mtDNA_lines[1]
    while "N" in mtDNA:
        index = mtDNA.index("N")
        # For some reason there is one N in the reference sequence
        if refseq[index] == "N":
            replace = "C"
        else:
            replace = refseq[index]
        mtDNA = mtDNA[:mtDNA.index("N")]+replace+mtDNA[mtDNA.index("N")+1:]
    
    seq1 = mtDNA[(break1-int(windowsize/2)):(break1+int(windowsize/2))]
    seq2 = mtDNA[break2-int(windowsize/2):(break2+int(windowsize))]
    return(dg(seq1+seq2, temp = 37.0))

output = []
for cellid, start, end in zip(meta.Title,meta.Start,meta.End):
    free_energy = alignedFreeEnergy("../data/"+str(cellid)+".fasta", start, end, 50)
    newrow = [cellid,start,end, free_energy]
    print(newrow)
    output.append(newrow)

df = pandas.DataFrame(output,columns=["CellID","Start","End","FreeEnergy"])
df.to_csv("../output/FreeEnergyResults50n.csv")





