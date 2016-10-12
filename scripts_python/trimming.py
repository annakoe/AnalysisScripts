from Bio import SeqIO
from Bio import Seq

def trim_toend(records, min_len):
    """Trims perfect adaptor sequences, checks read length.

    This is a generator function, the records argument should
    be a list or iterator returning SeqRecord objects.
    """
    for record in records:
        len_record = len(record) #cache this for later
        if len(record) < min_len:
            yield record
        elif len_record  >= min_len:
            #after trimming this will still be long enough
            yield (record[:-15])

original_reads = SeqIO.parse("inputfile.fa", "fasta")
trimmed_reads = trim_toend(original_reads, 15)
count = SeqIO.write(trimmed_reads, "trimmed.fa", "fasta")
print "Saved %i reads from fastq" % count
