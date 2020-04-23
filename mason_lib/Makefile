SCRIPTS := ../scripts
CUTOFF ?= 0.05

clean:
	rm -f mason.lib.fasta  mason.lib.err

# Note we do not include the full reverse primer as it has a restriction site
# we are searching for. Thus we include all but one base of the site (so that
# the last codon does not produce the site) and paste the rest of the constant
# sequence on later
mason.lib.fasta: mason.fasta re-sites.txt
	python $(SCRIPTS)/difflength-lib-gen.py \
	    --verbose \
	    $< \
	    --re-sites $(word 2, $^) \
	    --spacer 'gtcAGAAGAGCTTCTCGCAGTGac' \
	    --forward 'GAGTCTGAGCGGCGTTTA' \
	    --reverse 'ggagAGA' \
	    --yBsaI 'GACCGTCAGGCGAGCTAGGC' \
	    --revcomp-x \
	    --cutoff $(CUTOFF) \
	    --outerxfwdlist ACTGGTGCGTCGTCT AGCGAAACCGTGCGT ACCGGTTTCCACGCA GGGTTCGAGCGGGAG TGGGCGCCAAGAACC \
		--outerxrevlist GGCGAACACTTCCCG GGGACGCGGGAGAGT GTTCGCGCTTGGAGG CCTCTCTGCGCGCTA CTTGGGCGCAGTGGG \
		--outeryfwdlist ACTCGACGGCCTCTG CCCGGATCCCTTGCT GCTTCCGCCCGTAGA CGAGCCGTGGTTCCT ACGCCGAATCCCACA \
		--outeryrevlist CTTCAACGCGCGTGT GAGGCCAATGACGCG CCACCAGCGCGAATG AACACCCTCCCAGCG CCCGACGTTCCAGCA \
	    2> $(@:.fasta=.err) \
	    > $@

