SCRIPTS := ../scripts
CUTOFF ?= 0.1

clean:
	rm -f all-8k.lib.fasta all-8k.lib.err

# Note we do not include the full reverse primer as it has a restriction site
# we are searching for. Thus we include all but one base of the site (so that
# the last codon does not produce the site) and paste the rest of the constant
# sequence on later
R8000.lib.fasta: all-8k.fasta all-8k.pairs re-sites.txt
	python $(SCRIPTS)/lib-gen.py \
	    --verbose \
	    $< \
	    --pairs $(word 2, $^) \
	    --re-sites $(word 3, $^) \
	    --spacer 'GTCAgaagagcagtgAC' \
	    --forward 'AAGAGGGACGCAGCATTA' \
	    --reverse 'CCAGAga' \
	    --revcomp-x \
	    --cutoff $(CUTOFF) \
	    2> $(@:.fasta=.err) \
	    | awk '{print $$0"GACCTTCAGCGCGACCAGCC"}' \
	    > $@
