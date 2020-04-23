SCRIPTS := ../scripts
CUTOFF ?= 0.05

clean:
	rm -f 18k.lib.fasta  18k.lib.err

# Note we do not include the full reverse primer as it has a restriction site
# we are searching for. Thus we include all but one base of the site (so that
# the last codon does not produce the site) and paste the rest of the constant
# sequence on later
18k.lib.fasta: 18491_remove.fasta 18491_flipped_zipped_final_pairs.pairs re-sites.txt
	python $(SCRIPTS)/lib-gen.py \
	    --verbose \
	    $< \
	    --pairs $(word 2, $^) \
	    --re-sites $(word 3, $^) \
	    --spacer 'gtcAGAAGAGCCTCGCAGTGac' \
	    --forward 'GCTGGAGGCGAGGTTA' \
	    --reverse 'ggagAGA' \
	    --revcomp-x \
	    --cutoff $(CUTOFF) \
	    2> $(@:.fasta=.err) \
	    | awk '{print $$0"GACCAGAGTACCGCCGT"}' \
	    > $@
