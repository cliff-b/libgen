SCRIPTS := ../scripts

clean:
	rm -f R20.lib.csv

# Note we do not include the full reverse primer as it has a restriction site
# we are searching for. Thus we include all but one base of the site (so that
# the last codon does not produce the site) and paste the rest of the constant
# sequence on later
R20.lib.csv: R20.fasta ecoli_codon_usage.txt
	python $(SCRIPTS)/codon_opt_hb.py \
	$< \
	$^ \
	2> $(@:.fasta=.err) \
	> $@