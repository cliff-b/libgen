SCRIPTS := ../scripts

clean:
	rm -f R20.lib.csv

# smaall make file for relatively unrobust library generator
R20.lib.csv: R20.fasta ecoli_codon_usage.txt
	python $(SCRIPTS)/codon_opt_hb.py \
	$< \
	$(word 2, $^) \
	2> $(@:.csv=.err) \
	> $@