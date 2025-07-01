#!/bin/sh
# the next line restarts using wish \
exec cg source "$0" "$@"

# see where we are, add the lib dir in the original location to auto_path, load extension
if {![info exists argv]} {
	# interactive, take pwd
	set script [pwd]/readsaber
} else {
	set script [file join [pwd] [info script]]
	while 1 {
		if {[catch {set script [file join [file dir $script] [file readlink $script]]}]} break
	}
}
set testdir [file dir $script]
set readsaberdir [file dir $testdir]

puts "Using as testdir: $testdir"

source $testdir/tools.tcl

set ::env(PATH) $readsaberdir:$readsaberdir/bin:$::env(PATH)

cd $testdir

test polyt {basic} {
	file_write tmp/test.fastq [deindent {
		@seq1
		GCCGTTTTTTTTGCCGT
		+
		IIIIIIIIIIIIIIIII
		@seq2
		GCCGTTTTTTGCCGT
		+
		IIIIIIIIIIIIIII
		@seq3
		GTCGTCTTGTTGTCTTTTCTTTTTTGC
		+
		IIIIIIIIIIIIIIIIIIIIIIIIIII
		@seq4
		GTCGTCTTGTTGTCTTTTCTTTTTTGCGGTTTTTTTTTTGT
		+
		IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
	}]\n
	exec polyt < tmp/test.fastq
} {rname	rstart	rend	size	strand	chromosome
seq1	4	12	8	+	polyT
seq3	9	25	16	+	polyT
seq4	9	25	16	+	polyT
seq4	29	41	12	+	polyT}

test polyt {tabs/space in readnames} {
	file_write tmp/test.fastq [deindent {
		@seq1	test1
		GCCGTTTTTTTTGCCGT
		+
		IIIIIIIIIIIIIIIII
		@seq2	test2
		GCCGTTTTTTGCCGT
		+
		IIIIIIIIIIIIIII
		@seq3 test3
		GTCGTCTTGTTGTCTTTTCTTTTTTGC
		+
		IIIIIIIIIIIIIIIIIIIIIIIIIII
	}]\n
	exec polyt < tmp/test.fastq
} {rname	rstart	rend	size	strand	chromosome
seq1	4	12	8	+	polyT
seq3	9	25	16	+	polyT}

test readsaber {basic} {
	test_cleantmp
	exec readsaber_makerefdir tmp/genome_test data/genome_test.fa 2> tmp/test1_makeref.log
	file copy ../annotations/annotations_ontr10x.fa tmp/annot.fa
	file copy data/test1.fastq tmp/test1.fastq
	exec readsaber \
		-keepintermediate 1 \
		-addsequences 1 \
		-refseq tmp/genome_test \
		tmp/annot.fa \
		tmp/readannot-test.tsv \
		tmp/test1.fastq 2> tmp/test1.log
	exec diff tmp/readannot-test.tsv data/test1-readannot-test.tsv
	exec diff tmp/readannot-test_summary.tsv data/test1-readannot-test_summary.tsv
} {}

testsummarize
