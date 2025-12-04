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
		GTCGTCTTGTTGTCTTTTCTTTTTTTGC
		+
		IIIIIIIIIIIIIIIIIIIIIIIIIII
		@seq4
		GTCGTCTTGTTGTCTTTTCTTTTTTTGCGGTTTTTTTTTTGT
		+
		IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
	}]\n
	exec polyt 7 < tmp/test.fastq
} {rname	rstart	rend	size	numT	strand	chromosome
@seq1	4	12	8	8	+	polyT
@seq3	9	26	17	14	+	polyT
@seq4	9	26	17	14	+	polyT
@seq4	30	42	12	11	+	polyT}

test polyt {minT 11} {
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
		GTCGTCTTGTTGTCTTTTCTTTTTTTGC
		+
		IIIIIIIIIIIIIIIIIIIIIIIIIII
		@seq4
		GTCGTCTTGTTGTCTTTTCTTTTTTTGCGGTTTTTTTTTTGT
		+
		IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
	}]\n
	exec polyt 11 < tmp/test.fastq
} {rname	rstart	rend	size	numT	strand	chromosome
@seq3	9	26	17	14	+	polyT
@seq4	9	26	17	14	+	polyT
@seq4	30	42	12	11	+	polyT}

test readsaber {basic} {
	test_cleantmp
	exec readsaber_makerefdir tmp/genome_test data/genome_test.fa 2> tmp/test1_makeref.log
	file copy ../annotations/annotations_ontr10x.fa tmp/annot.fa
	file copy data/test1.fastq tmp/test1.fastq
	exec readsaber -stack 1 \
		-keepintermediate 1 \
		-addsequences 1 \
		-refseq tmp/genome_test \
		tmp/annot.fa \
		tmp/readannot-test.tsv \
		tmp/test1.fastq \
		2> tmp/test1.log
	exec diff tmp/readannot-test.tsv data/test1-readannot-test.tsv
	exec diff tmp/readannot-test_summary.tsv data/test1-readannot-test_summary.tsv
	exec diff tmp/readannot-test_shortsummary.tsv data/test1-readannot-test_shortsummary.tsv
} {}

test readsaber {multiple refseq} {
	test_cleantmp
	exec readsaber_makerefdir tmp/genome_test data/genome_test.fa 2> tmp/test1_makeref.log
	exec readsaber_makerefdir tmp/genome_test2 data/genome_test2.fa 2> tmp/test2_makeref.log
	file copy ../annotations/annotations_ontr10x.fa tmp/annot.fa
	file copy data/test1.fastq tmp/test1.fastq
	exec readsaber \
		-keepintermediate 1 \
		-addsequences 1 \
		-refseq tmp/genome_test \
		-refseq tmp/genome_test2 \
		tmp/annot.fa \
		tmp/readannot-test.tsv \
		tmp/test1.fastq 2> tmp/test1.log
	exec diff tmp/readannot-test.tsv data/test2-readannot-test.tsv
	exec diff tmp/readannot-test_summary.tsv data/test2-readannot-test_summary.tsv
} {}

test readsaber {-polyT 12} {
	test_cleantmp
	exec readsaber_makerefdir tmp/genome_test data/genome_test.fa 2> tmp/test1_makeref.log
	file copy ../annotations/annotations_ontr10x.fa tmp/annot.fa
	file copy data/test1.fastq tmp/test1.fastq
	exec readsaber -stack 1 -v 2 \
		-polyT 12 \
		-keepintermediate 1 \
		-addsequences 1 \
		-refseq tmp/genome_test \
		tmp/annot.fa \
		tmp/readannot-test.tsv \
		tmp/test1.fastq 2> tmp/test1.log
	exec diff tmp/readannot-test.tsv data/test-polyT12-readannot-test.tsv
	exec diff tmp/readannot-test_summary.tsv data/test-polyT12-readannot-test_summary.tsv
} {}

test readsaber {-polyT 0} {
	test_cleantmp
	exec readsaber_makerefdir tmp/genome_test data/genome_test.fa 2> tmp/test1_makeref.log
	file copy ../annotations/annotations_ontr10x.fa tmp/annot.fa
	file copy data/test1.fastq tmp/test1.fastq
	exec readsaber -stack 1 -v 2 \
		-polyT 0 \
		-keepintermediate 1 \
		-addsequences 1 \
		-refseq tmp/genome_test \
		tmp/annot.fa \
		tmp/readannot-test.tsv \
		tmp/test1.fastq 2> tmp/test1.log
	exec diff tmp/readannot-test.tsv data/test-polyT0-readannot-test.tsv
	exec diff tmp/readannot-test_summary.tsv data/test-polyT0-readannot-test_summary.tsv
} {}

test readsaber {-completeness 100} {
	test_cleantmp
	exec readsaber_makerefdir tmp/genome_test data/genome_test.fa 2> tmp/test1_makeref.log
	file copy ../annotations/annotations_ontr10x.fa tmp/annot.fa
	file copy data/test1.fastq tmp/test1.fastq
	exec readsaber -stack 1 -v 2 \
		-polyT 0 \
		-completeness 100 \
		-keepintermediate 1 \
		-addsequences 1 \
		-refseq tmp/genome_test \
		tmp/annot.fa \
		tmp/readannot-test.tsv \
		tmp/test1.fastq 2> tmp/test1.log
	exec diff tmp/readannot-test.tsv data/test-completeness100-readannot-test.tsv
	exec diff tmp/readannot-test_summary.tsv data/test-completeness100-readannot-test_summary.tsv
} {}

testsummarize
