#!/bin/sh
# the next line restarts using wish \
exec cg source "$0" "$@"

#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#


# see where we are, add the lib dir in the original location to auto_path, load extension
set script [file join [pwd] [info script]]
set scriptname [file tail $script]
while 1 {
	if {[catch {set script [file join [file dir $script] [file readlink $script]]}]} break
}
set readsaberdir [file dir $script]

proc formaterror {} {
	puts stderr "readsaber_makerefdir ?options? refdir genomesequence.fasta ?transcripts.gtf?"
	puts stderr "  with options: -organelles, -nolowgenecutoff"
	exit 1
}

# generic help on no args
if {[llength $argv] < 1} {
	formaterror
}

proc readsaber_makerefdir {args} {
	set organelles {}
	set transcripts {}
	set ontindex 1
	set pacbioindex 0
	set nolowgenecutoff 200000
	set groupchromosomes {}
	set pos 0
	foreach {key value} $args {
		switch $key {
			-organelles {
				set organelles $value
				incr pos 2
			}
			-nolowgenecutoff {
				set nolowgenecutoff $value
				incr pos 2
			}
			-groupchromosomes {
				set groupchromosomes $value
				incr pos 2
			}
			-ontindex {
				set ontindex $value
				incr pos 2
			}
			-pacbioindex {
				set pacbioindex $value
				incr pos 2
			}
			default break
		}
	}
	set args [lrange $args $pos end]
	if {[llength $args] > 3} {
		if {[string index [lindex $args 0] 0] == "-"} {
			puts stderr "error calling readsaber_makerefdir: unknown option "[lindex $args 0]", must be one of: -organelles, -transcripts"
		} else {
			formaterror
		}
		exit 1
	}
	foreach {refdir genomefasta transcripts} {{} {} {}} break
	foreach {refdir genomefasta transcripts} $args break
	if {$refdir eq ""} {
		puts stderr "error: argument refdir is obligatory"
		exit 1
	}
	if {$genomefasta eq ""} {
		puts stderr "error: argument genomefasta is obligatory"
		exit 1
	}
	if {[file exists $refdir]} {
		puts stderr "error: target refdir $refdir already exists"
		exit 1
	}
	set refdir [file normalize $refdir]
	set tail [file tail $refdir]
	set build $tail
	file mkdir $refdir
	file mkdir $refdir/extra

	# make ifas
	set result $refdir/genome_$tail.ifas
	puts "Making genome $result"
	catch_exec cg fas2ifas $genomefasta $result
	catch_exec samtools faidx $result

	# groupchromosomes
	if {$groupchromosomes ne ""} {
		puts "Making genome $result.groupchromosomes"
		groupchromosomes $result $groupchromosomes
	}

	# fullgenome
	unset -nocomplain a
	set rfile $refdir/extra/reg_${build}_fullgenome.tsv
	puts "Making $rfile"
	set data [file_read $result.fai]
	set o [open $rfile.temp w]
	puts $o chromosome\tbegin\tend
	list_foreach {chromosome len} [lrange [split [string trim $data] \n] 0 end] {
		set a($chromosome) 1
		puts $o $chromosome\t0\t$len
	}
	close $o
	file rename -force -- $rfile.temp $rfile

	if {$organelles ne ""} {
		puts "Writing $refdir/extra/reg_${build}_organelles.tsv"
		foreach organelle $organelles {
			if {![info exists a($organelle)]} {
				puts stderr "error making organelles file: $organelle is not a chromosome in the given genome"
				exit 1
			}
		}
		file_write $refdir/extra/reg_${build}_organelles.tsv chromosome\n[join $organelles \n]\n
	}

	# for cram
	catch_exec cg fasta2cramref $result $result.forcram

	# sequencedgenome
	set target $refdir/extra/reg_${build}_sequencedgenome.tsv.zst
	puts "Making $target"
	catch_exec cg calcsequencedgenome --stack 1 $result | cg zst > $target.temp
	file rename -force -- $target.temp $target

	if {$ontindex} {
		# minimap2_splice index
		puts "Making minimap2 ont index (can take large amount of memory)"
		catch_exec cg refseq_minimap2 $result splice
		catch_exec cg refseq_minimap2 $result splicesens
	}

	if {$pacbioindex} {
		# minimap2_splice index
		puts "Making minimap2 pacbio index (can take large amount of memory)"
		catch_exec cg refseq_minimap2 $result splice:hq
	}

	set target $refdir/genome_${build}.dict
	catch_exec samtools dict -o $target.temp $result
	file rename $target.temp $target

	if {[file exists $transcripts]} {
		# converting transcripts
		set root [file root [file tail $transcripts]]
		regsub ^gene_ $root {} root
		set target $refdir/gene_${build}_$root.tsv
		puts "Making $target"
		set ext [file extension [gzroot $transcripts]]
		if {$ext eq ".gtf"} {
			file copy -force $transcripts [file root $target].gtf
			catch_exec cg gtf2tsv $transcripts $target
		} elseif {$ext in ".gff .ggf2 .gff3"} {
			file copy $transcripts [file root $target]$ext
			catch_exec cg gff2tsv $transcripts $target
		} elseif {$ext eq ".tsv"} {
			file copy $transcripts $target
		} else {
			puts stderr "format of transcripts file $transcripts not supported: must be one of: gtf, tsv, gff"
			exit 1
		}
	
		set target $refdir/extra/reg_${build}_nolowgene200k.tsv.zst
		puts "Making $target"
		catch_exec cg distrreg_nolowgene $refdir $nolowgenecutoff
	}

	puts "Made $refdir (can now be used as refdir for readsaber)"
}

readsaber_makerefdir {*}$argv
