#!/bin/sh
# the next line restarts using tclsh \
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

set ::env(PATH) $readsaberdir/bin:$::env(PATH)

# replace help_get from genomecomb so `readsaber -h` works
auto_load help_get
proc help_get {action} {
	set help [file_read $::readsaberdir/help/$action.wiki]
	return $help
}

# main proc, using job system
proc readsaber_job {args} {
	upvar job_logdir job_logdir
	set cmdline [clean_cmdline cg readsaber {*}$args]
	set refseqs {}
	set refseqannot {transcript$postfix}
	set keepintermediate 0
	set addsequences 0
	set alimethod minimap2_ontshort
	set refalimethod minimap2_splicesens
	set ignoreN 3
	set polyT 7
	set completeness 50
	set threads 2
	if {[lindex $args] in "-h -help help"} {
		help readsaber
		exit 0
	}
	cg_options readsaber args {
		-refseq {lappend refseqs [refseq $value]}
		-refseqannot {set refseqannot $value}
		-alimethod {set alimethod $value}
		-refalimethod {set refalimethod $value}
		-keepintermediate {set keepintermediate $value}
		-addsequences {set addsequences $value}
		-ignoreN {set ignoreN $value}
		-polyT {
			if {![isint $value]} {error "$value is not an integer; -polyT value must be an integer"}
			set polyT $value
		}
		-completeness {set completeness $value}
		-threads {set threads $value}
		-version - -V {
			puts "v0.1.0"
			exit 0
		}
	} {annotationfile result fastq} 3 ... {
		test
	}
	set annotationfile [file_absolute $annotationfile]
	set result [file_absolute $result]
	set fastq [file_absolute $fastq]
	set fastqs [list [file_absolute $fastq]]
	foreach fastq $args {
		lappend fastqs [file_absolute $fastq]
	}
	# index annotationfile
#	job readannot_index-[file tail $annotationfile] -deps {
#		$annotationfile
#	} -targets {
#		$annotationfile.map-ont
#	} -vars {
#		annotationfile
#	} -code {
#		exec samtools faidx $annotationfile
#		catch_exec minimap2 -x map-ont -k 5 -w 1 -d $annotationfile.map-ont $annotationfile
#	}
	set root [file root [gzroot [file tail $result]]]
	set resultdir [file dir $result]
	job_logfile [file dir $result]/sc_readsaber-[file tail $result] [file dir $result] $cmdline \
		{*}[versions genomecomb samtools dbdir zstd os]

	set method $alimethod ; set alipreset {}
	regexp {(^[^_]+)_(.*)$} $alimethod temp alimethod alipreset
	refseq_${alimethod}_job $annotationfile $alipreset

	# mapping
	set workdir [gzroot $result].temp
	if {[file exists $workdir] && ![file isdir $workdir]} {
		file delete $workdir
	}
	shadow_mkdir $workdir
	if {!$keepintermediate} {
		job_cleanup_add_shadow $workdir
	}
	set alist {}
	foreach fastq $fastqs {
		set tail [file root [gzroot [file tail $fastq]]]
		set tempfastq $workdir/$root-$tail.fastq.gz
		job maketempfastq-$root-$tail -deps {
			$fastq
		} -targets {
			$tempfastq
		} -vars {
			fastq tempfastq
		} -code {
			if {[file extension $fastq] eq ".bam"} {
				catch_exec samtools fastq -T "RG,CB,QT,MI,MM,ML,Mm,Ml" $fastq | cg bgzip > $tempfastq.temp
			} elseif {[file extension $fastq] eq ".gz"} {
				mklink $fastq $tempfastq.temp
			} else {
				catch_exec cg bgzip -o $tempfastq.temp $fastq
			}
			file rename $tempfastq.temp $tempfastq
		}
		set toadd [list $workdir/$root-$tail.annot.ali.tsv.zst]
		map_job -nohardclips 1 -paired 0 -threads $threads \
			-method $alimethod -preset $alipreset \
			$workdir/$root-$tail.annot.ali.tsv.sam.zst \
			$annotationfile \
			annotation \
			$tempfastq
		job readannot_match-$root-$tail -cores $threads -deps {
			$workdir/$root-$tail.annot.ali.tsv.sam.zst
		} -targets {
			$workdir/$root-$tail.annot.ali.tsv.zst
		} -vars {
			fastq annotationfile workdir root tail tempfastq threads
		} -code {
			cg sam2tsv -f {AS ms cs} $workdir/$root-$tail.annot.ali.tsv.sam.zst | cg select -f {
				rname="@$qname"
				{rstart=if($strand eq "+",$qstart,$seqlen - $qend)}
				{rend=if($strand eq "+",$qend,$seqlen - $qstart)}
				{size=$qend - $qstart}
				strand chromosome AS ms begin end mapquality cigar *
			} | cg select -s {rname rstart rend} | cg zst > $workdir/$root-$tail.annot.ali.tsv.temp.zst
			file rename -force $workdir/$root-$tail.annot.ali.tsv.temp.zst $workdir/$root-$tail.annot.ali.tsv.zst
		}
		set refnr 0
		foreach refseq $refseqs {
			incr refnr
			if {$refnr == 1} {set postfix ""} else {set postfix $refnr}
			set target $workdir/$root-$tail.refseq$postfix.ali.tsv.zst
			lappend toadd $target
			if {[file extension $refseq] in ".cram .bam .sam"} {
				set useali $refseq
			} else {
				set useali $workdir/$root-$tail.refseq$postfix.ali.sam.zst
				set method $refalimethod ; set refalipreset {}
				regexp {(^[^_]+)_(.*)$} $refalimethod temp refalimethod refalipreset
				refseq_${refalimethod}_job $refseq $refalipreset
				map_job -nohardclips 1 -paired 0 -threads $threads \
					-method $refalimethod -preset $refalipreset \
					$useali \
					$refseq \
					refseq$postfix \
					$tempfastq
			}
			job readannot_match_ref$postfix-$root-$tail -deps {
				$useali
			} -targets {
				$target
			} -vars {
				useali refseqannot postfix
			} -code {
				cg sam2tsv -f {AS ms cs} $useali | cg select -f [string_change {
					rname="@$qname"
					{rstart=if($strand eq "+",$qstart,$seqlen - $qend)}
					{rend=if($strand eq "+",$qend,$seqlen - $qstart)}
					{size=$qend - $qstart}
					{postfix="@POSTFIX@"}
					strand
					{chromosome=if($chromosome ne "*","@ANNOT@","*")}
					AS ms begin end mapquality cigar *
				} [list @ANNOT@ $refseqannot @POSTFIX@ $postfix] ] -q {$chromosome ne "*"} \
					| cg select -s {rname rstart rend} \
					| cg zst > $target.temp.zst
				file rename -force $target.temp.zst $target
			}
		}
		if {$polyT} {
			lappend toadd $workdir/$root-$tail.polyt.ali.tsv.zst
			job readannot_polyt-$root-$tail -deps {
				$tempfastq
			} -targets {
				$workdir/$root-$tail.polyt.ali.tsv.zst
			} -vars {
				fastq workdir root tail tempfastq polyT
			} -code {
				catch_exec cg zcat $tempfastq | polyt $polyT \
				    | cg zst > $workdir/$root-$tail.polyt.ali.tsv.temp.zst
				file rename -force $workdir/$root-$tail.polyt.ali.tsv.temp.zst $workdir/$root-$tail.polyt.ali.tsv.zst
			}
		}
		set target $workdir/$root-$tail.annot.tsv.zst
		job readannot_schema-$root-$tail \
		-deps [list {*}$toadd $annotationfile] \
		-targets {
			$target
		} -rmtargets {
			$tempfastq
		} -vars {
			fastq workdir root tail addsequences ignoreN polyT toadd tempfastq completeness annotationfile
		} -code {
			# read annotation sizes
			set f [gzopen $annotationfile]
			unset -nocomplain minsizea
			while 1 {
				if {[gets $f name] == -1} break
				set name [string range $name 1 end]
				set seq [gets $f]
				set size [string length $seq]
				set minsizea($name) [expr {int($completeness*$size/100.0)}]
			}
			gzclose $f

			# file delete $tempfastq
			set concat $workdir/$root-$tail.concat.ali.tsv.zst
			if {[llength $toadd] > 1} {
				set tempfile [tempfile]
				cg cat -m 1 {*}$toadd \
					| cg select -s {rname rstart -rend} | cg zst > $concat.temp.zst
			} else {
				cg select -overwrite 1 -s {rname rstart -rend} $workdir/$root-$tail.annot.ali.tsv.zst $concat.temp.zst
			}
			file rename -force $concat.temp.zst $concat
			set src $concat

			catch {gzclose $f}; catch {gzclose $o}
			set f [gzopen $src]
			set o [wgzopen $target.temp.zst]
			set oheader rname\treadsize\tschema\tshortschema\tschema2
			if {$addsequences} {append oheader \tsequences}
			puts $o $oheader
			set header [tsv_open $f]
			set poss [list_cor $header {rname rstart rend strand chromosome AS ms qstart qend seq}]
			set curname {}
			set curpos 0
			set curseq {}
			set todo {}
			if {[gets $f line] == -1} {error "error reading $src: no data"}
			set data [list_sub [split $line \t] $poss]
			set todo [list $data]
			set curname [lindex $data 0]
			set lastone 0
			while 1 {
				if {[gets $f line] == -1} {
					if {$lastone} break
					set lastone 1
				}
				set nextdata [list_sub [split $line \t] $poss]
				set nextrname [lindex $nextdata 0]
				# foreach {rname rstart rend strand chromosome AS ms seq} $data break
				if {$nextrname eq $curname} {
					lappend todo $nextdata
					continue
				}
				if {![llength $todo]} continue
				#
				# we have a full todo (all hits on one read), process
				set schema {}
				set shortschema {}
				set schema2 {}
				set sequences {}
				set curpos 0
				foreach {rname rstart rend strand chromosome AS ms qstart qend seq} [lindex $todo 0] break
				set readsize [string length $seq]
				if {$readsize == 0} {
					set seq [lindex [list_sub $todo -exclude [list_find [list_subindex $todo 4] polyT]] 0 9]
					set readsize [string length $seq]
				}
				# make schema
				if {$addsequences} {
					if {$strand eq "-"} {set seq [seq_complement $seq]}
				}
				set pos [lsearch -exact [list_subindex $todo 4] *]
				if {$pos != -1} {
					set todo [list_sub $todo -exclude $pos]
				}
				# list_subindex $todo 4
				foreach line $todo {
					foreach {rname rstart rend strand chromosome AS ms qstart qend} $line break
					# putsvars rname rstart rend strand chromosome AS ms qstart qend
					if {[regexp ^chr $chromosome]} {
						set chromosome transcript
						set strand ~
					} elseif {[regexp ^transcript $chromosome]} {
						set strand ~
					} elseif {[info exists minsizea($chromosome)] && [isint $qstart] && [isint $qend]} {
						if {$qend-$qstart < $minsizea($chromosome)} continue
					}
					set emptysize [expr {$rstart-$curpos}]
					if {$emptysize >= $ignoreN} {
						lappend schema ~ N
						lappend schema2 ~ N_$emptysize
						if {$addsequences} {
							lappend sequences ~ N [string range $seq $curpos [expr {$rstart-1}]]
						}
					}
					if {$chromosome eq [lindex $schema end]} {
						# ignore same chromosome
					} elseif {$rend <= $curpos} {
						# ignore full overlap
					} elseif {$rstart < $curpos} {
						# overlap
						lappend schema $strand $chromosome
						lappend shortschema $strand $chromosome
						lappend schema2 $strand ${chromosome}_[expr {$rend-$curpos}]
						if {$addsequences} {
							lappend sequences $strand ${chromosome} [string range $seq $curpos [expr {$rend-1}]]
						}
					} else {
						# no overlap, large empty is already done
						lappend schema $strand $chromosome
						lappend shortschema $strand $chromosome
						lappend schema2 $strand ${chromosome}_[expr {$rend-$rstart}]
						if {$addsequences} {
							lappend sequences $strand ${chromosome} [string range $seq $rstart [expr {$rend-1}]]
						}
					}
					set curpos $rend
				}
				# add anything left
				set clen [string length $seq]
				set remainder [expr {$clen-$curpos}]
				if {$remainder > 3} {
					lappend schema ~ N
					lappend schema2 ~ N_$remainder
					if {$addsequences} {
						lappend sequences ~ N [string range $seq $curpos end]
					}							
				}
				# write to output
				if {$addsequences} {
					puts $o [string range $curname 1 end]\t$readsize\t[join $schema]\t[join $shortschema]\t[join $schema2]\t$sequences
				} else {
					puts $o [string range $curname 1 end]\t$readsize\t[join $schema]\t[join $shortschema]\t[join $schema2]
				}
				set curname $nextrname
				set curpos 0
				set todo [list $nextdata]
			}

			gzclose $f; gzclose $o
			file rename -force $target.temp.zst $target
		}
		lappend alist $target
	}

	job readannot_merge-[file tail $result] -deps $alist -targets {
		$result
		$resultdir/${root}_summary.tsv
		$resultdir/${root}_shortsummary.tsv
	} -vars {
		alist root resultdir result
	} -code {
		cg cat {*}$alist {*}[compresspipe $result] > $result.temp2
		file rename -force $result.temp2 $result
		cg select -s -count -g schema -gc count,percent,q1(readsize),avg(readsize),q3(readsize) $result > $resultdir/${root}_summary.tsv.temp2
		file rename -force $resultdir/${root}_summary.tsv.temp2 $resultdir/${root}_summary.tsv
		cg select -s -count -g shortschema -gc count,percent,q1(readsize),avg(readsize),q3(readsize) $result > $resultdir/${root}_shortsummary.tsv.temp2
		file rename -force $resultdir/${root}_shortsummary.tsv.temp2 $resultdir/${root}_shortsummary.tsv
	}
}

proc readsaber {args} {
	# pick up options like -stack and -v
	set args [parse_generic_args readsaber $args]
	# process common cg job_options (-d, -dsubmit, ...)
	set args [job_init {*}$args]
	# run main command
	readsaber_job {*}$args
	# if needed (e.g. for -d <num>), wait until jobs are finished
	job_wait
}

readsaber {*}$argv
