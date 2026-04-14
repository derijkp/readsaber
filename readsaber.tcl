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

proc shortenname {string} {
	join [list_remdup [split $string _-]] _
}

# main proc, using job system
proc readsaber_job {args} {
	upvar job_logdir job_logdir
	set cmdline [clean_cmdline cg readsaber {*}$args]
	set refseqs {}
	set refseqannots {}
	set keepintermediate 0
	set addsequences 0
	set alimethod minimap2_short
	set refalimethod minimap2_splice
	set ignoreN 3
	set remN {}
	set refminsize 7
	set polyT 16
	set completeness 50
	set threads 2
	if {[lindex $args] in "-h -help help"} {
		help readsaber
		exit 0
	}
	cg_options readsaber args {
		-refseq {lappend refseqs [refseq $value]}
		-refseqannot - -refseqname {lappend refseqannots $value}
		-alimethod {set alimethod $value}
		-refalimethod {set refalimethod $value}
		-refminsize {set refminsize $value}
		-keepintermediate {set keepintermediate $value}
		-addsequences {set addsequences $value}
		-ignoreN {set ignoreN $value}
		-remN {set remN $value}
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
	set aliextraopts {}
	if {$remN eq ""} {
		set remN $ignoreN
	} elseif {$remN < $ignoreN} {
		set remN $ignoreN
	}
	set annotationfile [file_absolute $annotationfile]
	set result [file_absolute $result]
	set fastq [file_absolute $fastq]
	set fastqs [list [file_absolute $fastq]]
	foreach fastq $args {
		lappend fastqs [file_absolute $fastq]
	}
	set len [llength $refseqs]
	set lena [llength $refseqannots]
	if {$lena < $len} {
		lappend refseqannots {*}[list_fill [expr {$len-$lena}] {transcript$postfix}]
	} elseif {$lena > $len} {
		set refseqannots [lrange $refseqannots 0 [expr {$len-1}]]
	}
	set root [file root [gzroot [file tail $result]]]
	set resultdir [file dir $result]
	job_logfile [file dir $result]/sc_readsaber-[file tail $result] [file dir $result] $cmdline \
		{*}[versions genomecomb samtools dbdir zstd os]

	set method $alimethod ; set alipreset {}
	regexp {(^[^_]+)_(.*)$} $alimethod temp alimethod alipreset
	set annotrefseq [refseq_${alimethod}_job $annotationfile $alipreset]
	set refalipreset {}
	regexp {(^[^_]+)_(.*)$} $refalimethod temp refalimethod refalipreset
	unset -nocomplain refa
	foreach refseq $refseqs refseqannot $refseqannots {
		set refa($refseq) [refseq_${refalimethod}_job $refseq $refalipreset]
	}

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

	set endtargets [list $result $resultdir/${root}_summary.tsv $resultdir/${root}_shortsummary.tsv]
	foreach fastq $fastqs {
		set tail [file root [gzroot [file tail $fastq]]]
		if {[string length $root-$tail.polyt.ali.tsv.temp.zst] >= 252} {
			set tail [shortenname $tail]
		}
		set perfastqfinal $workdir/$root-$tail.final.tsv.zst
		set tempfastq $workdir/$root-$tail.fastq.gz
		job maketempfastq-$root-$tail -skip $endtargets -skip $perfastqfinal -deps {
			$fastq
		} -targets {
			$tempfastq
		} -vars {
			fastq tempfastq
		} -code {
			if {[file extension $fastq] eq ".bam"} {
				catch_exec samtools fastq -T "RG,CB,QT,MI,MM,ML,Mm,Ml" $fastq | cg bgzip > $tempfastq.temp
			} elseif {[file extension $fastq] eq ".gz"} {
				mklink -absolute 1 $fastq $tempfastq.temp
			} else {
				catch_exec cg bgzip -o $tempfastq.temp $fastq
			}
			file rename $tempfastq.temp $tempfastq
		}
		set toadd [list $workdir/$root-$tail.annot.ali.tsv.zst]
		map_job -skip $endtargets -skip $perfastqfinal \
			-nohardclips 1 -paired 0 -threads $threads \
			-method $alimethod -preset $alipreset \
			-extraopts $aliextraopts \
			-sort nosort \
			$workdir/$root-$tail.annot.ali.tsv.sam.zst \
			$annotationfile \
			annotation \
			$tempfastq
		job readannot_match-$root-$tail -skip $endtargets -skip $perfastqfinal -cores $threads -deps {
			$workdir/$root-$tail.annot.ali.tsv.sam.zst
		} -targets {
			$workdir/$root-$tail.annot.ali.tsv.zst
		} -vars {
			fastq workdir root tail tempfastq threads
		} -code {
			cg sam2tsv -f {AS ms cs} $workdir/$root-$tail.annot.ali.tsv.sam.zst | cg select -f {
				rname="@$qname"
				{rstart=if($strand eq "+",$qstart,$seqlen - $qend)}
				{rend=if($strand eq "+",$qend,$seqlen - $qstart)}
				{size=$qend - $qstart}
				strand chromosome AS ms begin end mapquality cigar *
				{priority=1}
			} | cg zst > $workdir/$root-$tail.annot.ali.tsv.temp.zst
			file rename -force $workdir/$root-$tail.annot.ali.tsv.temp.zst $workdir/$root-$tail.annot.ali.tsv.zst
		}
		set priority 2
		if {$polyT} {
			lappend toadd $workdir/$root-$tail.polyt.ali.tsv.zst
			job readannot_polyt-$root-$tail -skip $endtargets -skip $perfastqfinal -deps {
				$tempfastq
			} -targets {
				$workdir/$root-$tail.polyt.ali.tsv.zst
			} -vars {
				fastq workdir root tail tempfastq polyT priority
			} -code {
				catch_exec cg zcat $tempfastq | polyt $polyT \
					| cg select -f [list * "priority=$priority"] \
					| cg zst > $workdir/$root-$tail.polyt.ali.tsv.temp.zst
				file rename -force $workdir/$root-$tail.polyt.ali.tsv.temp.zst $workdir/$root-$tail.polyt.ali.tsv.zst
			}
		}
		set refnr 0
		foreach refseq $refseqs refseqannot $refseqannots {
			incr priority
			incr refnr
			if {$refnr == 1} {set postfix ""} else {set postfix $refnr}
			set target $workdir/$root-$tail.refseq$postfix.ali.tsv.zst
			lappend toadd $target
			if {[file extension $refseq] in ".cram .bam .sam"} {
				set useali $refseq
			} else {
				set useali $workdir/$root-$tail.refseq$postfix.${refalimethod}_$refalipreset.ali.sam.zst
				map_job -skip $endtargets -skip $perfastqfinal \
					-nohardclips 1 -paired 0 -threads $threads \
					-method $refalimethod -preset $refalipreset \
					-sort nosort \
					$useali \
					$refseq \
					refseq$postfix \
					$tempfastq
			}
			job readannot_match_ref$postfix-$root-$tail -skip $endtargets -skip $perfastqfinal -deps {
				$useali
			} -targets {
				$target
			} -vars {
				useali refseqannot postfix priority
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
					{priority=@priority@}
				} [list @ANNOT@ $refseqannot @POSTFIX@ $postfix @priority@ $priority] ] -q {$chromosome ne "*"} \
					| cg zst > $target.temp.zst
				file rename -force $target.temp.zst $target
			}
		}
		set concat $workdir/$root-$tail.concat.ali.tsv.zst
		job readannot_concat-$root-$tail -skip $endtargets \
		-deps $toadd \
		-targets {
			$concat
		} -vars {
			toadd concat
		} -code {
			if {[llength $toadd] > 1} {
				set tempfile [tempfile]
				cg cat -m 1 {*}$toadd \
					| cg select -s {rname rstart -rend} | cg zst > $concat.temp.zst
			} else {
				cg select -overwrite 1 -s {rname rstart -rend} $workdir/$root-$tail.annot.ali.tsv.zst $concat.temp.zst
			}
			file rename -force $concat.temp.zst $concat
		}

		set target $perfastqfinal
		job readannot_schema-$root-$tail \
		-deps [list {*}$concat $annotationfile] \
		-targets {
			$target
		} -rmtargets {
			$tempfastq
		} -vars {
			concat addsequences ignoreN remN polyT completeness annotationfile refminsize
		} -code {
			# read annotation sizes
			set f [gzopen $annotationfile]
			unset -nocomplain minsizea
			while 1 {
				if {[gets $f name] == -1} break
				set name [lindex [string range $name 1 end] 0]
				set seq [gets $f]
				set size [string length $seq]
				set minsizea($name) [expr {int($completeness*$size/100.0)}]
			}
			gzclose $f

			# file delete $tempfastq
			set src $concat

			catch {gzclose $f}; catch {gzclose $o}
			set f [gzopen $src]
			set o [wgzopen $target.temp.zst]
			set oheader rname\treadsize\tschema\tshortschema\tschema2
			if {$addsequences} {append oheader \tsequences}
			puts $o $oheader
			set header [tsv_open $f]
			set poss [list_cor $header {rstart rend strand chromosome qstart qend priority seq}]
			set rnamepos [lsearch $header rname]
			set curname {}
			set curpos 0
			set curseq {}
			set todo {}
			if {[gets $f line] == -1} {error "error reading $src: no data"}
			set split [split $line \t]
			set data [list_sub $split $poss]
			set todo [list $data]
			set curname [lindex $split $rnamepos]
			set lastone 0

			while 1 {
				if {[gets $f line] == -1} {
					if {$lastone} break
					set lastone 1
				}
				set split [split $line \t]
				set nextdata [list_sub $split $poss]
				set nextrname [lindex $split $rnamepos]
				# foreach {rname rstart rend strand chromosome seq} $data break
				if {$nextrname eq $curname} {
					lappend todo $nextdata
					continue
				}
				if {![llength $todo]} continue
				# puts [join [list_subindex $todo 0 1 2 3 4 5 6] \n]
				#
				# we have a full todo (all hits on one read), process
				# get sequence from first (if it does not have sequence, from another)
				foreach {rstart rend seqstrand chromosome qstart qend priority seq} [lindex $todo 0] break
				#
# if {$curname eq "@seq7"} {error testing_interrupt}
				# find a non-empty sequence
				# -------------------------
				set seqs [list_subindex $todo 7]
				set strands [list_subindex $todo 2]
				foreach seq $seqs strand $strands {
					if {$seq ni {* {}}} break
				}
				if {$strand eq "-"} {
					set seq [seq_complement $seq]
				}
				# recalculate rstart and rend where needed (* sequence was given in alignment)
				# ----------------------------------------
				set readsize [string length $seq]
				if {"*" in $seqs} {
					set pos 0
					set changes 0
					foreach oseq $seqs strand $strands {
						if {$oseq eq "*" && $strand eq "-"} {
							set line [lindex $todo $pos]
							set rstart [expr {$readsize - [lindex $line 5]}]
							set rend [expr {$readsize - [lindex $line 4]}]
							lset line 0 $rstart
							lset line 1 $rend
							lset todo $pos $line
							incr changes
						}
						incr pos
					}
					if {$changes} {
						set todo [lsort -index 0 -integer [lsort -index 1 -integer -decreasing $todo]]
					}
				}
				# remove unaligned
				set pos [lsearch -exact [list_subindex $todo 3] *]
				if {$pos != -1} {
					set todo [list_sub $todo -exclude $pos]
				}

				#
				# cut overlaps (store in annots)
				# ------------
				# list_subindex $todo 4
				set annots {}
				set curpos 0
				# go over (sorted) todo
				foreach {pstart pend pstrand pchromosome pqstart pqend ppriority} {0 0 {} {} 0 0 0} break
				set pos 0
				set len [llength $todo]
				while {$pos < $len} {
					set line [lindex $todo $pos]
					incr pos
					foreach {rstart rend strand chromosome qstart qend priority} $line break
					# putsvars rname rstart rend strand chromosome qstart qend
					if {[regexp ^chr $chromosome]} {
						# set chromosome transcript
						set strand ~
					} elseif {[regexp ^transcript $chromosome]} {
						set strand ~
					} elseif {[info exists minsizea($chromosome)] && [isint $qstart] && [isint $qend]} {
						if {$qend-$qstart < $minsizea($chromosome)} continue
					}
# if {$chromosome eq "polyT"} error
					# putsvars sequences line chromosome pchromosome strand pstrand ppriority priority rstart rend pstart pend
					set Nsize [expr {$rstart-$pend}]
					if {$chromosome eq $pchromosome && $strand eq $pstrand && $Nsize <= $ignoreN && $ppriority != 1 && $priority != 1} {
						# ignore same chromosome if not annot (priority 1)
						# -> current will get prolonged
						# annot will always be separated (to allow finding duplicate annots)
						if {$rend > $pend} {set pend $rend}
					} elseif {$rend <= $pend} {
						# ignore full overlap unless it is higher priority
						if {$priority < $ppriority} {
							if {$chromosome eq "polyT" && [expr {$rstart-$pstart}] >= 2 && [expr {$pend-$rend}] >= 2} {
								# special case: a polyT in the middle of a transcript is ignored
								# but is kept if it is at the ends
							} else {
								if {$pend > $rend} {
									set temp [list [list $rend $pend $pstrand $pchromosome {} {} [expr {$pend - $rend}] $ppriority] \
										{*}[lrange $todo $pos end]]
									set temp [lsort -index 0 -integer [lsort -index 1 -integer -decreasing $temp]]
									set todo [list {*}[lrange $todo 0 [expr {$pos-1}]] {*}$temp]
								}
								set psize [expr {$rstart-$pstart}]
								if {$psize > 0} {
									if {$ppriority > 2 && $psize < $refminsize} {set pchromosome N ; set pstrand ~}
									lappend annots [list $pstart $rstart $pstrand $psize $pchromosome $ppriority]
								}
								set pend $rstart
								foreach {pstart pend pstrand pchromosome pqstart pqend ppriority} \
									[list $rstart $rend $strand $chromosome $qstart $qend $priority] break
							}
						}
					} else {
						if {$rstart < $pend} {
							# overlap 
							if {$priority < $ppriority} {
								set pend $rstart
							} else {
								set rstart $pend
							}
							set Nsize 0
						}
						if {$ppriority ne "0"} {
							set psize [expr {$pend-$pstart}]
							if {$psize > 0} {
								if {$ppriority > 2 && $psize < $refminsize} {set pchromosome N ; set pstrand ~}
								lappend annots [list $pstart $pend $pstrand $psize $pchromosome $ppriority]
							}
						}
						if {$Nsize > 0} {
							lappend annots [list $pend $rstart ~ $Nsize N]
						}
						foreach {pstart pend pstrand pchromosome pqstart pqend ppriority} \
							[list $rstart $rend $strand $chromosome $qstart $qend $priority] break
					}
					set ppriority $priority
				}
				# Add previous (lagging)
				set psize [expr {$pend-$pstart}]
				if {$psize > 0} {
					if {$ppriority > 2 && $psize < $refminsize} {set pchromosome N ; set pstrand ~}
					lappend annots [list $pstart $pend $pstrand $psize $pchromosome $ppriority]
				}
				# add empty if there are left
				set clen [string length $seq]
				set remainder [expr {$clen-$pend}]
				if {$remainder > 0} {
					lappend annots [list $pend $clen ~ $remainder N]
				}
				# join $annots \n
				#
				# join Ns in annots (some created from ref)
				# -----------------
				set temp {}
				set pline {}
				set ptype {}
				foreach line $annots {
					set type [lindex $line 4]
					if {$type eq "N"} {
						if {$ptype eq "N"} {
							set nend [lindex $line 1]
							lset pline 1 $nend
							lset pline 3 [expr {$nend - [lindex $pline 0]}]
						} else {
							set pline $line ; set ptype N
						}
					} else {
						if {[llength $pline]} {
							lappend temp $pline
							set pline {}
						}
						lappend temp $line
					}
					set ptype $type
				}
				if {[llength $pline]} {lappend temp $pline}
				set annots $temp
				# join $annots \n
				#
				# make schemas
				# ------------
				set schema {}
				set shortschema {}
				set schema2 {}
				set sequences {}
				foreach line $annots {
					foreach {start end strand size chromosome priority} $line break
					if {$chromosome ne "N" || $size >= $remN} {
						lappend schema $strand $chromosome
					}
					if {$chromosome ne "N"} {
						lappend shortschema $strand $chromosome
					}
					lappend schema2 $strand ${chromosome}_$size
					if {$addsequences} {
						lappend sequences $strand ${chromosome} [string range $seq $start [expr {$end-1}]]
					}
				}
				# putsvars schema schema2 shortschema sequences
				#
				# write to output
				if {$addsequences} {
					puts $o [string range $curname 1 end]\t$readsize\t[join $schema]\t[join $shortschema]\t[join $schema2]\t$sequences
				} else {
					puts $o [string range $curname 1 end]\t$readsize\t[join $schema]\t[join $shortschema]\t[join $schema2]
				}
				set curname $nextrname
				set pend 0
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
