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

# procs
# -----
proc shortenname {string} {
	join [list_remdup [split $string _-]] _
}

proc simplifyschema {schema args} {
	set removeelements {}
	set trimelements {}
	if {[llength $args]} {
		foreach {removeelements trimelements} $args break
	}
	if {[llength $removeelements]} {set useremoveelements 1} else {set useremoveelements 0}
	if {[llength $trimelements]} {set usetrimelements 1} else {set usetrimelements 0}
	set simpleschema {}
	unset -nocomplain a
	set prev ""
	foreach {strand type} $schema {
		if {$useremoveelements && $type in $removeelements} continue
		if {$prev ne "" && $type eq $prev} continue
		incr a($strand)
		lappend simpleschema $type
		set prev $type
	}
	if {$usetrimelements} {
		set pos 0
		foreach el $simpleschema {
			if {$el ni $trimelements} break
			incr pos
		}
		set simpleschema [lrange $simpleschema $pos end]
		set pos [llength $simpleschema]
		foreach el [list_reverse $simpleschema] {
			incr pos -1
			if {$el ni $trimelements} break
		}
		set simpleschema [lrange $simpleschema 0 $pos]
	}
	if {[get a(-) 0] < [get a(+) 0]} {
		set simpleschema [list_reverse $simpleschema]
	}
	return $simpleschema
}

# main proc, using job system
# ---------------------------
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
	set simplify_remove {N}
	set simplify_trim {}
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
		-simplify_remove {set simplify_remove $value}
		-simplify_trim {set simplify_trim $value}
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
	job_logfile [file dir $result]/readsaber-[file tail $result] [file dir $result] $cmdline \
		{*}[versions genomecomb samtools dbdir zstd os]

	set method $alimethod ; set alipreset {}
	regexp {(^[^_]+)_(.*)$} $alimethod temp alimethod alipreset
#	if {$alimethod eq "minimap2"} {
#		lappend aliextraopts -p 0.5 -N 500
#	}
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

	set endtargets [list $result $resultdir/${root}_summary.tsv $resultdir/${root}_shortsummary.tsv $resultdir/${root}_simplesummary.tsv]
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
			simplify_remove simplify_trim
		} -procs {
			simplifyschema
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
			set o [wgzopen $target.temp.zst]
			set oheader rname\treadsize\tschema\tshortschema\tsimpleschema\tschema2
			if {$addsequences} {append oheader \tsequences}
			set f [gzopen $src]
			puts $o $oheader
			set header [tsv_open $f]
			set poss [list_cor $header {rstart rend strand chromosome qstart qend mapquality priority seq}]
			set seqpos 8
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
				foreach {rstart rend seqstrand chromosome qstart qend mapquality priority seq} [lindex $todo 0] break
				#
# if {$curname eq "@00ab4554-1acd-43b0-b4ae-f7aa44c2745d"} {error testing_interrupt}
				# find a non-empty sequence
				# -------------------------
				set seqs [list_subindex $todo $seqpos]
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
				foreach {pstart pend pstrand pchromosome pqstart pqend pmapquality ppriority} {0 0 {} {} 0 0 0 0} break
				set pos 0
				set len [llength $todo]
				while {$pos < $len} {
					set line [lindex $todo $pos]
					incr pos
					foreach {rstart rend strand chromosome qstart qend mapquality priority} $line break
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
					# putsvars sequences line chromosome pchromosome strand pstrand pmapquality mapquality ppriority priority rstart rend pstart pend
					# puts [list $pstart $pend $pstrand $psize $pchromosome $ppriority]
					# puts [list $rstart $rend $strand $size $chromosome $priority]
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
									set temp [list [list $rend $pend $pstrand $pchromosome {} {} [expr {$pend - $rend}] $pmapquality $ppriority] \
										{*}[lrange $todo $pos end]]
									set temp [lsort -index 0 -integer [lsort -index 1 -integer -decreasing $temp]]
									set todo [list {*}[lrange $todo 0 [expr {$pos-1}]] {*}$temp]
									set len [llength $todo]
								}
								set psize [expr {$rstart-$pstart}]
								if {$psize > 0} {
									if {$ppriority > 2 && $psize < $refminsize} {set pchromosome N ; set pstrand ~}
									lappend annots [list $pstart $rstart $pstrand $psize $pchromosome $pmapquality $ppriority]
								}
								set pend $rstart
								foreach {pstart pend pstrand pchromosome pqstart pqend pmapquality ppriority} \
									[list $rstart $rend $strand $chromosome $qstart $qend $mapquality $priority] break
							}
						} elseif {$priority == $ppriority 
							&& [isint $pmapquality] && [isint $mapquality] && $mapquality > $pmapquality
						} {
							set pchromosome $chromosome
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
								lappend annots [list $pstart $pend $pstrand $psize $pchromosome $mapquality $ppriority]
							}
						}
						if {$Nsize > 0} {
							lappend annots [list $pend $rstart ~ $Nsize N]
						}
						foreach {pstart pend pstrand pchromosome pqstart pqend pmapquality ppriority} \
							[list $rstart $rend $strand $chromosome $qstart $qend $mapquality $priority] break
					}
					set ppriority $priority
					set pmapquality $mapquality
				}
				# puts [list $pstart $pend $pstrand $psize $pchromosome $ppriority]
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
				set simpleschema [simplifyschema $schema $simplify_remove $simplify_trim]
				# putsvars schema schema2 shortschema sequences
				#
				# write to output
				if {$addsequences} {
					puts $o [string range $curname 1 end]\t$readsize\t$schema\t$shortschema\t$simpleschema\t$schema2\t$sequences
				} else {
					puts $o [string range $curname 1 end]\t$readsize\t$schema\t$shortschema\t$simpleschema\t$schema2
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
		$resultdir/${root}_simplesummary.tsv
	} -vars {
		alist root resultdir result
	} -code {
		cg cat {*}$alist {*}[compresspipe $result] > $result.temp2
		file rename -force $result.temp2 $result
		cg select -s -count -g schema -gc count,percent,q1(readsize),avg(readsize),q3(readsize) $result > $resultdir/${root}_summary.tsv.temp2
		file rename -force $resultdir/${root}_summary.tsv.temp2 $resultdir/${root}_summary.tsv
		cg select -s -count -g shortschema -gc count,percent,q1(readsize),avg(readsize),q3(readsize) $result > $resultdir/${root}_shortsummary.tsv.temp2
		file rename -force $resultdir/${root}_shortsummary.tsv.temp2 $resultdir/${root}_shortsummary.tsv
		cg select -s -count -g simpleschema -gc count,percent,q1(readsize),avg(readsize),q3(readsize) $result > $resultdir/${root}_simplesummary.tsv.temp2
		file rename -force $resultdir/${root}_simplesummary.tsv.temp2 $resultdir/${root}_simplesummary.tsv
	}

	set specialelements {}
	job readsaber_graphs-[file tail $result] -procs {readsaber_graph} -deps {
		$annotationfile
		$resultdir/${root}_summary.tsv
		$resultdir/${root}_simplesummary.tsv
	} -targets {
		$resultdir/${root}_summary.png
		$resultdir/${root}_simplesummary.png
	} -vars {
		root resultdir result specialelements annotationfile
	} -code {
		readsaber_graph \
			-usesimpleschema 0 \
			-annotationfile $annotationfile \
			-specialelements $specialelements \
			$resultdir/${root}_summary.tsv \
			$resultdir/${root}_summary.png
		readsaber_graph \
			-usesimpleschema 1 \
			-annotationfile $annotationfile \
			-specialelements $specialelements \
			$resultdir/${root}_simplesummary.tsv \
			$resultdir/${root}_simplesummary.png
	}
}

proc readsaber_graph {args} {
	set usesimpleschema 0
	set annotationfile {}
	set specialelements {N transcript polyT}
	cg_options readsaber_graph args {
		-usesimpleschema {set usesimpleschema $value}
		-annotationfile {set annotationfile $value}
		-specialelements {
			set specialelements [list_remdup [list N {*}$value]]
		}
	} {file result} 2 2 {
		Make readsaber graph
	}
	if {$annotationfile eq ""} {
		set annotnames {}
	} else {
		set annotnames [split [string trim [cg fasta2tsv $annotationfile | cg select -sh /dev/null -f id]] \n]
	}
	if {$usesimpleschema} {
		set elements [list_remove [list_subindex [lrange [split [cg select -g -simpleschema $file] \n] 1 end] 0] + - ~]
	} else {
		set elements [list_remove [list_subindex [lrange [split [cg select -g -schema $file] \n] 1 end] 0] + - ~]
	}
	set specialelements [list_union $specialelements $elements]
	set specialelements [list_lremove $specialelements $annotnames]
	R -vars {file annotationfile result usesimpleschema} -listvars {specialelements} {
	
		library(tidyverse)
		library(patchwork)

		annotation <- read_tsv(file)%>%
			filter(percent >= 1)
		if (usesimpleschema) {
			annotation$simpleschema[is.na(annotation$simpleschema)] = ""
			annotation$schema = annotation$simpleschema
		}
		layout <- read_tsv(annotationfile, col_names = FALSE)
		fasta <- tibble(name = layout$X1[seq(1,nrow(layout)-1, by = 2)],
			seq = layout$X1[seq(2,nrow(layout), by = 2)]) %>%
			mutate(l = nchar(seq),
			name = str_remove(name, ">")) %>%
			dplyr::select(-seq)
		offset <- mean(annotation$count) * 0.5

		color_palette <- c(
			"#0072B2",  # Blue
			"#E69F00",  # Orange
			"#56B4E9",  # Sky Blue
			"#D55E00",  # Vermilion
			"#009E73",  # Bluish Green
			"#CCBB44",  # Gold
			"#999999",  # Gray
			"#CC79A7",  # Reddish Purple
			"#F0E442",  # Yellow
			"#BB5500",  # Burnt sienna
			"#DDCC77",  # Mustard
			"#BBBBBB",  # Silver
			"#EE6677",  # Rose
			"#661100",  # Brown
			"#AA3377",  # Magenta
			"#228833",  # Forest green
			"#882255",  # Dark Red
			"#117733",  # Olive
			"#AA4499",  # Coral
			"#DDDDDD",  # Light Gray
			"#44AA99",  # Teal
			"#332288"  # Indigo
		)
		all_elements <- c(specialelements, fasta$name)
		color_map <- setNames(
			color_palette,
			all_elements
		)

		p1 <- annotation %>%
			mutate(schema = factor(schema, levels = rev(schema))) %>%
			filter(percent >= 1) %>%
			mutate(offset = count * 0.02) %>%
			ggplot(aes(count, schema)) +
			geom_histogram(stat = "identity", fill = "black") +
			geom_text(aes(label = sprintf("%6.2f%%", percent)), x = 0, hjust = 0, family = "Courier New") +
			theme_minimal() +
			xlab("number of reads") +
			ylab("") +
			theme(axis.text.y = element_blank(),
				panel.grid.major.y = element_blank(),
			plot.margin = margin(5, 5, 5, 5)) +
			scale_x_reverse(limits = c(max(annotation$count) + max(annotation$count)*0.2, 0)) +
			coord_cartesian(clip = "off")
		if (!usesimpleschema) {
			p2 <- annotation %>%
				mutate(schema = factor(schema, levels = rev(schema)))  %>%
				dplyr::select(schema) %>%
				mutate(split = str_split(schema, " [~+-] "),
					symbols = str_extract_all(schema, "[~+-]")) %>%
				unnest(c(split, symbols)) %>%
				rowwise() %>%
				mutate(split = str_remove(split, "[~+-] ")) %>% 
				left_join(fasta, by = c("split" = "name")) %>%
				mutate(l = ifelse(split == "N" | split == "polyT", 10, 
					ifelse(split %in% specialelements, 50, l))) %>%
				group_by(schema) %>%
				mutate(start = cumsum(c(0, rev(rev(l)[-1]))),
					end = start + l, 
					middle = (start + end) / 2) %>%
				ggplot() +
				geom_segment(aes(y = schema, yend = schema, x = start, xend =end, col = split), size = 6) +
				geom_segment(data = . %>% filter(symbols == "+"), aes(y = schema, yend = schema, x = start, xend =end), col = "white", arrow = arrow(length = unit(0.5,"cm")), linewidth = 1) +
				geom_segment(data = . %>% filter(symbols == "-"), aes(y = schema, yend = schema, x = end, xend =start), col = "white", arrow = arrow(length = unit(0.5,"cm")), linewidth = 1) +
				geom_text(aes(label = split, y = schema, x = middle), size = 3, fontface = "bold", color =   "black") +
				theme_void() +
				theme(legend.position = "none",plot.margin = margin(5, 5, 5, 2)) +
				scale_color_manual(values = color_map)
		} else {
			p2 <- annotation %>%
				mutate(schema = factor(schema, levels = rev(schema)))  %>%
				dplyr::select(schema) %>%
				mutate(split = str_split(schema, " ")) %>%
				unnest(c(split)) %>%
				rowwise() %>%
				mutate(split = str_remove(split, "[~+-] ")) %>% 
				left_join(fasta, by = c("split" = "name")) %>%
				mutate(l = ifelse(split == "N" | split == "polyT", 10, 
					ifelse(split %in% specialelements, 50, l))) %>%
				group_by(schema) %>%
				mutate(start = cumsum(c(0, rev(rev(l)[-1]))),
					end = start + l, 
					middle = (start + end) / 2) %>%
				ggplot() +
				geom_segment(aes(y = schema, yend = schema, x = start, xend =end, col = split), size = 6) +
				geom_text(aes(label = split, y = schema, x = middle), size = 3, fontface = "bold", color =   "black") +
				theme_void() +
				theme(legend.position = "none",plot.margin = margin(5, 5, 5, 2)) +
				scale_color_manual(values = color_map)
		}
		n_rows <- sum(annotation$percent >= 1)
		row_height <- 0.25  # height in inches per row, adjust to taste
		plot_height <- n_rows * row_height + 1
		p1 + p2 + plot_layout(widths = c(3, 5))
		ggsave(result,width = 15, height = plot_height, dpi = 300, bg = "white")
	}
}

proc readsaber {args} {
	# pick up options like -stack and -v
	# check for subcommands
	switch [lindex $args 0] {
		graph {
			if {[lindex $args 1] in "-h -help help"} {
				help readsaber_graph
				exit 0
			}
			set args [parse_generic_args readsaber [lrange $args 1 end]]
			readsaber_graph {*}$args
			return
		}
	}
	set args [parse_generic_args readsaber $args]
	# process common cg job_options (-d, -dsubmit, ...)
	set args [job_init {*}$args]
	# run main command
	readsaber_job {*}$args
	# if needed (e.g. for -d <num>), wait until jobs are finished
	job_wait
}

readsaber {*}$argv
