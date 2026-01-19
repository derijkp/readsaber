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
		-refseqannot {lappend refseqannots $value}
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

	foreach fastq $fastqs {
		set tail [file root [gzroot [file tail $fastq]]]
		if {[string length $root-$tail.polyt.ali.tsv.temp.zst] >= 252} {
			set tail [shortenname $tail]
		}
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
				mklink -absolute 1 $fastq $tempfastq.temp
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
			fastq workdir root tail tempfastq threads
		} -code {
			cg sam2tsv -f {AS ms cs} $workdir/$root-$tail.annot.ali.tsv.sam.zst | cg select -f {
				rname="@$qname"
				{rstart=if($strand eq "+",$qstart,$seqlen - $qend)}
				{rend=if($strand eq "+",$qend,$seqlen - $qstart)}
				{size=$qend - $qstart}
				strand chromosome AS ms begin end mapquality cigar *
				{priority=1}
			} | cg select -s {rname rstart rend} | cg zst > $workdir/$root-$tail.annot.ali.tsv.temp.zst
			file rename -force $workdir/$root-$tail.annot.ali.tsv.temp.zst $workdir/$root-$tail.annot.ali.tsv.zst
		}
		set priority 1
		if {$polyT} {
			incr priority
			lappend toadd $workdir/$root-$tail.polyt.ali.tsv.zst
			job readannot_polyt-$root-$tail -skip $endtargets -skip $workdir/$root-$tail.annot.tsv.zst -deps {
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
				set useali $workdir/$root-$tail.refseq$postfix.ali.sam.zst
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
					| cg select -s {rname rstart rend} \
					| cg zst > $target.temp.zst
				file rename -force $target.temp.zst $target
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
			set poss [list_cor $header {rname rstart rend strand chromosome qstart qend seq priority}]
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
				# foreach {rname rstart rend strand chromosome seq} $data break
				if {$nextrname eq $curname} {
					lappend todo $nextdata
					continue
				}
				if {![llength $todo]} continue
				#
				# we have a full todo (all hits on one read), process
				# get sequence from first (if it does not have sequence, from another)
				foreach {rname rstart rend seqstrand chromosome qstart qend seq priority} [lindex $todo 0] break
# if {$curname eq "@bd6f13c0-9e0d-4ff4-9841-6cb59fe6dbd4"} error
				set readsize [string length $seq]
				if {$readsize == 0} {
					set sline [lindex [list_sub $todo -exclude [list_find [list_subindex $todo 4] polyT]] 0]
					set seq [lindex $sline 9]
					set seqstrand [lindex $sline 3]
					set readsize [string length $seq]
				}
				if {$seqstrand eq "+"} {
					set seqa(+) $seq
					set seqa(-) [seq_complement $seq]
				} else {
					set seqa(-) $seq
					set seqa(+) [seq_complement $seq]
				}
				set seqa(~) $seqa(+)
				set seqa() $seqa(+)
				# remove unaligned
				set pos [lsearch -exact [list_subindex $todo 4] *]
				if {$pos != -1} {
					set todo [list_sub $todo -exclude $pos]
				}

				# list_subindex $todo 4
				# make schema
				set schema {}
				set shortschema {}
				set schema2 {}
				set sequences {}
				set curpos 0
				# go over (sorted) todo
				foreach {pname pstart pend pstrand pchromosome pqstart pqend ppriority} {{} 0 0 {} {} 0 0 0} break
				foreach line $todo {
					foreach {rname rstart rend strand chromosome qstart qend temp priority} $line break
					# putsvars rname rstart rend strand chromosome qstart qend
					if {[regexp ^chr $chromosome]} {
						# set chromosome transcript
						set strand ~
					} elseif {[regexp ^transcript $chromosome]} {
						set strand ~
					} elseif {[info exists minsizea($chromosome)] && [isint $qstart] && [isint $qend]} {
						if {$qend-$qstart < $minsizea($chromosome)} continue
					}
					# putsvars sequences line chromosome pchromosome strand pstrand ppriority priority rstart rend pstart pend
					set Nsize [expr {$rstart-$pend}]
					if {$chromosome eq $pchromosome && $strand eq $pstrand && $Nsize <= $ignoreN && $ppriority != 1 && $priority != 1} {
						# ignore same chromosome if not annot (priority 1)
						# -> current will get prolonged
						# annot will always be separated (to allow finding duplicate annots)
						set pend $rend
					} elseif {$rend <= $pend} {
						# ignore full overlap
					} else {
						if {$rstart < $pend} {
							# overlap 
							if {$priority < $ppriority} {
								set pend $rstart
							} else {
								set rstart $pend
							}
							set Nsize 0
						} else {
							if {$Nsize >= $ignoreN} {
								set Nstring [string range $seqa($pstrand) $pend [expr {$rstart-1}]]
							}
						}
						if {$ppriority ne "0"} {
							lappend schema $pstrand $pchromosome
							lappend shortschema $pstrand $pchromosome
							lappend schema2 $pstrand ${pchromosome}_[expr {$pend-$pstart}]
							if {$addsequences} {
								lappend sequences $pstrand ${pchromosome} [string range $seqa($pstrand) $pstart [expr {$pend-1}]]
							}
						}
						if {$Nsize >= $ignoreN} {
							lappend schema ~ N
							lappend schema2 ~ N_$Nsize
							if {$addsequences} {
								lappend sequences ~ N $Nstring
							}
						}
						foreach {pname pstart pend pstrand pchromosome pqstart pqend ppriority} \
							[list $rname $rstart $rend $strand $chromosome $qstart $qend $priority] break
					}
					set ppriority $priority
				}

				# Add previous (lagging)
				if {[llength $todo]} {
					lappend schema $pstrand $pchromosome
					lappend shortschema $pstrand $pchromosome
					lappend schema2 $pstrand ${pchromosome}_[expr {$pend-$pstart}]
					if {$addsequences} {
						lappend sequences $pstrand ${pchromosome} [string range $seqa($pstrand) $pstart [expr {$pend-1}]]
					}
				}
				# add empty if there are left
				set clen [string length $seq]
				set remainder [expr {$clen-$pend}]
				if {$remainder > 3} {
					lappend schema ~ N
					lappend schema2 ~ N_$remainder
					if {$addsequences} {
						lappend sequences ~ N [string range $seqa($pstrand) $pend end]
					}							
				}
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
