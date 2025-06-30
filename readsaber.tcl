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
	set refseq {}
	set keepintermediate 0
	set addsequences 0
	set ignoreN 3
	if {[lindex $args] in "-h -help help"} {
		help readsaber
		exit 0
	}
	cg_options readsaber args {
		-refseq {set refseq [refseq $value]}
		-keepintermediate {set keepintermediate $value}
		-addsequences {set addsequences $value}
		-ignoreN {set ignoreN $value}
		-version - -V {
			puts "v0.1.0"
			exit 0
		}
	} {annotationfile result fastq} 3 ... {
		test
	}
	set ref [file_absolute $annotationfile]
	set result [file_absolute $result]
	set fastq [file_absolute $fastq]
	set fastqs [list [file_absolute $fastq]]
	foreach fastq $args {
		lappend fastqs [file_absolute $fastq]
	}
	# index ref
	job readannot_index-[file tail $ref] -deps {
		$ref
	} -targets {
		$ref.map-ont
	} -vars {
		ref
	} -code {
		exec samtools faidx $ref
		catch_exec minimap2 -x map-ont -k 5 -w 1 -d $ref.map-ont $ref
	}
	set root [file root [gzroot [file tail $result]]]
	set resultdir [file dir $result]
	job_logfile [file dir $result]/sc_readsaber-[file tail $result] [file dir $result] $cmdline \
		{*}[versions genomecomb samtools dbdir zstd os]

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
		job readannot_match-$tail -deps {
			$fastq
		} -targets {
			$workdir/$root-$tail.ali.tsv.zst
		} -vars {
			fastq ref refseq workdir root tail
		} -code {
			set tempfile [tempfile].fastq.gz
			if {[file extension $fastq] eq ".bam"} {
				catch_exec samtools fastq -T "RG,CB,QT,MI,MM,ML,Mm,Ml" $fastq | gzip > $tempfile
			} else {
				mklink $fastq $tempfile
			}
			catch_exec minimap2 -Y -a -x map-ont -t 4 -n 1 -m 1 -k 5 -w 1 -s 20 $ref $tempfile \
			    | cg zst > $workdir/$root-$tail.ali.sam.zst 2>@ stderr
			cg sam2tsv -f {AS ms cs} $workdir/$root-$tail.ali.sam.zst | cg select -f {
				rname=$qname
				{rstart=if($strand eq "+",$qstart,$seqlen - $qend)}
				{rend=if($strand eq "+",$qend,$seqlen - $qstart)}
				{size=$qend - $qstart}
				strand chromosome AS ms begin end mapquality cigar *
			} | cg select -s {rname rstart rend} | cg zst > $workdir/$root-$tail.ali.tsv.temp.zst
			file rename -force $workdir/$root-$tail.ali.tsv.temp.zst $workdir/$root-$tail.ali.tsv.zst
		}
		if {$refseq ne ""} {
			job readannot_match_ref-$tail -mem [map_mem_minimap2 "" 4 map-ont $refseq] -deps {
				$refseq
			} -targets {
				$workdir/$root-$tail.refseq.ali.tsv.zst
			} -vars {
				fastq ref refseq workdir root tail
			} -code {
				set tempfile [tempfile].fastq.gz
				if {[file extension $fastq] eq ".bam"} {
					catch_exec samtools fastq -T "RG,CB,QT,MI,MM,ML,Mm,Ml" $fastq | gzip > $tempfile
				} else {
					mklink $fastq $tempfile
				}
				catch_exec minimap2 -P -Y -a -x map-ont --splice -t 4 $refseq $tempfile \
				    | cg zst > $workdir/$root-$tail.refseq.ali.sam.zst 2>@ stderr
				cg sam2tsv -f {AS ms cs} $workdir/$root-$tail.refseq.ali.sam.zst | cg select -f {
					rname=$qname 
					{rstart=if($strand eq "+",$qstart,$seqlen - $qend)}
					{rend=if($strand eq "+",$qend,$seqlen - $qstart)}
					{size=$qend - $qstart}
					strand chromosome AS ms begin end mapquality cigar *
				} -q {$chromosome ne "*"} | cg select -s {rname rstart rend} | cg zst > $workdir/$root-$tail.refseq.ali.tsv.temp.zst
				file rename -force $workdir/$root-$tail.refseq.ali.tsv.temp.zst $workdir/$root-$tail.refseq.ali.tsv.zst
			}
		}
		set target $workdir/$root-$tail.tsv.zst
		job readannot_schema-$tail -deps {
			$workdir/$root-$tail.ali.tsv.zst
			($workdir/$root-$tail.refseq.ali.tsv.zst)
		} -targets {
			$target
		} -vars {
			fastq ref refseq workdir root tail addsequences ignoreN
		} -code {
			set concat $workdir/$root-$tail.concat.ali.tsv.zst
			if {$refseq ne ""} {
				cg cat $workdir/$root-$tail.ali.tsv.zst $workdir/$root-$tail.refseq.ali.tsv.zst \
					| cg select -s {rname rstart rend} | cg zst > $concat.temp.zst
			} else {
				cg select -overwrite 1 -s {rname rstart rend} $workdir/$root-$tail.ali.tsv.zst $concat.temp.zst
			}
			file rename -force $concat.temp.zst $concat
			set src $concat

			catch {gzclose $f}; catch {gzclose $o}
			set f [gzopen $src]
			set o [wgzopen $target.temp.zst]
			set oheader rname\tschema\tschema2
			if {$addsequences} {append oheader \tsequences}
			puts $o $oheader
			set header [tsv_open $f]
			set poss [list_cor $header {rname rstart rend strand chromosome AS ms seq}]
			set curname {}
			set curpos 0
			set curseq {}
			set todo {}
			if {[gets $f line] == -1} {error "error reading $src: no data"}
			set data [list_sub [split $line \t] $poss]
			set todo [list $data]
			set curname [lindex $data 0]
			while 1 {
				if {[gets $f line] == -1} break
				set nextdata [list_sub [split $line \t] $poss]
				set nextrname [lindex $nextdata 0]
				# foreach {rname rstart rend strand chromosome AS ms seq} $data break
				if {$nextrname eq $curname} {
					lappend todo $nextdata
					continue
				}
				if {![llength $todo]} continue
				set schema {}
				set schema2 {}
				set sequences {}
				set curpos 0
				foreach {rname rstart rend strand chromosome AS ms seq} [lindex $todo 0] break
				if {$addsequences} {
					if {$strand eq "-"} {set seq [seq_complement $seq]}
				}
				set pos [lsearch -exact [list_subindex $todo 4] *]
				if {$pos != -1} {
					set todo [list_sub $todo -exclude $pos]
				}
				# list_subindex $todo 4
				foreach line $todo {
					foreach {rname rstart rend strand chromosome AS ms} $line break
					# putsvars rname rstart rend strand chromosome AS ms
					if {[regexp ^chr $chromosome]} {
						set chromosome transcript
						set strand ~
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
						lappend schema2 $strand ${chromosome}_[expr {$rend-$curpos}]
						if {$addsequences} {
							lappend sequences $strand ${chromosome} [string range $seq $curpos [expr {$rend-1}]]
						}
					} else {
						# no overlap, large empty is already done
						lappend schema $strand $chromosome
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
					puts $o $curname\t[join $schema]\t[join $schema2]\t$sequences
				} else {
					puts $o $curname\t[join $schema]\t[join $schema2]
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
	} -vars {
		alist root resultdir result
	} -code {
		cg cat {*}$alist {*}[compresspipe $result] > $result.temp2
		file rename -force $result.temp2 $result
		cg select -s -count -g schema -gc count,percent $result > $resultdir/${root}_summary.tsv.temp2
		file rename -force $resultdir/${root}_summary.tsv.temp2 $resultdir/${root}_summary.tsv
	}
}

proc readsaber {args} {
	# process common cg job_options (-d, -dsubmit, ...)
	set args [job_init {*}$args]
	# run main command
	readsaber_job {*}$args
	# if needed (e.g. for -d <num>), wait until jobs are finished
	job_wait
}

readsaber {*}$argv
