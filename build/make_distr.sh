#!/bin/sh
# the next line restarts using wish \
exec cg source "$0" ${1+"$@"}

# currently, this script can only update a genomecomb distribution
# it expects dirtcl to be present at $dest

if {[llength $argv] > 2} {
	error "format is make_distr.sh ?-dest dest?"
}

#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

# find appdir
# -----------
set script [file join [pwd] [info script]]
while 1 {
	if {[catch {set script [file join [pwd] [file readlink $script]]}]} break
}
if {[file isdir $script]} {
	set scriptdir $script
} else {
	set scriptdir [file dir $script]
}
if {[file tail $scriptdir] eq "build"} {
	set appdir [file dir $scriptdir]
} else {
	set appdir $scriptdir
}

# parse arguments
# ---------------

set version 0.1.0

set dest ""
set genomecomb ""
set arch linux-x86_64
foreach {key value} $argv {
	switch $key {
		-genomecomb {set genomecomb $value}
		-dest {set dest $value}
		default {error "unknown option, must be one off: -dest, -genomecomb"}
	}
}
if {$arch eq "linux-x86_64"} {
	set archsuffix ""
} else {
	set archsuffix "-linux-ix86"
}

if {$dest eq ""} {
	set dest $::env(HOME)/build/readsaber-$version-$arch
}
if {$genomecomb eq ""} {
	set genomecomb [file dir $dest]/genomecomb-0.112.0-$arch
}
if {![file exists $genomecomb]} {
	error "genomecomb binary directory ($genomecomb) not found"
	if {![file exists $genomecomb/cg]} {
		error "genomecomb binary directory ($genomecomb) exists, but is not a genomecomb appdir"
	}
}

# make distribution
# -----------------
puts "build readsaber distribution in $dest (based on $genomecomb)"
if {[file exists $dest]} {
	file delete -force $dest.old
	file rename $dest $dest.old
}
file copy $genomecomb $dest
cd $appdir
file copy -force readsaber readsaber.tcl readsaber_makerefdir readsaber_makerefdir.tcl README.md help $dest
file copy -force {*}[glob bin/*] $dest/bin
