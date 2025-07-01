package require Extral
catch {tk appname test}

package require pkgtools
namespace import -force pkgtools::*
package require Extral

set test_cleantmp 1

# pkgtools::testleak 100

set keeppath $::env(PATH)

file mkdir $::testdir
file mkdir $::testdir/tmp

proc test_cleantmp {} {
	foreach file [list_remove [glob -nocomplain $::testdir/tmp/* $::testdir/tmp/.*] $::testdir/tmp/.. $::testdir/tmp/.] {
		catch {file attributes $file -permissions ugo+xw}
		rm --recursive 1 $file
	}
	foreach file [list_remove [glob -nocomplain tmp/* tmp/.*] tmp/.. tmp/.] {
		catch {file attributes $file -permissions ugo+xw}
		rm --recursive 1 $file
	}
	cg indexclean
}

proc test {args} {
	set testdir $::testdir
	file mkdir $::testdir
	file mkdir $::testdir/tmp
	cd $::testdir
	if {![file exists $::testdir/data]} {
		mklink $::appdir/tests/data $::testdir/data
	}
	if {[get ::test_cleantmp 1]} {test_cleantmp}
	catch {job_init}
	set description [lindex $args 1]
	append description " ($::testdir)"
	lset args 1 $description
	pkgtools::test {*}$args
	cd $::appdir/tests
	return {}
}
