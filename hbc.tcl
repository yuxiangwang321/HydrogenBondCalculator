# HydrogenBondCalculator, version 20230523
# lengthofDA, angleofDHA, and lengthofHA are preserved variables which can be used in the construction of weight function
# note; angleofDHA is measured by setting the position of hydrogen atom as the middle point
# command example 1: source hbc.tcl; hbc test.psf test.dcd xz 0 0 1 {{0 0 10} {23 13 23} {23 1 13}} {type OP HPA OPA} 4.0 90.0 1
# mode = 1, HB density
# mode = 2, for every hydrogen bond, times a weight function: ($cutoff-$lengthofDA)/$cutoff

proc hbc {file_1 file_2 projection firstFrame lastFrame step cell atomSelection cutoff angle mode} { 

if {[string equal $file_1 $file_2]} {
	mol new $file_1 first $firstFrame last $lastFrame step $step waitfor all
} else {
	mol new $file_1
	mol addfile $file_2 first $firstFrame last $lastFrame step $step waitfor all
}

set binNumX [lindex $cell 2 0]
set binNumY [lindex $cell 2 1]
set binNumZ [lindex $cell 2 2]

set xl [lindex $cell 0 0]
set xh [lindex $cell 1 0]
set yl [lindex $cell 0 1]
set yh [lindex $cell 1 1]
set zl [lindex $cell 0 2]
set zh [lindex $cell 1 2]

set dx [expr {([lindex $cell 1 0] - [lindex $cell 0 0])/$binNumX}]
set dy [expr {([lindex $cell 1 1] - [lindex $cell 0 1])/$binNumY}]
set dz [expr {([lindex $cell 1 2] - [lindex $cell 0 2])/$binNumZ}]
set binArea [expr {$dx*$dy*$dz}]


set hbAll [open hbAll.dat w]
puts $hbAll "# Command: hbc $file_1 $file_2 $projection $firstFrame $lastFrame $step $cell $atomSelection $cutoff $angle $mode"
puts $hbAll "# Note: The first two columns are X-Y, or X-Z, or Y-Z coordinates for plane, the third column are the hb values."
puts $hbAll "# dx=$dx dy=$dy dz=$dz"

set hbNumRecord [open hbNumRecord.dat w]
puts $hbNumRecord "# frame	num"
set hbNumTotal 0

set nFrame [molinfo top get numframes]
for {set f 0} {$f<$nFrame} {incr f} {
	puts $hbAll "# "
	puts $hbAll "# frame $f"
	puts "progress:  frame $f / $nFrame"
	set hbList [measure hbonds $cutoff $angle [atomselect top "$atomSelection" frame $f] ]
	set hbNum [llength [lindex $hbList 0]]
	set hbNumTotal [expr {$hbNumTotal+$hbNum}]
	puts $hbNumRecord "$f $hbNum"

	for {set i 0} {$i<$hbNum} {incr i 1} {
		set donorAtom [atomselect top "index [lindex $hbList 0 $i]" frame $f]
		set acceptorAtom [atomselect top "index [lindex $hbList 1 $i]" frame $f]
		set hydrogenAtom [atomselect top "index [lindex $hbList 2 $i]" frame $f]
		set donorAtomID [lindex $hbList 0 $i]
		set acceptorAtomID [lindex $hbList 1 $i]
		set hydrogenAtomID [lindex $hbList 2 $i]
		set posDonor [$donorAtom get {x y z}]
		set posAcceptor [$acceptorAtom get {x y z}]
		set posHydrogen [$hydrogenAtom get {x y z}]
		set posXofHB [expr {([lindex $posHydrogen 0 0] + [lindex $posAcceptor 0 0])/2.0}]
		set posYofHB [expr {([lindex $posHydrogen 0 1] + [lindex $posAcceptor 0 1])/2.0}]
		set posZofHB [expr {([lindex $posHydrogen 0 2] + [lindex $posAcceptor 0 2])/2.0}]
		set lengthofDA [measure bond [list $donorAtomID $acceptorAtomID]  frame $f]
		set lengthofHA [measure bond [list $hydrogenAtomID $acceptorAtomID] frame $f]
		set angleofDHA [measure angle [list $donorAtomID $hydrogenAtomID $acceptorAtomID] frame $f]
		puts $hbAll "# D-H...A index $donorAtomID $acceptorAtomID $hydrogenAtomID"
		puts $hbAll "# angle = $angleofDHA"
		puts $hbAll "# lengthofDA = $lengthofDA"
		puts $hbAll "# posDonor = $posDonor"
		
		if {$mode == 1} {
			set weight 1.0
			
		} elseif {$mode == 2} {
			set weight [expr {($cutoff-$lengthofDA)/$cutoff}]
			
		} elseif {$mode ==3} {
# ---------------------- user can define their weight function here ----------------------
			# set weight xxxxxx		
		}
		set strengthHB [expr {1.0*$weight}]
		# change the positions to bin
		set binX [expr {int($posXofHB/$dx)*$dx}]
		set binY [expr {int($posYofHB/$dy)*$dy}]
		set binZ [expr {int($posZofHB/$dz)*$dz}]
		
# ---------------------------------------- HB on XZ plane ----------------------------------------
		if {[string equal -nocase $projection xz]} {	
			puts $hbAll "$binX	$binZ	$strengthHB"	
			
# ---------------------------------------- HB on YZ plane ----------------------------------------		
		} elseif {[string equal -nocase $projection yz]} {
			puts $hbAll "$binY	$binZ	$strengthHB"	
			
# ---------------------------------------- HB on XY plane ----------------------------------------
		} elseif {[string equal -nocase $projection xy]} {
			puts $hbAll "$binX	$binY	$strengthHB"	
			
# ---------------------------------------- HB on X axis ----------------------------------------
		} elseif {[string equal -nocase $projection x]} {	
			puts $hbAll "$binX	$strengthHB"
			
# ---------------------------------------- HB on Y axis ----------------------------------------		
		} elseif {[string equal -nocase $projection y]} {
			puts $hbAll "$binY	$strengthHB"	
			
# ---------------------------------------- HB on Z axis ----------------------------------------
		} elseif {[string equal -nocase $projection z]} {
			puts $hbAll "$binZ	$strengthHB"		
		}
	}
}

puts $hbNumRecord "# "
puts $hbNumRecord "# hbNumTotal = $hbNumTotal"
puts $hbNumRecord "# hbNumAvg = [expr {$hbNumTotal/$f}]"
close $hbNumRecord

puts $hbAll "# "
puts $hbAll "# empty bin"
# add the coordinates of all bins with 0 value
# ---------------------------------------- HB on XZ plane ----------------------------------------
if {[string equal -nocase $projection xz]} {
	for {set k 0} {$k<$binNumX} {incr k 1} {
		for {set m 0} {$m<$binNumZ} {incr m 1} {
			puts $hbAll "[expr {$xl+$dx*$k}]	[expr {$zl+$dz*$m}]	0.0"	
		}
	}
# ---------------------------------------- HB on YZ plane ----------------------------------------		
} elseif {[string equal -nocase $projection yz]} {
	for {set k 0} {$k<$binNumY} {incr k 1} {
		for {set m 0} {$m<$binNumZ} {incr m 1} {
			puts $hbAll "[expr {$yl+$dy*$k}]	[expr {$zl+$dz*$m}]	0.0"	
		}
	}	
	
# ---------------------------------------- HB on XY plane ----------------------------------------
} elseif {[string equal -nocase $projection xy]} {
	for {set k 0} {$k<$binNumX} {incr k 1} {
		for {set m 0} {$m<$binNumY} {incr m 1} {
			puts $hbAll "[expr {$xl+$dx*$k}]	[expr {$yl+$dy*$m}]	0.0"	
		}
	}			
# ---------------------------------------- HB on X axis ----------------------------------------
} elseif {[string equal -nocase $projection x]} {	
	for {set k 0} {$k<$binNumX} {incr k 1} {
		puts $hbAll "[expr {$xl+$dx*$k}]	0.0"	
	}
# ---------------------------------------- HB on Y axis ----------------------------------------		
} elseif {[string equal -nocase $projection y]} {
	for {set k 0} {$k<$binNumY} {incr k 1} {
		puts $hbAll "[expr {$yl+$dy*$k}]	0.0"	
	}	
# ---------------------------------------- HB on Z axis ----------------------------------------
} elseif {[string equal -nocase $projection z]} {
	for {set k 0} {$k<$binNumZ} {incr k 1} {
		puts $hbAll "[expr {$zl+$dz*$k}]	0.0"	
	}		
}
		
close $hbAll 

# refine the data in file hbAll.dat to file $filename
set inputFile [open "hbAll.dat" r]
set binnedData [dict create]

if {[string equal -nocase $projection xy] || [string equal -nocase $projection yz] || [string equal -nocase $projection xz]} {
	while {[gets $inputFile line] >= 0} {
		if {[string match "#*" $line]} {
			continue
		} 
		set coordinates [split $line]
		set pos_1 [lindex $coordinates 0]
		set pos_2 [lindex $coordinates 1]
		set value [lindex $coordinates 2]	
		
		set bin [list $pos_1 $pos_2]
		if {[dict exists $binnedData $bin]} {
			set existingvalue [dict get $binnedData $bin]
			dict set binnedData $bin [expr {$existingvalue + $value}]
		} else {
			dict set binnedData $bin $value
		}
	}
} elseif {[string equal -nocase $projection x] || [string equal -nocase $projection y] || [string equal -nocase $projection z]} {
	while {[gets $inputFile line] >= 0} {
		if {[string match "#*" $line]} {
			continue
		} 
		set coordinates [split $line]
		set pos_1 [lindex $coordinates 0]
		set value [lindex $coordinates 1]	
		set bin $pos_1
		if {[dict exists $binnedData $bin]} {
			set existingvalue [dict get $binnedData $bin]
			dict set binnedData $bin [expr {$existingvalue + $value}]
		} else {
			dict set binnedData $bin $value
		}
	}
}

close $inputFile

set filename [format "hb-$file_1-$projection-%d-%d-%d-%.1f-%.1f-mode_%d" $firstFrame $lastFrame $step $cutoff $angle $mode]
set outputFile [open "$filename.dat" w]

dict for {bin value} $binnedData {
	puts $outputFile "$bin [expr {$value/$nFrame/$binArea}]"
}

close $outputFile

}  