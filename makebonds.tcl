# for drawing the protein in VMD. 
# in VMD console, after opening the structure/trajectory, type "source <path to makebonds.tcl>", 
# then "drawbonds <15bonds.txt or 30bonds.txt", depending on which condition>
proc drawbonds { fi } {
  set fp [open $fi r]
  set file_data [read $fp]
  close $fp
  
  #  Process data file
  set data [split $file_data "\n"]
  foreach line $data {
    set l [split $line " "]

    set r1 [lindex $l 0]
    set n1 [lindex $l 1]
    set r2 [lindex $l 2]
    set n2 [lindex $l 3]
    

    set a [atomselect top "resid $r1 and name $n1"] 
    set aindex [$a get index]

    set a1 [atomselect top "resid $r2 and name $n2"] 
    set a1index [$a1 get index]

    topo addbond $aindex $a1index
  }
}