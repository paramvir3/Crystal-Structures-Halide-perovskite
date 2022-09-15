
#!/usr/bin/tclsh
if {[catch {package require topotools 1.1} ver]} {
  vmdcon -error "$ver. This script requires at least TopoTools v1.1. Exiting..."
    quit
}
if {[catch {package require pbctools 2.3} ver]} {
  vmdcon -error "$ver. This script requires at least pbctools v2.3. Exiting..."
    quit
}
set fname slab.pdb
# check for presence of coordinate file
if {! [file exists $fname]} {
  vmdcon -error "Required file '$fname' not available. Exiting..."
    quit
}

mol new $fname autobonds no waitfor all

proc distance {new_pos x y z} {
  set x1 [lindex $new_pos 0]
    set bc [pbc get]
    set lx [lindex $bc 0 0]
    set ly [lindex $bc 0 1]
    set lz [lindex $bc 0 2]
    set y1 [lindex $new_pos 1]
    set z1 [lindex $new_pos 2]
    set dx [expr $x - $x1]
    set dy [expr $y - $y1]
    set dz [expr $z - $z1]
    if {$dx >   $lx * 0.5} { set dx [expr $dx - $lx] }
  if {$dx <= -$lx * 0.5} { set dx [expr $dx + $lx] }
  if {$dy >   $ly * 0.5} { set dy [expr $dy - $ly] }
  if {$dy <= -$ly * 0.5} { set dy [expr $dy + $ly] }
  if {$dz >   $lz * 0.5} { set dz [expr $dz - $lz] }
  if {$dz <= -$lz * 0.5} { set dz [expr $dz + $lz] }
  set dr [expr ($dx)**2+($dy)**2+($dz)**2]
    return $dr
}
############### SOLUTE ##################
set selH [atomselect top {name H}]
$selH set type H
$selH set mass 1.00800
$selH set charge 0.4648
set numH [$selH num]
set posH [$selH get {x y z}]
set indexH [$selH get index]

set selHe [atomselect top {name He}]
$selHe set type He
$selHe set mass 1.00800
$selHe set charge 0.22250000 
set numHe [$selHe num]
set posHe [$selHe get {x y z}]
set indexHe [$selHe get index]

set selC [atomselect top {name C}]
$selC set type C
$selC set mass 12.0100
$selC set charge 0.56710000 
set numC [$selC num]
set posC [$selC get {x y z}]
set indexC [$selC get index]

set selN [atomselect top {name N}]
$selN set type N
$selN set mass 14.0100
$selN set charge -0.82450000 
set numN [$selN num]
set posN [$selN get {x y z}]
set indexN [$selN get index]

set selPb [atomselect top {name Pb}]
$selPb set type Pb
$selPb set mass 207.2
$selPb set charge 2.00
set numPb [$selPb num]
set posPb [$selPb get {x y z}]
set indexPb [$selPb get index]

set selI [atomselect top {name I}]
$selI set type I
$selI set mass 126.90
$selI set charge -1.00
set numI [$selI num]
set posI [$selI get {x y z}]
set indexI [$selI get index]
#######################################
set sel [atomselect top all]
set natoms [molinfo top get numatoms]
set bondDistanceHN 1.2
set bondDistanceNC 2.5
set bondDistanceCHa 1.5
set selText all
set selTemp [atomselect top $selText]
set pos [$selTemp get {x y z}]
set index [$selTemp get index]
$selTemp delete
set bondDistance2HN [expr $bondDistanceHN*$bondDistanceHN]
set bondDistance2NC [expr $bondDistanceNC*$bondDistanceNC]
set bondDistance2CHa [expr $bondDistanceCHa*$bondDistanceCHa]
set blist {}

############ SOLUTE ##################
foreach r $posC ind1 $indexC {
# Select neighboring atoms.
  foreach {x y z} $r {break}
  foreach new_pos $posHe ind2 $indexHe {
    set dr [ distance $new_pos $x $y $z ]
      if {($dr < $bondDistance2CHa)} { lappend blist [list $ind1 $ind2]}
  }
}
foreach r $posN ind1 $indexN {
# Select neighboring atoms.
  foreach {x y z} $r {break}
  foreach new_pos $posH ind2 $indexH {
    set dr [ distance $new_pos $x $y $z ]
      if {($dr < $bondDistance2HN)} { lappend blist [list $ind1 $ind2]}
  }
}
foreach r $posN ind1 $indexN {
# Select neighboring atoms.
  foreach {x y z} $r {break}
  foreach new_pos $posC ind2 $indexC {
    set dr [ distance $new_pos $x $y $z ]
      if {($dr < $bondDistance2NC)} { lappend blist [list $ind1 $ind2]}
  }
}

return $blist
set reD {}
foreach b $blist {
  set bPerm [list [lindex $b 1] [lindex $b 0]]
    set match [lsearch $reD $bPerm]
    if {$match == -1} {lappend reD $b}
}
return $reD
topo setbondlist type $reD
topo retypebonds 
vmdcon -info "assigned [topo numbondtypes] bond types to [topo numbonds] bonds:"
vmdcon -info "bondtypes: [topo bondtypenames]"
topo guessangles  
vmdcon -info "assigned [topo numangletypes] angle types to [topo numangles] angles:"
vmdcon -info "angletypes: [topo angletypenames]"
topo guessdihedrals
vmdcon -info "assigned [topo numdihedraltypes] dihedral types to [topo numdihedrals] dihedrals:"
vmdcon -info "dihedraltypes: [topo dihedraltypenames]"
topo guessimpropers
vmdcon -info "assigned [topo numimpropertypes] improper types to [topo numimpropers] impropers:"
vmdcon -info "impropertypes: [topo impropertypenames]"
mol reanalyze top
pbc get
# we use a high-level tool from to multiply the system.
#TopoTools::replicatemol top 4 4 4

# and write out the result as a lammps data file.
topo writelammpsdata data.FAPI full
#topo writegmxtop topol.top

# for easier testing and visualization, we
# also write out copies in .pdb and .psf format.
#animate write pdb solvated.pdb
#animate write psf solvated.psf
# final sanity check the whole system has to be neutral.
$sel delete
set sel [atomselect top all]
set totq [vecsum [join [$sel get charge] { }]]
if {[expr {abs($totq)}] > 0.0005} {
  vmdcon -warning "Total system not neutral: $totq"
}

# done. now exit vmd
quit

