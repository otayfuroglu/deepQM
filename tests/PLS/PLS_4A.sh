export GMX_MAXBACKUP=-1
workdir=$(pwd)

## we need to remove pbc and center protein in BOX
echo 1 0 | gmx trjconv -f md_0.xtc -s min_steep_0.tpr -n index.ndx -o trjmol.xtc -pbc mol -ur compact -center

#### we need to create ordered pdb file for indexing 24-PL + 4A 17-SOL
echo 24 17 | gmx trjorder -f trjmol.xtc -s md_0.tpr -n index.ndx -o ordered.pdb -da 0 -r 0.6 -dt 9000

#### we need to convert xtc file so that SOL is ordered around protein
echo 24 17 | gmx trjorder -f trjmol.xtc -s md_0.tpr -n index.ndx -o ordered.xtc -nshell nshell.xvg -da 0 -r 0.4 -dt 100

### we need to print number of SOL in nshell!
awk '{ total += $2; count++ } END { print int(total/count) }' nshell.xvg
natoms=$(awk '{ total += $2; count++ } END { print 3*int(total/count) }' nshell.xvg)
echo "$natoms"
nwater=$((natoms/3))

### we need to print number of atoms in PL complex
echo 24 | gmx editconf -f md_0.gro -o pl.gro -n index.ndx
filename='pl.gro'
pro_atm=$(cat ${filename} | wc -l)
pro_atm=$((pro_atm-3))
solstart=$((pro_atm+1))
total=$((pro_atm+natoms))
echo "There are " $total " atoms in the protein"

#--> take avg
gmx make_ndx -f ordered.pdb -o index.ndx -n index.ndx <<EOF
a1-$total
q
EOF

## we need to convert xtc file to new xtc file including only PL + 4A SOL from ordered trajectory
echo a_1-"$total" | gmx trjconv -f ordered.xtc -s min_steep_0.tpr -o pls.xtc -n index.ndx

## we need a reefrence structure including only PL + 4A SOL from ordered trajectory
echo a_1-"$total" | gmx trjconv -f ordered.xtc -s min_steep_0.tpr -o pls.gro -n index.ndx -dump 0

## below is to create index file including PL + 4A SOL as 0-System and 1|13 as 17-PL
gmx make_ndx -f pls.gro -o index_4A.ndx<<EOF
q
EOF

## we have to modify topoloji file fixing SOL number
head -n -3 "$workdir"/"$pro"/topol.top > topol_sol_number.top

echo "SOL    $nwater" >> topol_sol_number.top

mkdir pdb_pro_lig_sol

# we need to convert xtc file to new xtc file including only PL + 4A SOL from ordered trajectory
echo 0 | gmx trjconv -f pls.xtc -s min_steep_0.tpr -o pdb_pro_lig_sol/trjmol.pdb -n index_4A.ndx -sep
