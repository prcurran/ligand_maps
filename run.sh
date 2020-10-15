#!/bin/bash
conda activate hotspots

#python ligand_overlay.py "P31749" "3.5" "3CQW" "A" -o "results/AKT1" -b "CQW999"
#python ligand_pharmacophore.py 'results/AKT1/ligand_overlay.mol2' '3CQW' -o 'results/AKT1'
#python ligand_overlay.py "P00811" "3.5" "1L2S" "A" -o "results/AMPC" -b "STC1115"
#python ligand_pharmacophore.py 'results/AMPC/ligand_overlay.mol2' '1L2S' -o 'results/AMPC'
#python ligand_overlay.py "Q00534" "3.5" "3NUX" "A" -o "results/CP3A4" -b "3NV900"
#python ligand_pharmacophore.py 'results/CP3A4/ligand_overlay.mol2' '3NUX' -o 'results/CP3A4'
#python ligand_overlay.py "P04150" "3.5" "3BQD" "A" -o "results/GCR" -b "DAY301"
#python ligand_pharmacophore.py 'results/GCR/ligand_overlay.mol2' '3BQD' -o 'results/GCR'
python ligand_overlay.py "P03367" "3.5" "1XL2" "A" -o "results/HIVPR" -b "1891001"
python ligand_pharmacophore.py 'results/HIVPR/ligand_overlay.mol2' '1XL2' -o 'results/HIVPR'
python ligand_overlay.py "P04585" "3.5" "3LAN" "A" -o "results/HIVRT" -b "KBT561"
python ligand_pharmacophore.py 'results/HIVRT/ligand_overlay.mol2' '3LAN' -o 'results/HIVRT'
#python ligand_overlay.py "P52732" "3.5" "3CJO" "A" -o "results/KIF11" -b "K301"
#python ligand_pharmacophore.py 'results/KIF11/ligand_overlay.mol2' '3CJO' -o 'results/KIF11'