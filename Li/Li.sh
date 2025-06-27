#!/bin/bash
#SBATCH -N1 -t 0-60:00 -n 1 --mem 15000 --licenses sps --account nemo

source /sps/nemo/scratch/chauveau/software/gitlab-snfee-goliviero/this_snfee.sh

ref_correction="1561-2115"
run_num=(1560 1565 1570 1575 1579 1583 1590 1595 1603 1606 1611 1615 1619 1624 1630 1634 1642 1645 1649 1655 1661 1665 1670 1674 1682 1685 1690 1694 1698 1705 1710 1714 1723 1726 1731 1735 1740 1747 1755 1759 1772 1776 1780 1785 1790 1795 1802 2014 2018 2022 2025 2030 2034 2038 2043 2060 2065 2069 2073 2081 2087 2091 2107 2110 2114 2120 2124 2129 2134 2138 2152 2155)
#run_num=(1595 1603 1606 1611 1615 1619 1624 1630 1634 1642 1645 1649 1655 1661 1665 1670 1674 1682 1685 1690 1694 1698 1705 1710 1714 1723 1726 1731 1735 1740 1747 1755 1759 1772 1776 1780 1785)
mkdir entree
mkdir entree/root/
mkdir sortie
mkdir sortie/SN_Li
mkdir sortie/ref_Li
mkdir sortie/ref_Li/Fit_Ampl_Ref
mkdir sortie/ref_Li/fit_Li

#Create Li files
for num in "${run_num[@]}"; do
    if [[ ! -f "entree/root/snemo_run-${num}_LI.root" ]]; then
	./Li_system -i /sps/nemo/snemo/snemo_data/raw_data/RTD/v3/snemo_run-${num}_rtd.data.gz -o entree/root/snemo_run-${num}_LI.root -r $num
    else echo "file snemo_run-${num}_LI.root already exist"
    fi
done

#Analyse Li runs
./Li "$ref_correction" "${run_num[@]}"
