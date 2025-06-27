#!/bin/bash
#SBATCH -N1 -t 0-60:00 -n 1 --mem 12000 --licenses sps --account nemo


source /sps/nemo/scratch/chauveau/software/gitlab-snfee-goliviero/this_snfee.sh

pdf_num=1519
ref_num=(1561 1566 1571 1576 1580 1584 1591 1596 1607 1612 1616 1620 1625 1631 1635 1646 1650 1656 1662 1666 1671 1675 1686 1691 1695 1699 1706 1711 1715 1727 1732 1736 1741 1748 1756 1760 1773 1777 1781 1786 1791 1796 1803 2015 2019 2023 2026 2031 2035 2040 2044 2057 2061 2066 2070 2074 2082 2088 2092 2111 2115 2121 2130 2135 2139 2156 2161 2165 2169 2175 2180 2185 2256 2300)
#ref_num=(1561 1566 1571 1576 1580 1584 1591 1596 1607 1612 1616 1620 1625 1631 1635 1646 1650 1656 1662 1666 1671 1675 1686 1691 1695 1699 1706 1711 1715 1727 1732 1736 1741 1748 1756 1760 1773 1777 1781 1786 1791 1796 1803 2015 2019 2023 2026 2031 2035 2040 2044 2057 2061 2066 2070 2074 2082 2088 2092 2111 2115)
#ref_num=(1561)
 
mkdir entree
mkdir entree/root/
mkdir entree/Modele/
mkdir sortie
mkdir sortie/root_file_final
mkdir sortie/root_file_final/alpha_fit


#To create the pdf
if [[ ! -f "entree/Modele/Modele_OM_${pdf_num}.root" ]]; then
	./charge_amplitude_energie_RHD -i /sps/nemo/snemo/snemo_data/raw_data/RHD/v3/snemo_run-${pdf_num}_crate-1_rhd.data.gz -o entree/Modele/Modele_OM_${pdf_num}.root -r $pdf_num
else
    echo "file entree/Modele/Modele_OM_${pdf_num}.root already exist"
fi

#Divide the pdf spectra analysis for the 5 ref OMs
if [[ ! -f "entree/Modele/Modele_OM_712_${pdf_num}.root" ]]; then
    ./Modele_OM ${pdf_num}
fi

 
#To create OM ref files
for num in ${ref_num[@]}; do
    if [[ ! -f "entree/root/histo_ref_${num}.root" ]]; then
	echo "./charge_amplitude_energie -i /sps/nemo/snemo/snemo_data/raw_data/RTD/v3/snemo_run-${num}_rtd.data.gz -o entree/root/histo_ref_${num}.root -r ${num}"
	./charge_amplitude_energie -i /sps/nemo/snemo/snemo_data/raw_data/RTD/v3/snemo_run-${num}_rtd.data.gz -o entree/root/histo_ref_${num}.root -r ${num}
    else echo "file histo_ref_${num}.root already exist"
    fi
done

./corr_om_ref_exe ${pdf_num} "${ref_num[@]}"

