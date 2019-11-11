all: E001_pgd1 E002_arp2lfi E003_mesonh E004_diag E006_ncl 

E001_pgd1 :
	cd 001_pgd1	       && get_chimie_files
	cd 001_pgd1	       && run_prep_pgd_xyz
E002_arp2lfi :
	cd 002_arp2lfi	       && run_arp2lfi_xyz
E003_mesonh:
	cd 003_mesonh          && run_mesonh_xyz
E004_diag:
	cd 004_diag            && run_diag_xyz
E006_ncl:
	cd 006_ncl	       && run_ncl

clean:
	cd 001_pgd1            && clean_prep_pgd_xyz
	cd 002_arp2lfi	       && clean_arp2lfi_xyz
	cd 003_mesonh	       && clean_mesonh_xyz
	cd 004_diag  	       && clean_diag_xyz
	cd 006_ncl	       && clean_ncl

