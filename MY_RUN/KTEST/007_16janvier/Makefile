all: E001_pgd1 E002_pgd2 E003_nest E004_arp2lfi E005_spa_mod1_mod2 E006_preal E008_run2 E009_diag E010_conv2dia E011_diaprog  E011_ncl E012_spectre

all_ncl_nc4: E001_pgd1 E002_pgd2 E003_nest E004_arp2lfi E005_spa_mod1_mod2 E006_preal E008_run2 E009_diag E011_ncl_nc4

all_ncl_lfi2cdf: E001_pgd1 E002_pgd2 E003_nest E004_arp2lfi E005_spa_mod1_mod2 E006_preal E008_run2 E009_diag E010_conv2dia E011_ncl_lfi2cdf

all_ncl_extractdia: E001_pgd1 E002_pgd2 E003_nest E004_arp2lfi E005_spa_mod1_mod2 E006_preal E008_run2 E009_diag E010_conv2dia E011_ncl_extractdia

E001_pgd1 :
	cd 001_pgd1	       && get_pgd_files 
	cd 001_pgd1	       && run_prep_pgd_xyz 
E002_pgd2:
	cd 002_pgd2	       && run_prep_pgd_xyz 
E003_nest:
	cd 003_nest            && run_prep_nest_pgd_xyz
E004_arp2lfi :
	cd 004_arp2lfi	       && run_arp2lfi_xyz
E005_spa_mod1_mod2:
	cd 005_spa_mod1_mod2   && run_spawning_xyz
E006_preal:
	cd 006_preal	       && run_preal_xyz
E008_run2:
	cd 008_run2	       && run_mesonh_xyz
E009_diag:
	cd 009_diag	       && run_diag_xyz
E010_conv2dia:
	cd 010_conv2dia	       && run_conv2dia
E011_diaprog:
	cd 011_diaprog	       && run_diaprog
E011_ncl:
	cd 011_ncl             && run_ncl
E011_ncl_nc4:
	cd 011_ncl_nc4         && run_ncl
E011_ncl_lfi2cdf:
	cd 011_ncl_lfi2cdf     && run_ncl
E011_ncl_extractdia:
	cd 011_ncl_extractdia  && run_ncl
E012_spectre:
	cd 012_spectre	       && run_spectre_xyz

clean:
	cd 001_pgd1	       && clean_prep_pgd_xyz
	cd 002_pgd2	       && clean_prep_pgd_xyz
	cd 003_nest            && clean_prep_nest_pgd_xyz
	cd 004_arp2lfi	       && clean_arp2lfi_xyz
	cd 005_spa_mod1_mod2   && clean_spawning_xyz
	cd 006_preal	       && clean_preal_xyz
	cd 008_run2	       && clean_mesonh_xyz
	cd 009_diag	       && clean_diag_xyz
	cd 010_conv2dia	       && clean_conv2dia
	cd 011_diaprog	       && clean_diaprog
	cd 011_ncl_lfi2cdf     && clean_ncl
	cd 011_ncl_extractdia  && clean_ncl
	cd 011_ncl_nc4         && clean_ncl
	cd 011_ncl             && clean_ncl
	cd 012_spectre	       && clean_spectre_xyz
