all: E001_pgd1  E002_ecmwf2lfi E003_run E004_restart E005_conv2dia E006_diaprog E007_ncl

E001_pgd1 :
	cd 001_pgd1	       && run_prep_pgd_xyz 
E002_ecmwf2lfi :
	cd 002_ecmwf2lfi       && run_ecmwf2lfi_xyz
E003_run:
	cd 003_run	       && run_mesonh_xyz
E004_restart:
	cd 004_restart	       && run_mesonh_xyz
E005_conv2dia:
	cd 005_conv2dia	       && run_conv2dia
E006_diaprog:
	cd 006_diaprog	       && run_diaprog
E007_ncl:
	cd 007_ncl 	       && run_ncl
clean:
	cd 001_pgd1	       && clean_prep_pgd_xyz
	cd 002_ecmwf2lfi       && clean_ecmwf2lfi_xyz
	cd 003_run	       && clean_mesonh_xyz
	cd 004_restart	       && clean_mesonh_xyz		
	cd 005_conv2dia	       && clean_conv2dia	
	cd 006_diaprog	       && clean_diaprog
	cd 007_ncl	       && clean_ncl		
