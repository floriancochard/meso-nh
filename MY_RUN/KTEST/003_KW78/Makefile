all:
	cd 001_prep_ideal_case && run_prep_ideal_case_xyz
	cd 002_mesonh          && run_mesonh_xyz
	cd 003_convdia         && run_conv2dia
	cd 004_diaprog         && run_diaprog
	cd 005_ncl             && run_ncl
	cd 006_diag            && run_diag_xyz
	cd 007_ncl_diag        && run_ncl

clean:
	cd 001_prep_ideal_case && clean_prep_ideal_case_xyz
	cd 002_mesonh          && clean_mesonh_xyz
	cd 003_convdia         && clean_conv2dia
	cd 004_diaprog         && clean_diaprog
	cd 005_ncl             && clean_ncl
	cd 006_diag            && clean_diag_xyz
	cd 007_ncl_diag        && clean_ncl
