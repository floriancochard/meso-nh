all:
	cd 1_INPUT_MNH/ && get_pgd_files
	cd 1_INPUT_MNH/ && run_prep_mesonh_xyz
	cd 2_INPUT_TOY/ && run_prep_toy
	cd A_RUN_MNH_TOY/ && run_mesonh_xyz
	cd A_RUN_MNH_TOY/ && run_plot
clean:
	cd 1_INPUT_MNH/ && clean_prep_mesonh_xyz
	cd 2_INPUT_TOY/ && clean_prep_toy
	cd A_RUN_MNH_TOY/ && clean_mesonh_xyz
	cd A_RUN_MNH_TOY/ && clean_plot
