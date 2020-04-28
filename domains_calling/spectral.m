input_file = 'matrix.tab'
resolution = 1e4 % 10kb

thresholds = 0:0.01:1
for i = 1:length(thresholds)

	HiC_wCent = HiC_load_mat(input_file)
	[HiC,~,idx_cent] = HiC_remove_cent(HiC_wCent)

	TADs_iter = TAD_Laplace(HiC, thresholds(i))
	
	output_file = ['output_TADs_' num2str(thresholds(i)) '.txt']
	saveTADs(output_file, TADs_iter, resolution, idx_cent)
end
