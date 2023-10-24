library(structToolbox)


D = iris_DatasetExperiment()
D$sample_meta$class = D$sample_meta$Species

M = mean_centre() + PLSDA(number_components=2,factor_name='class')
M = model_apply(M,D)



