add_library(fabm_models_deb OBJECT
            model_library.F90
            population.F90
           )

add_dependencies(fabm_models_deb fabm_base)

#-----------------------------------------------------------------------
# Copyright (C) 2015-2016 - PML
#-----------------------------------------------------------------------
