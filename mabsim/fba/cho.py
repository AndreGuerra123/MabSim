import os, sys;
import cobra as cb;
from d3flux import flux_map as fm;

model = cobra.io.read_sbml_model('cho.xml')

map = fm(model, display_name_format=lambda x: str(x.id), figsize=(300,250),

      flux_dict={rxn.id: None for rxn in model.reactions})

print(map.data)

