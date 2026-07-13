# Joint typology

This folder contains the scripts used to code discrete joint characters,
assign joint-type categories and summarize observed versus theoretically
possible joint-character combinations.

Recommended order:

1. `assign_joint_type_categories.R`
   standardizes the joint-character table and assigns mechanical joint-type
   categories per specimen.
2. `audit_joint_character_combinations.R`
   enumerates theoretical character combinations and checks which combinations
   are observed or absent in the dataset.
3. `summarize_joint_character_combinations.R`
   creates summary tables and plots for the observed and missing
   joint-character combinations.
