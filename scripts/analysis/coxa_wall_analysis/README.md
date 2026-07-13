# Coxa-wall analyses

This folder contains the scripts used to check and analyse the coxal wall
opening character. The files are not alternative versions of the same analysis;
they represent different steps of the coxa-wall workflow.

Recommended order:

1. `update_coxal_wall_opening_coding.R`
   Standardizes the binary coxal wall opening character used in the downstream
   analyses.
2. `update_coxa_wall_character_coding.R`
   Updates the broader coxa-wall character table used for joint-character
   summaries.
3. `analyze_coxa_size_association.R`
   Tests whether coxal wall opening is associated with coxa size or coxal wall
   thickness.

The scripts assume the input tables distributed with the study and use
placeholder project roots rather than workstation-specific paths.
