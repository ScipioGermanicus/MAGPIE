![Version](https://img.shields.io/badge/version-0.1.1-blue)

# MAG-based PICRUSt2 Integrated Extension (MAGPIE) (work in progress)

MAGPIE is an integrated extension for PICRUSt2 that enables creating bacterial and archaeal default files from MAGs, effectively building a custom reference database for PICRUSt2. This enables building databases that are more specialised for certain environments and can improve the accuracy of PICRUSt2-based predictions.

MAGPIE is completely based on and faithful to the PICRUSt2 developers' instructions on how to recreate the necessary reference files

# Work in progress
This tool is not currently finished. It can create the required reference files for running PICRUSt2 based on a defined set of MAGs, however, the tool is not yet fully automated.

Note that you will need to obtain CheckM quality statistics, GTDB-Tk classifications, ssu models, and EggNOG annotations beforehand, as these steps are not yet integrated into MAGPIE. It is therefore necessary to point MAGPIE to existing reference files when running it.
You also need to add the resulting output files manually to an existing PICRUSt2 conda environment, replacing the old reference files.

Feel free to message me with questions or requests for further developing this tool.
