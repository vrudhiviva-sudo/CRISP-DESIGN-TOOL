# CRISP-DESIGN-TOOL
CRISPR design tools are used to convert guide RNAs (gRNAs) into targeted DNA sequences for gene editing. These tools help in designing gRNAs that can specifically target and edit the desired genomic sites, enabling precise genetic modifications in living organisms. 
CRISPR-Cas design tools are computer software platforms and bioinformatics tools used to facilitate the design of guide RNAs (gRNAs) for use with the CRISPR/Cas gene editing system. System was originally discovered to be an acquired immune response mechanism used by archaea and bacteria. It has since been adopted for use as a tool in the genetic engineering of higher organisms.
From a raw DNA sequence to a ranked list of guides, off-target checks and ready-to-order cloning oligos, explained step by step. Designing an appropriate gRNA is an important element of genome editing with the CRISPR/Cas system. A gRNA can and at times does have unintended interactions ("off-targets") with other locations of the genome of interest. For a given candidate gRNA, these tools report its list of potential off-targets in the genome thereby allowing the designer to evaluate its suitability prior to embarking on any experiments.Scientists have also begun exploring the mechanics of the CRISPR/Cas system and what governs how good, or active, a gRNA is at directing the Cas nuclease to a specific location of the genome of interest.[1][2] As a result of this work, new methods of assessing a gRNA for its 'activity' have been published,[1][2] and it is now best practice to consider both the unintended interactions of a gRNA as well as the predicted activity of a gRNA at the design stage.

CRISPR in four ideas: The genome is a book
                    : The gRNA is the address
                    : Cas9 is the scissors
                    : The PAM is the doorbell

Six steps inside the app
    1.Find the PAM
    2.Filter the weak ones
    3.Score on-target
    4.Search for off-targets
    5.Score specificity
    6.Rank and export

Project Structure
          app.py – Web application interface.
          
