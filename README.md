# funcnorm
Functional Normalization Toolbox: Inter-subject alignment of cortical anatomy using FMRI response and connectivity correlations

Intro:

The code was developed to generate warps from the cortical meshes and FMRI data for 10 subjects viewing a movie (Raiders of the Lost Ark) then tested on FMRI data from these same subjects performing other tasks. When warps are applied to novel data, they bring into alignment areas with similar patterns of activity as measured by reduced variability and average correlation. The methods and their generality are still a matter of active research. The software comes with no guarantees, and your mileage will vary.

Mechanics:

The alignment starts where FreeSurfer [http://freesurfer.net] and SUMA [http://afni.nimh.nih.gov/afni/suma] leave off:
- the inputs are spherical triangulated meshes (one per hemisphere per subject) aligned by FreeSurfer's 'recon-all' command using high-resolution anatomical scans from each subject, averaged across 3 or more sessions, resampled uniformly by SUMA's MapIcosahedron at 2mm resolution (~36k vertices per mesh) so that vertices are in 1:1 anatomic alignment across subjects.
- the algorithm then aligns subject pairs, by maximizing correlation while varying the position of each vertex, using a cost function that minimizes folding and stretching
- the warps are further constrained such that the group geometric average location of each vertex remains at its original (anatomically-aligned) location.
- the outputs are warps for each subject (both left and right hemisphere), which can be used to align future data sets.
- method one ("funcnorm", Sabuncu et al 2010) uses the inter-subject correlation of FMRI response directly to drive the alignment.
- method two ("connnorm", Conroy et al 2013) uses functional connectivity (the pattern of within-subject correlations) to drive the alignment.

Going forward "funcnorm" is meant as an umbrella term for both methods. The first version added to the repository (tagged 'v1.0' in git, and available via [https://github.com/aligners/funcnorm/archive/v1.0.zip]) corresponds to the state of the code at the time of the 2013 Conroy et al publication.

See also the project home page:

https://pni.princeton.edu/pni-software-tools/princeton-functional-normalization-toolbox

References:

- Conroy, BR, Singer BD, Guntupalli JS, Ramadge PJ and Haxby JV (2013). Inter-subject alignment of human cortical anatomy using functional connectivity. NeuroImage. 81, 400-411. [http://www.ncbi.nlm.nih.gov/pubmed/23685161]

- Sabuncu MR, Singer BD, Conroy BR, Ramadge PJ and Haxby JV (2010). Function-based inter-subject alignment of cortical anatomy. Cerebral Cortex, 20, 130-140. [http://www.ncbi.nlm.nih.gov/pubmed/19420007]
