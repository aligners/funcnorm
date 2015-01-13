# funcnorm
Functional Normalization Toolbox: Inter-subject alignment of cortical anatomy using FMRI response and connectivity correlations

Intro:

The code was developed to generate warps from the FMRI data of subjects viewing a movie (Raiders of the Lost Ark) and then used to generalize to fit FMRI data from these same subjects performing other tasks. When warps are applied to novel data, they bring into alignment areas with similar patterns of activity as measured by reduced variability and average correlation.

Mechanics:

The alignment starts where FreeSurfer [http://freesurfer.net] and SUMA [http://afni.nimh.nih.gov/afni/suma] leave off:
- the inputs are spherical triangulated meshes (one per hemisphere per subject) aligned by FreeSurfer based on a common atlas of cortial curvature landmarks, resampled uniformly by SUMA's MapIcosahedron at 2mm resolution (~36k vertices per mesh) so that vertices are in 1:1 alignment across subjects.
- the algorithms then align subject pairs, by maximizing correlation while manipulating the position of each vertex, using a cost function that minimizes folding and stretching
- The warps are further constrained such that the group geometric average location of each vertex remains at its original (anatomically-aligned) location.

- method one ("funcnorm", Sabuncu et al 2010) uses the inter-subject correlation of FMRI response directly to drive the alignment. Subject must view the same stimulus, 
- method two ("connnorm", Conroy et al 2013) uses functional connectivity (the pattern of within-subject correlations) to drive the alignment.

Going forward "funcnorm" is meant as an umbrella term for both methods. The same code base contains both methods but the code is undergoing validation testing to confirm this (January 2015).

More soon.

See also the project home page at:

https://pni.princeton.edu/pni-software-tools/princeton-functional-normalization-toolbox

References:

- Conroy, BR, Singer BD, Guntupalli JS, Ramadge PJ and Haxby JV (2013). Inter-subject alignment of human cortical anatomy using functional connectivity. NeuroImage. 81, 400-411. [http://www.ncbi.nlm.nih.gov/pubmed/23685161]

- Sabuncu MR, Singer BD, Conroy BR, Ramadge PJ and Haxby JV (2010). Function-based inter-subject alignment of cortical anatomy. Cerebral Cortex, 20, 130-140. [http://www.ncbi.nlm.nih.gov/pubmed/19420007]
