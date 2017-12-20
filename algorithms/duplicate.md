---
permalink: /algorithms/duplicate.html
layout: default
title: Samtools - Duplicate Marking
highlighting: yes
---
## Duplicate Marking

Duplicates are defined as having multiple templates whose aligned 5'
coordinates match.  For a paired-end template this requires both 
primary reads to have matching 5' coordinates.  Coordinates are
based on the unclipped position of the read against the reference.
Reads also need to match orientation.  When a duplicate is detected,
the overall highest quality template is kept and all others have
the duplicate flag set.

For primary reads, this definition is the same as used in 
Picard (v2.10.3) and Biobambam2 (bamstreamingmarkduplicates v2.0.57).
None of these tools use supplementary data when determining a
duplicate. However unlike Picard, supplementary reads in a duplicate
template do not have their flags modified in Samtools by default.

![Duplicate example](../images/duplicate_example.png)
