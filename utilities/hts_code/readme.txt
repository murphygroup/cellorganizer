|-----------------------------------------------------------------------|
| Hybrid Texture Synthesis MATLAB package                               |
|                                                                       |
| Date:   05/19/2003                                                    |
|                                                                       |
| Author: Andrew Nealen                                                 |
|         Discrete Geometric Modeling Group                             |
|         Technische Universitaet Darmstadt, Germany                    |
|                                                                       |
| Note:   This is part of the prototype MATLAB implementation           |
|         accompanying our paper/my thesis                              |
|                                                                       |
|         Hybrid Texture Synthesis. A. Nealen and M. Alexa              |
|         Eurographics Symposium on Rendering 2003                      |
|                                                                       |
|         Hybrid Texture Synthesis. A. Nealen                           |
|         Diplomarbeit (MSc), TU Darmstadt, 2003                        |
|                                                                       |
|         See the paper/thesis for further details.                     |
|-----------------------------------------------------------------------|



1. Introduction
---------------
What follows is a brief description on the intended usage of our 
prototype MATLAB implementation of Hybrid Texture Synthesis. Please 
note, that this is NON OPTIMIZED RESEARCH CODE! We have tried to keep 
all comments in the code consistent, yet there still my be the one or 
other misleading comment and/or badly named variable. All development
was carried out using MATLAB version 6.1.

Requirements:
-------------
- Matlab 6.1
- Image processing toolkit
- Optional: OpenTSTools Matlab package for use of build_knn.m. 
  available under the GPL from http://www.physik3.gwdg.de/tstool/.


2. Files
--------
This package contains the following .m files.

build_knn.m                   - builds an optional k-nearest 
                                neighbors datastructure for
                                k-coherence search
buildalphamask.m              - constructs the alpha-mask for boundary 
                                feathering
buildimagemask.m              - constructs the image mask used for 
                                best patch search
buildtraversalmap.m           - constructs the traversal map for 
                                overlap re-synthesis
circshiftpatch.m              - circular shift of a geometrically
                                defined patch
errorimage.m                  - constructs the errorimage for a 
                                given mask and input texture
errorimage_dst.m              - destination importance map weighted
                                error image
errorimage_sqdiff.m           - src/dst squared difference 
                                importance map weighted error image
errorimage_src.m              - source importance map weighted
                                error image
errorimage_sum.m              - src/dst sum importance map 
                                weighted error image
errorsurface.m                - computes the error surface
frequencymap.m                - a frequency map in [0.2,1]
gen_input.m                   - generates input to hybridsynthesize.m
genquadmeshoffset.m           - generates a 2D quad-mesh 
gentoroidalquadmeshoffset.m   - generates a toroidal, 2D quad-mesh
gray2rgb.m                    - simple image conversion
hybridsynthesize.m            - the wrapper function for 
                                hybridsynthesizerec.m
hybridsynthesizerec.m         - the main algorithm
overlap_resynthesis_exh.m     - exhaustive O(kNlogN) per-pixel
                                re-synthesis
overlap_resynthesis_knn.m     - k-coherence O(rkn) per-pixel re-synthesis

Copy all the files in the zip archive to a path where Matlab can find them.


3. Basic usage notes
--------------------
We generally start off by constructing a k-nearest neighbors (knn)
datastructure from an input texture, which, depending on the
texture size, chosen neighborhood size, and variation in the
texture, can take a while. As we can save this structure to a .mat
file for future use, it generally pays off. NOTE: you will need to
have the OpenTSTools MATLAB package installed to use
build_knn.m, available under the GPL from

http://www.physik3.gwdg.de/tstool/. 

Otherwise just perform exhaustive search (see further below).

On to the MATLAB command line:

>
> T = imread('my_cool_texture.tif');   % load a texture
> [knn,dist] = build_knn(T,3,10);      % 7x7 neighborhood, k = 10 (knn)
> save knn_for_cool_texture knn dist   % save the knn/dist information
>

All other necessary texture information is held in a 'data' structure,
which is generated using gen_input.m (type 'help gen_input' for more
on how to modify the switches of the algorithm).

>
> data = gen_input(T,256,256,0,1,'nowrap',0.02,0.04,32,32,'simple');
>

So we want to generate a 256x256 texture with no error tolerance, 
feathering enabled, no wrapping of the input texture, 0.02 pixel tolerance,
0.04 patch tolerance, an initial patch size of 32x32, and the simple
distance metric for overlapping patches. All other necessary values
are set to defaults. The structure generated by gen_input can and
should be modified to one's liking. See the comments in the .m files
and the paper/thesis for a more detailed description.

So, finally, we want (fast, k-coherence search) RESULTS:

>
> I = hybridsynthesize(data, 7, knn);
>

Or optionally, perform an exhaustive per-pixel search

>
> I = hybridsynthesize(data, 7);
>

In both cases we use an initial pixel overlap for patch placement of 7
pixels. Using a value somewhere around patch_cornerlength/6 seems 
resonable, but can strongly vary from texture to texure. Experimentation
is key.


4. Special cases
----------------
'Patch Based Sampling' (PBS) and 'Non-parametric Sampling' (NPS) are 
special cases of 'Hybrid Texture Synthesis' and can be mimicked with 
the following input data:

PBS (patch_delta_max = pixel_delta_max = 1.0)

>
> data = gen_input(T,256,256,0.01,1,'nowrap',1.0,1.0,32,32,'simple');
> I = hybridsynthesize(data, 6);  & uses a 6 pixel overlap (= 32/6)
>

NPS (no feathering, single pixel grid)

>
> data = gen_input(T,256,256,0.01,0,'nowrap',1,0,1,1,'simple');
> I = hybridsynthesize(data, 4);  & uses a 9x9 neighborhood
>


5. Settings
-----------

In many values in 'data' can and should be modified. by typing

>
> data = gen_input(...);
> data.VARIABLE = VALUE;
>

the VALUE of VARIABLE can be changed. Here the list of variables/meanings:

VARIABLE                 MEANING

rows / cols              pixels rows/columns in result.

errtol                   candidate-patch-set error tolerance in [0,1]
                         0 selects only best match, 1 selects randomly´.
                         example: setting patch_delta_max = 1 and errtol = 1
                         results in random patch placement with no splits.

feather                  binary value. set to 0/1 to disable/enable
                         overlap feathering

mode                     wrapping of input texture. possible
                         strings: 'wrap', 'nowrap',

pixel_delta_max          pixel error tolerance.

patch_delta_max}         patch error tolerance.

psr / psc                pixels rows/columns per patch in the initial grid

metric                   patch selection metric. possible
                         strings: 'simple', 'src', 'dst', 'sum', 'sqdiff'.

exh_overlap /            pixel strip around target pixel for exhaustive /
knn_overlap              k-coherence search during overlap re-synthesis. 
                         results in a nxn neighborhood with n=2*overlap+1.

k                        the 'k' in k-coherence search. defaults to 5.

min_patchsize            smallest allowed patch size in pixels. default is 8.

display_resynthesis      set to 0 to disable algorithm visualization.

verbatim                 set to 1 for some command line output.


6. Contact
----------
I will gladly answer any questions concerning the basic algorithm.
I am aware, that some of the loops in this code could be tightened/
vectorized, yet i didn't find the time to do so. It is just meant
as a prototype in any case, and MATLAB is great for this task! You
can reach me at

andy@nealen.com


7. Disclaimer
-------------
Repeating myself here... This is research/prototype code. We will take
no responsibility concerning its use by others and simply published it
to give better insight to the inner workings of the algorithm.


Other than that, have fun fiddling with it! and let me know
if you have any constructive comments/suggestions.


Cheers,
Andy

05/19/2003

__________________________________________________
Andrew Nealen
Discrete Geometric Modeling Group
Technische Universitaet Darmstadt, Germany
http://www.dgm.informatik.tu-darmstadt.de
http://www.nealen.com
andy@nealen.com


VERSION HISTORY

ver 1.00 = Released May 19th 2003

	 - Initial version


ver 1.01 = released September 16th 2003

         - Removed a nasty bug in overlap_resynthesis_exh.m, 
           which only allowed for input textures with identical
           width and height.
