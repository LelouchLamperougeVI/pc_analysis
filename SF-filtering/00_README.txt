README.txt

Two-dimensional spatial frequency filtering by FFT using Matlab

Authors: Kota S. Sasaki and Izumi Ohzawa
Graduate School of Frontier Biosciences, Osaka University
kota@fbs.osaka-u.ac.jp, ohzawa@fbs.osaka-u.ac.jp

License: BSD license, http://creativecommons.org/licenses/BSD/

Date: 2007-12-20

------------------------------------------

SFFiltering1script.m:

This Matlab script band-pass filters a 2-dimensional image using FFT.
Computed results are displayed in a window.  A captured version of this
is in results/SFFiltering1script.png.
Top left: the original image (gray scaled)
Top right: amplitude spectrum (logarithmically scaled amplitude)
Bottom left: 2-d filter (circularly symmetric Gaussian band-pass filter)
Bottom right: filtered image.

SFFiltering2script.m:

This Matlab script applies an oriented band-pass filter on a 2-dimensional image using FFT.
Computed results are displayed in a window.  A captured version of this
is in results/SFFiltering2script.png.
The filter shape is equivalent to that of a filter by a single simple cell in V1.
(Think of the results as an image encoded by a single type of simple cell
tuned to one spatial frequency band and one orientation band.)

Top left: the original image (gray scaled)
Top right: amplitude spectrum (logarithmically scaled amplitude)
Bottom left: 2-d filter (oriented Gaussian band-pass filter)
Bottom right: filtered image.


GetLuminanceImage.m:

Reads in an image (may be color) and converts it into a gray scale image
for computation.

Myff2.m:

My version of 2D FFT for better handling.
