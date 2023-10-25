Common problems with fitting
============================

Here, I have listed some common issues that can make the fit results worse. If
you, the reader, has other issues which you expect other people to run into,
feel free to contribute to this document!

Bad pixels in my maps
---------------------

If during the fitting, a message is printed out for the pixel in question ("Too
many bad values"), then the pixel is skipped by design. This is done when more
than 30% (number subject to change) of the data points are nan or exactly zero.
If the data are fine

This usually means that the fitting went wrong for those particular spectra.
Here is a list of things to check
- Zero or very small uncertainties can completely throw off the minimizer
- If the initial guess is too close, but in the wrong way, it can get stuck in a
  local minimum.
